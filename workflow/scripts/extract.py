"""Safely extract archives and prepare a source-member index for R stages."""
import csv
import argparse
import gzip
import json
from pathlib import Path, PurePosixPath
import shutil
import zipfile
import yaml
from download import sha256


def extract(source_id):
    raw = Path("data/rawdata/downloads") / source_id
    receipt = json.loads((raw / "receipt.json").read_text())
    payload = raw / "payload"
    if sha256(payload) != receipt["sha256"]:
        raise ValueError("Source checksum changed: " + source_id)
    destination = Path("data/procdata/extracted") / source_id
    destination.mkdir(parents=True, exist_ok=True)
    paths = []
    if zipfile.is_zipfile(payload):
        with zipfile.ZipFile(payload) as archive:
            for item in archive.infolist():
                member = PurePosixPath(item.filename)
                if member.is_absolute() or ".." in member.parts or (item.external_attr >> 16) & 0o170000 == 0o120000:
                    raise ValueError("Unsafe ZIP member: " + item.filename)
                # Figures and manuscripts are preserved in the source archive.
                if item.is_dir() or member.suffix.lower() not in (".csv", ".tsv", ".txt", ".xls", ".xlsx", ".yml", ".json"):
                    continue
                if source_id == "rnaseq" and member.name.startswith("rnaseq_all_data_"):
                    continue  # The same measurements are in the three assay matrices.
                target = destination.joinpath(*member.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(item) as src, target.open("wb") as out:
                    shutil.copyfileobj(src, out, length=1024 * 1024)
                paths.append(target)
    else:
        with payload.open("rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
        filename = receipt["filename"]
        if compressed and filename.endswith(".gz"):
            filename = filename[:-3]
        genes_only = filename.endswith(".gtf")
        if genes_only:
            filename = filename[:-4] + ".genes.gtf"
        target = destination / filename
        opener = gzip.open if compressed else open
        with opener(payload, "rb") as src, target.open("wb") as out:
            if genes_only:
                # Gene-level curation needs only gene records. Keep the complete
                # compressed source unchanged, avoiding multi-GB transcript copies.
                for line in src:
                    fields = line.split(b"\t", 3)
                    if len(fields) > 2 and fields[2] == b"gene":
                        out.write(line)
            else:
                shutil.copyfileobj(src, out, length=1024 * 1024)
        paths.append(target)
    return [{"source_id": source_id, "member": str(p.relative_to(destination)), "path": str(p)} for p in paths]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_id")
    args = parser.parse_args()
    config = yaml.safe_load(Path("config/pipeline.yaml").read_text())
    if args.source_id not in {s["id"] for s in config["sources"]}:
        raise ValueError("Unknown source id")
    rows = extract(args.source_id)
    destination = Path("data/procdata/extracted") / args.source_id
    with (destination / "manifest.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_id", "member", "path"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
