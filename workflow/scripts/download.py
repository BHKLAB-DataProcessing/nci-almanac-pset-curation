"""Fetch one configured public source, with an atomic payload and SHA-256 receipt."""
import argparse
import csv
import email.message
import hashlib
import gzip
import json
import os
from pathlib import Path
import shutil
import time
import urllib.parse
import urllib.request
import yaml


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_content(path, source):
    expected = source.get("content_sha256")
    if expected:
        # Sanger's generated gzip header includes a changing timestamp.
        # Pin the scientific CSV bytes, while receipts still hash the whole download.
        digest = hashlib.sha256()
        with gzip.open(path, "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError("Decompressed source checksum differs from configured release: " + source["id"])


def download(source, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    payload, receipt = directory / "payload", directory / "receipt.json"
    expected_hash = source.get("sha256")
    if receipt.exists() and payload.exists():
        old = json.loads(receipt.read_text())
        if old["download_url"] == source["download_url"] and sha256(payload) == old["sha256"]:
            if expected_hash and old["sha256"] != expected_hash:
                raise ValueError("Cached source checksum differs from configured release: " + source["id"])
            verify_content(payload, source)
            return old
    url = source["download_url"]
    if urllib.parse.urlparse(url).scheme != "https":
        raise ValueError("Sources must use HTTPS: " + url)
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "BHKLab-PSet-curation/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response:
                disposition = email.message.Message()
                disposition["content-disposition"] = response.headers.get("Content-Disposition", "")
                filename = Path(disposition.get_filename() or urllib.parse.unquote(urllib.parse.urlparse(response.url).path).split("/")[-1] or source["id"]).name
                if "text/html" in response.headers.get("Content-Type", "").lower():
                    raise ValueError("Received an HTML page instead of source data: " + url)
                temporary = payload.with_suffix(".part")
                with temporary.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
                expected = response.headers.get("Content-Length")
                if expected and temporary.stat().st_size != int(expected):
                    raise IOError("Incomplete download: " + url)
                result = {**source, "filename": filename, "effective_url": response.url,
                          "bytes": temporary.stat().st_size, "sha256": sha256(temporary),
                          "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
            if expected_hash and result["sha256"] != expected_hash:
                temporary.unlink()
                raise ValueError("Downloaded source checksum differs from configured release: " + source["id"])
            verify_content(temporary, source)
            temporary.replace(payload)
            receipt_temporary = receipt.with_suffix(".part")
            receipt_temporary.write_text(json.dumps(result, indent=2) + "\n")
            receipt_temporary.replace(receipt)
            return result
        except ValueError:
            raise
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_id")
    args = parser.parse_args()
    config = yaml.safe_load(Path("config/pipeline.yaml").read_text())
    matching = [s for s in config["sources"] if s["id"] == args.source_id]
    if len(matching) != 1:
        raise ValueError("Expected exactly one configured source: " + args.source_id)
    print(json.dumps(download(matching[0], Path("data/rawdata/downloads") / args.source_id)))
