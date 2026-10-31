# Read-only scientific QC helpers. No fitting, acquisition or PSet mutation.
setDTthreads(1L)
qc_dir <- 'docs/qc-data'
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)
qc_table <- function(x, name, show=TRUE, n=Inf) {
  x <- as.data.table(x); fwrite(x, file.path(qc_dir, paste0(name,'.tsv')), sep='\t')
  if(show) cat('<div class="qc-table">',knitr::kable(head(x,n), digits=4, format='html', table.attr='class="table table-striped table-condensed"'),'</div>\n')
  invisible(x)
}
qc_plot <- function(p, name, width=9, height=5) {
  ggsave(file.path(qc_dir,paste0(name,'.png')),p,width=width,height=height,dpi=140,bg='white')
  print(p); invisible(p)
}
finite_cor <- function(x,y,method='pearson') {
  k<-is.finite(x)&is.finite(y)
  if(sum(k)<3L || sd(x[k])==0 || sd(y[k])==0) return(NA_real_)
  cor(x[k],y[k],method=method)
}
paired_stats <- function(x,y,label) {
  k<-is.finite(x)&is.finite(y); d<-y[k]-x[k]
  data.table(comparison=label,contexts=length(x),paired=sum(k),reference_only=sum(is.finite(x)&!is.finite(y)),
    recomputed_only=sum(!is.finite(x)&is.finite(y)),neither=sum(!is.finite(x)&!is.finite(y)),
    Pearson=finite_cor(x,y),Spearman=finite_cor(x,y,'spearman'),
    mean_bias=if(length(d))mean(d) else NA_real_,median_bias=if(length(d))median(d) else NA_real_,
    MAE=if(length(d))mean(abs(d)) else NA_real_,RMSE=if(length(d))sqrt(mean(d*d)) else NA_real_,
    p95_abs_error=if(length(d))unname(quantile(abs(d),.95)) else NA_real_,
    descriptive_lower=if(length(d)>1)mean(d)-1.96*sd(d) else NA_real_,
    descriptive_upper=if(length(d)>1)mean(d)+1.96*sd(d) else NA_real_)
}
plot_subset <- function(x,n=40000L) {if(nrow(x)<=n)x else {set.seed(20261027);x[sample.int(nrow(x),n)]}}
agreement_plots <- function(x,y,label,stem,units='Response fraction') {
  d<-data.table(reference=x,recomputed=y)[is.finite(reference)&is.finite(recomputed)]
  if(!nrow(d))return(invisible(NULL)); d<-plot_subset(d)
  qc_plot(ggplot(d,aes(reference,recomputed))+geom_point(alpha=.15,size=.5,color='#276983')+
    geom_abline(slope=1,intercept=0,linetype=2,color='#ad4c35')+coord_equal()+
    labs(title=label,subtitle='Identity line; statistics use all finite pairs, points capped at 40,000',x=paste('Reference:',units),y=paste('Recomputed:',units)),paste0(stem,'-identity'))
  s<-paired_stats(x,y,label)
  qc_plot(ggplot(d,aes((reference+recomputed)/2,recomputed-reference))+geom_point(alpha=.15,size=.5,color='#276983')+
    geom_hline(yintercept=0,linetype=2)+geom_hline(yintercept=s$mean_bias,color='#ad4c35')+
    geom_hline(yintercept=c(s$descriptive_lower,s$descriptive_upper),linetype=3,color='#ad4c35')+
    labs(title=paste(label,'— difference versus mean'),subtitle='Red: mean bias and descriptive ±1.96 SD; repeated contexts are not independent replicates',x=paste('Pair mean:',units),y='Recomputed − reference'),paste0(stem,'-difference'))
}
qc_distribution <- function(x,cols,stem,labels=NULL) {
  z<-rbindlist(lapply(cols,function(n)data.table(metric=if(is.null(labels))n else labels[[n]],value=x[[n]])))
  qc_table(z[,.(rows=.N,finite=sum(is.finite(value)),missing=sum(!is.finite(value)),
    median=if(any(is.finite(value)))median(value[is.finite(value)]) else NA_real_,
    p05=if(any(is.finite(value)))quantile(value[is.finite(value)],.05) else NA_real_,
    p95=if(any(is.finite(value)))quantile(value[is.finite(value)],.95) else NA_real_),by=metric],paste0(stem,'-coverage'))
  z<-plot_subset(z[is.finite(value)],100000L)
  if(nrow(z))qc_plot(ggplot(z,aes(value))+geom_histogram(bins=60,fill='#276983',color='white')+facet_wrap(~metric,scales='free')+
    labs(x='Metric value',y='Sampled contexts',subtitle='Plot capped at 100,000 contexts; table uses all rows'),stem,height=6)
}
