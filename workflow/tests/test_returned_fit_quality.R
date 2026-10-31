suppressPackageStartupMessages(library(data.table))
source('workflow/qc/returned_fit_quality.R')
d<-c(.1,1,10);p<-.2+(.9-.2)/(1+d)
fit<-data.table(sampleid='s',Hill_slope=1,E0=.9,E_inf=.2,EC50=1,Rsquare=1,n_doses=3L,fit_status='ok')
x<-data.table(sampleid='s',treatment1dose=d,viability=p)
q<-returned_curve_quality(x,fit,'sampleid','test')
stopifnot(abs(q$returned_variance_Rsquare-1)<1e-12,abs(q$returned_SSE_Rsquare-1)<1e-12)
# The package's variance-based statistic ignores a constant prediction bias.
x[,viability:=p-.1];q<-returned_curve_quality(x,fit,'sampleid','test')
stopifnot(abs(q$returned_variance_Rsquare-1)<1e-12,q$returned_SSE_Rsquare<1,abs(q$mean_bias_observed_minus_prediction+.1)<1e-12)
# Stored clipped asymptotes can disagree with the optimizer's reported statistic.
fit[,E0:=1];x[,treatment1dose:=c(.3,1,10)];x[,viability:=.2+(1.2-.2)/(1+treatment1dose)]
q<-returned_curve_quality(x,fit,'sampleid','test')
stopifnot(q$asymptote_at_boundary,abs(q$difference_reported_minus_returned)>1e-3)
x[,viability:=.5];q<-returned_curve_quality(x,fit,'sampleid','test')
stopifnot(is.na(q$returned_variance_Rsquare),is.na(q$returned_SSE_Rsquare))
bad<-copy(fit);bad[,n_doses:=4L]
stopifnot(inherits(try(returned_curve_quality(x,bad,'sampleid','test'),silent=TRUE),'try-error'))
cat('PASS: exact curve, bias-sensitive SSE diagnostic, clipped-parameter discrepancy, flat-response missingness and dose-count mismatch.\n')
