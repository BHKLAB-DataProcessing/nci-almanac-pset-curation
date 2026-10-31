# Evaluate stored parameters at the actual modeling inputs; never fit or mutate a PSet.
returned_curve_quality <- function(input, fits, keys, assay) {
  values <- c('Hill_slope','E0','E_inf','EC50','Rsquare','n_doses','fit_status')
  stopifnot(all(c(keys,values)%in%names(fits)),!anyDuplicated(fits,by=keys))
  fit <- copy(fits[,c(keys,values),with=FALSE])
  setnames(fit,c('Rsquare','n_doses'),c('reported_Rsquare','stored_n_doses'))
  input <- input[is.finite(treatment1dose)&treatment1dose>0&is.finite(viability)]
  before <- nrow(input)
  z <- merge(input,fit,by=keys,all.x=TRUE,sort=FALSE)
  stopifnot(nrow(z)==before,!anyNA(z$fit_status))
  z[,model_response:=pmin(1,pmax(0,viability))]
  z[,prediction:=NA_real_]
  z[is.finite(Hill_slope)&is.finite(E0)&is.finite(E_inf)&is.finite(EC50)&EC50>0,
    prediction:=E_inf+(E0-E_inf)*plogis(-Hill_slope*(log(treatment1dose)-log(EC50)))]
  q <- z[,.(model_points=.N,input_n_doses=uniqueN(treatment1dose),
    finite_predictions=sum(is.finite(prediction)),response_variance=var(model_response),
    residual_variance=var(model_response-prediction),
    SSE=sum((model_response-prediction)^2),
    mean_bias_observed_minus_prediction=mean(model_response-prediction)),by=keys]
  out <- merge(fit,q,by=keys,all.x=TRUE,sort=FALSE)
  stopifnot(nrow(out)==nrow(fit),
    !any(is.na(out$input_n_doses)&out$stored_n_doses>0),
    all(out$input_n_doses==out$stored_n_doses,na.rm=TRUE))
  out[,`:=`(returned_variance_Rsquare=NA_real_,returned_SSE_Rsquare=NA_real_)]
  out[model_points>1&finite_predictions==model_points&is.finite(response_variance)&response_variance>0,
    `:=`(returned_variance_Rsquare=1-residual_variance/response_variance,
      returned_SSE_Rsquare=1-SSE/((model_points-1)*response_variance))]
  out[,`:=`(assay=assay,asymptote_at_boundary=E0%in%c(0,1)|E_inf%in%c(0,1),
    difference_reported_minus_returned=reported_Rsquare-returned_variance_Rsquare)]
  out
}

collect_returned_fit_quality <- function(raw,tre,dataset) {
  mean_finite <- function(x)if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
  keys <- c('sampleid','tech_rep','plateid','treatment1id')
  if(dataset=='NCI_ALMANAC')keys<-c(setdiff(keys,'plateid'),'NSC1')
  if(dataset=='GDSC_MATRIX')keys<-c(keys,'matrixid')
  dosekeys <- c(keys,'treatment1dose')
  if(dataset=='GDSC_MATRIX') {
    left<-raw[,c(dosekeys,'viability_1'),with=FALSE];setnames(left,'viability_1','viability')
    rightkeys<-sub('treatment1','treatment2',dosekeys,fixed=TRUE)
    right<-raw[,c(rightkeys,'viability_2'),with=FALSE]
    setnames(right,c(rightkeys,'viability_2'),c(dosekeys,'viability'))
    mono<-unique(rbind(left,right));rm(left,right)
    mono<-mono[,.(viability=mean_finite(viability)),by=dosekeys]
    wide<-tre[['monotherapy_profiles']]
    values<-c('Hill_slope','E0','E_inf','EC50','Rsquare','n_doses','fit_status')
    fits<-unique(rbindlist(lapply(1:2,function(side){
      from<-c(sub('treatment1id',paste0('treatment',side,'id'),keys,fixed=TRUE),paste0(values,'_',side))
      z<-wide[,..from];setnames(z,from,c(keys,values));z
    })))
  } else {
    selected<-is.na(raw$treatment2id)
    if(dataset=='GDSC_ANCHOR')selected<-raw$source_is_single==1&raw$source_annotation_complete==1
    mono<-raw[which(selected),.(viability=mean_finite(viability)),by=dosekeys]
    fits<-tre[[if(dataset=='HMCL')'profiles' else 'monotherapy_profiles']]
  }
  ans<-returned_curve_quality(mono,fits,keys,'monotherapy')
  rm(mono,fits);invisible(gc())
  if(dataset=='GDSC_ANCHOR') {
    combo_keys<-c('sampleid','tech_rep','plateid','treatment1id','treatment2id','treatment2dose','anchor_doses')
    combo<-raw[!is.na(treatment2id)&source_annotation_complete==1,
      c(combo_keys,'treatment1dose','viability'),with=FALSE]
    ans<-rbind(ans,returned_curve_quality(combo,tre[['profiles']],combo_keys,'combination'),fill=TRUE)
  }
  ans
}
