rm(list=ls())
set.seed(121)

library(rstan)
#library(parallel)
#library(doParallel)
#library(doRNG)
#source("~/Sim/RFunctionsMultiTreatMVBPR.R")
source("~/Sim/FunctionsMVBPRtoinclude.r")



library(rstan)



STANcode="
functions {
  // for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + lgamma(sum(y)+1) + sum(lgamma(alpha + to_vector(y)))
      - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha))-sum(lgamma(to_vector(y)+1));
    }
  }

data {
  int<lower=1> N; // total number of observations
  int<lower=1> M; // total number of observations in the prediction set
  int<lower=2> J; // number of categories of y
  int<lower=2> P; // number of predictor levels
  matrix[N,P] X; // predictor design matrix
  matrix[M,P] Xp; // predictor design matrix
  int <lower=0> Y[N,J]; // data // response variable
  //real<lower=0> sd_prior;
  real<lower=0> sd_prior;
  real<lower=0> psi;
}

parameters {
  matrix[P, J] beta_raw; // coefficients (raw)
  vector[J] beta0; // intercept
  matrix<lower=0>[P,J] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[J] tau; // global shrinkage
}

transformed parameters{
  matrix[P,J] beta; // coefficients
  matrix<lower=0>[P,J] lambda; // local shrinkage
  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda;
}

model {
  // prior:
  //for(k in 1:J){
  //  beta0[k] ~ normal(0, 10);
  //}
  beta0 ~ normal(0, 10);

  for (k in 1:P) {
    //for (l in 1:J) {
      //tau[l] ~ cauchy(0.1, 1); // flexible
      lambda_tilde[k,] ~ cauchy(0, 1);
      beta_raw[k,] ~ normal(0,sd_prior);
    //}
  }
  tau ~ cauchy(0.1, 1);

  for (i in 1:N) {
    vector[J] logits;
    for (j in 1:J){
      logits[j] = beta0[j]+X[i,] * beta[,j];
    }
    Y[i,] ~ dirichlet_multinomial(softmax(logits)*(1-psi)/psi);
  }
}

generated quantities {
  matrix[M, J] pipred;
  for (i in 1:M) {
    vector[J] logits;
    for (j in 1:J){
      logits[j] = beta0[j]+Xp[i,] * beta[,j];
    }
    pipred[i,] = transpose(softmax(logits));
  }
}"



MODEL=rstan::stan_model(model_code = STANcode)
StanDR=function(ldata, modelstan=MODEL){

#modelstan <-








  #extract y and X and W train for all treatment
  X=Y=W=Xpred=Wpred=trt=NULL
  init=1
  Ymat=matrix(0, nrow = sum(sapply(ldata$fitdata,  function(x) {length(x$Y)})),
              ncol = max(sapply(ldata$fitdata,  function(x) {max(x$Y+1)})))
  for(t in 1:length(ldata$fitdata) ){
    Y=c(Y,ldata$fitdata[[t]]$Y)
    nt=length(Y)
    X=rbind(X, ldata$fitdata[[t]]$XDiscPred)
    W=rbind(W, ldata$fitdata[[t]]$WProgn)

    trt=c(trt,rep(t-1,length(ldata$fitdata[[t]]$Y)))

    #transform y in matrix


  }

  Ymat[cbind(1:nt, Y+1)]=1
  init=init+nt
  #train
  df <- data.frame(X, W, trt = trt)

  # Rename columns
  w_terms <- paste0("w", seq_len(ncol(W)))
  d_terms <- paste0("d", seq_len(ncol(X)))

  colnames(df)[1:ncol(X)] <-d_terms
  colnames(df)[(ncol(X)+1):(ncol(X)+ncol(W))] <- w_terms



  Xpred=rbind(Xpred, ldata$postpredictivedata$XDiscPred)
  Wpred=rbind(Wpred, ldata$postpredictivedata$WProgn)
  npred=nrow(Xpred)
  #trtpred are fictitious, the observations are repeated, the 1st part of df will serves to predict treatment 1
  #2nd to predict treatment 2 and so on...
  trtpred=rep((1:length(ldata$fitdata))-1,each=npred)

  dfpred <- data.frame(Xpred, Wpred)
  dfpred=do.call(rbind, replicate(length(ldata$fitdata), dfpred, simplify = FALSE))
  dfpred=data.frame(dfpred, trt=factor(trtpred))

  colnames(dfpred)[1:ncol(X)] <- d_terms
  colnames(dfpred)[(ncol(X)+1):(ncol(X)+ncol(W))] <- w_terms



  df2=rbind(df,dfpred)
  df2[grep("^d\\d+$|^trt$", names(df2))] <- lapply(df2[grep("^d\\d+$|^trt$", names(df2))], as.factor)

  df=df2 [1:nrow(df),]
  dpred=df2[(nrow(df)+1):nrow(df2),]

  nameswtoremove=d_terms[unlist(lapply(df[grep("^d\\d+$", names(df))], function(x){length(levels(x))}))==1]

  d_termsr=setdiff(d_terms,nameswtoremove)

  # Build formula: all d* + trt + all w* + all w*:trt
  if(length(levels(df$trt))>1){
    fml <- as.formula(paste(
      "~ (", paste(d_termsr, collapse = " + "), ") * trt +",
     paste(w_terms, collapse = " + ")
  ))
    }else{
    fml <- as.formula(paste(
      "~ (", paste(d_termsr, collapse = " + "), ")   +",
      paste(w_terms, collapse = " + ")))
  }

  # Design matrix intercept
  Xmat = model.matrix(fml, df)[,-1]
dim(Xmat)





  #dfpred[grep("^d\\d+$|^trt$", names(dfpred))] <- lapply(dfpred[grep("^d\\d+$|^trt$", names(dfpred))], as.factor)
  Xpred = model.matrix(fml, dpred)[,-1]

  dim(Xpred)

  ##sampling
  Data = list(N = nrow(Xmat),  P = ncol(Xmat),
                 M = nrow(Xpred),
                 X = Xmat,
                 Xp = Xpred,
                 Y = Ymat, J = ncol(Ymat),
                 sd_prior = 1.0, psi = 0.01)

  time=Sys.time()
  WU=1000
  ITm=2000
  Stanfit=sampling(modelstan,data=Data, cores = 2, iter = ITm,
                   chains = 2, verbose = F, warmup = WU, seed = 121#,
        #           control = list(max_treedepth = 50, adapt_delta = 0.9999)
                )
  time=Sys.time()-time
#    fit <- rstan::stan(file = "R/dmhs-scripts/model.stan", data = ss_data, cores = 1, iter = 1000,
 #                   chains = 1, verbose = T, warmup = 200, seed = 121,
  #                   control = list(max_treedepth = 15, adapt_delta = 0.995))#995))


  # Sample ypred

params=rstan::extract(Stanfit)
  Ypred=t(apply(params$pipred, c(1,2), SampleMult1))-1

#return(Stanfit)
  return(list(Fit=Stanfit,ypred=Ypred,Time=time))
  }

