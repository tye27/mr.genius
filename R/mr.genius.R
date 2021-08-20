#' Simulate dataset
#'
#' @param m Number of SNPs
#' @param n Sample size
#' @param beta0 True causal effect
#' @param gamma A parameter that controls the magnitude of heteroscedasticity, i.e., the identification strength
#' @param case Simulation scenarios used in Table 1 in Ye et al., (2021)
#'
#' @return A list
#' \describe{
#' \item{z}{A n*m matrix of SNPs.}
#' \item{a}{A n-dimentional vector for the exposure.}
#' \item{y}{A n-dimentional vector for the outcome.}
#' }
#' @references Ting Ye, Zhonghua Liu, Baoluo Sun, and Eric Tchetgen Tchetgen (2021). GENIUS-MAWII: For Robust Mendelian Randomization with Many Weak Invalid Instruments.\url{https://arxiv.org/abs/2107.06238}.
#'
#' @import stats
#' @export
#' @examples
#' df<-data_gen(m=100,n=1e5,beta=0.4,gamma=0.1)
data_gen<-function(m,n,beta0,gamma,case="case2"){
  z_coef<-data.frame(z_x=rep(gamma,m))
    z<-matrix(nrow=n,ncol=m)
    for(j in 1:m){
      tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
      z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
    }
  if(case=="case2"){
    u<-rnorm(n,0,1)
    a<-z%*% rep(1,m) +u+rnorm(n,0,z %*% z_coef$z_x) # *a* must be heteroscedestic
    y<-beta0*a+ z%*% rep(1,m)+2*u+rnorm(n,0,2)
  }else if (case=="case1"){
    u<-rnorm(n,0,1)
    a<-z%*% rep(1,m) +u+rnorm(n,0,sqrt(z %*% z_coef$z_x)) # *a* must be heteroscedestic
    y<-beta0*a+ z%*% rep(1,m)+2*u+rnorm(n,0,2)
  }
  return(list(z=z, a=a, y=y))
}

#' Simulate a dataset with covariates
#'
#' @param m Number of SNPs
#' @param n Sample size
#' @param beta0 True causal effect
#' @param gamma A parameter that controls the magnitude of heteroscedasticity, i.e., the identification strength
#'
#' @return A list
#' \describe{
#' \item{z}{A n*m matrix of SNPs.}
#' \item{a}{A n-dimentional vector for the exposure.}
#' \item{y}{A n-dimentional vector for the outcome.}
#' \item{x}{A n-dimentional vector for the covariate.}
#' }
#' @references Ting Ye, Zhonghua Liu, Baoluo Sun, and Eric Tchetgen Tchetgen (2021). GENIUS-MAWII: For Robust Mendelian Randomization with Many Weak Invalid Instruments.\url{https://arxiv.org/abs/2107.06238}.
#'
#' @import stats
#' @export
#' @examples
#' df<-data_gen_x(m=100,n=1e5,beta=0.4,gamma=0.1)
data_gen_x<-function(m,n,beta0,gamma){
  z_coef<-data.frame(z_x=rep(gamma,m))
  z<-matrix(nrow=n,ncol=m)
  for(j in 1:m){
    tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
    z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
  }
  u<-rnorm(n,1,1)
  x<-rnorm(n,0,1)
  a<-z%*% rep(1,m) +u*x+rnorm(n,0,z %*% z_coef$z_x) # *a* must be heteroscedestic
  y<-beta0*a+ z%*% rep(1,m)*2+2*u*x+rnorm(n,0,2)
  return(list(z=z, a=a, y=y, x=x))
}

data_gen_mis<-function(m,n,beta0,gamma, xi_a,xi_y,prop){
  z_coef<-data.frame(z_x=rep(gamma,m))
  z<-matrix(nrow=n,ncol=m)
  for(j in 1:m){
    tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
    z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
  }
  xi_a_vec<-c(rep(xi_a,round(m*prop)),rep(0,m-round(m*prop)))
  xi_y_vec<-c(rep(xi_y,round(m*prop)),rep(0,m-round(m*prop)))
  #Case 2
    u<-rnorm(n,0,1)
    a<-z%*% rep(1,m) + u+u*z%*% xi_a_vec +rnorm(n,0,z %*% z_coef$z_x) # *a* must be heteroscedestic
    y<-beta0*a+ z%*% rep(1,m)+2*u+u*z%*% xi_y_vec+rnorm(n,0,2)
  return(list(z=z, a=a, y=y))
}


gf_yc_ac<-function(beta, df){ #  cont. outcome, cont. exposure
  a.pred<-lm(a~z,data=df)$fitted.values
  z.demean<-apply(df$z,2,function(x) x-mean(x)) # n*m matrix
  y.pred<-lm(y~z,data=df)$fitted.values
  delta<-(df$a-a.pred)*(df$y-y.pred)-beta*(df$a-a.pred)^2
  g<-as.vector(delta-mean(delta))*z.demean# n*m matrix # -mean(delta) is important
  return(g)
}

gf_yc_ac_x<-function(beta, df,omega.hat,theta.hat){ #  cont. outcome, cont. exposure, with X
  a.pred<-lm(a~z+x,data=df)$fitted.values
  z.demean<-apply(df$z,2,function(z) z-lm(z~df$x)$fitted.values) # n*m matrix
  y.pred<-lm(y~z+x,data=df)$fitted.values
  delta<-(df$a-a.pred)*(df$y-y.pred)-beta*(df$a-a.pred)^2
  g<-as.vector(delta-omega.hat+beta*theta.hat)*z.demean# n*m matrix # -mean(delta) is important
  return(g)
}

cue_obj_gf<-function(param, df, x=c(FALSE,TRUE),omega.hat=NULL, theta.hat=NULL){ # CUE objective function
  n<-length(df$a)
  if(x==FALSE){
    g<-gf_yc_ac(param, df)
  }else if(x==TRUE){
    g<-gf_yc_ac_x(param, df,omega.hat, theta.hat)
  }
  g_hat<-matrix(nrow=dim(g)[2], ncol=1)
  g_hat[,1]<-apply(g,2,mean)
  Omega_hat<-t(g)%*% g /n  # m*m matrix
  obj<-t(g_hat)%*% solve(Omega_hat) %*% g_hat/2
  return(obj)
}

deriv_gf<-function(param, df,x=c(FALSE,TRUE),omega.hat=NULL, theta.hat=NULL){ # calculates the derivativ of CUE objective function, used for rootsolve
  n<-length(df$a)
  if(x==FALSE){
    g<-gf_yc_ac(param,df)
    g_hat<-matrix(nrow=dim(g)[2], ncol=1)
    g_hat[,1]<-apply(g,2,mean)
    Omega_hat<-t(g)%*%g/n
    z.demean<-apply(df$z,2,function(x) x-mean(x)) # n*m matrix
    a.pred<-lm(a~z,data=df)$fitted.values
    G<-as.vector(-(df$a-a.pred)^2+mean((df$a-a.pred)^2 ))*z.demean  # n*m matrix
  }else if(x==TRUE){
    g<-gf_yc_ac_x(param,df)
    g_hat<-matrix(nrow=dim(g)[2], ncol=1)
    g_hat[,1]<-apply(g,2,mean)
    Omega_hat<-t(g)%*%g/n
    z.demean<-apply(df$z,2,function(z) z-lm(z~df$x)$fitted.values) # n*m matrix
    a.pred<-lm(a~z+x,data=df)$fitted.values
    G<-as.vector(-(df$a-a.pred)^2+theta.hat)*z.demean  # n*m matrix
  }
  G_hat<-matrix(nrow=dim(G)[2], ncol=1)
  G_hat[,1]<-apply(G,2,mean)
  Omega_hat_1deriv<-(t(G)%*% g+t(g)%*% G)/n
  res<-(t(G_hat)- t(g_hat) %*% solve(Omega_hat) %*% Omega_hat_1deriv/2) %*% solve(Omega_hat) %*% g_hat
  return(res)
}

cue_var_gf<-function(param,df,x=c(FALSE,TRUE),omega.hat=NULL, theta.hat=NULL){ # variance estimator and strength measures
  n<-length(df$a)
  if(x==FALSE){
    g<-gf_yc_ac(param,df)
    g_hat<-matrix(nrow=dim(g)[2], ncol=1)
    g_hat[,1]<-apply(g,2,mean)
    Omega_hat<-t(g)%*%g/n
    z.demean<-apply(df$z,2,function(x) x-mean(x)) # n*m matrix
    a.pred<-lm(a~z,data=df)$fitted.values
    G<-as.vector(-(df$a-a.pred)^2+mean((df$a-a.pred)^2 ))* z.demean  # n*m matrix
  }else if(x==TRUE){
    g<-gf_yc_ac_x(param,df,omega.hat, theta.hat)
    g_hat<-matrix(nrow=dim(g)[2], ncol=1)
    g_hat[,1]<-apply(g,2,mean)
    Omega_hat<-t(g)%*%g/n
    z.demean<-apply(df$z,2,function(z) z-lm(z~df$x)$fitted.values) # n*m matrix
    a.pred<-lm(a~z+x,data=df)$fitted.values
    G<-as.vector(-(df$a-a.pred)^2+theta.hat)*z.demean  # n*m matrix
  }
  G_hat<-matrix(nrow=dim(G)[2], ncol=1)
  G_hat[,1]<-apply(G,2,mean)
  Omega_hat_1deriv<-(t(G)%*% g+t(g)%*% G)/n
  Omega_hat_2deriv<-2*(t(G) %*% G)/n
  H_hat<-t(G_hat) %*% solve(Omega_hat) %*% G_hat -
    2*t(G_hat) %*% solve(Omega_hat) %*% Omega_hat_1deriv %*% solve(Omega_hat) %*% g_hat+
    t(g_hat) %*% solve(Omega_hat) %*% Omega_hat_1deriv %*% solve(Omega_hat) %*% Omega_hat_1deriv %*% solve(Omega_hat) %*% g_hat-
    1/2* t(g_hat) %*% solve(Omega_hat) %*% Omega_hat_2deriv %*% solve(Omega_hat) %*% g_hat
  D_hat<-G_hat-(t(G)%*% g) %*% solve(Omega_hat) %*% g_hat/n # m*1
  V_hat<-solve(H_hat) %*% t(D_hat) %*% solve(Omega_hat) %*% D_hat %*% solve(H_hat)
  var_est<-V_hat[1,1]/n # variance estimator for CUE
  strength_full<-as.numeric(min(eigen(n*H_hat)$values))
  strength_simple1<-as.numeric(min(eigen(n* t(G_hat) %*% solve(Omega_hat) %*% G_hat)$values)) # strength measure for standard GMM
  strength_simple2<-as.numeric(min(eigen(n* t(G_hat) %*% G_hat)$values))
  return(list(var_est=var_est,strength_full=strength_full, strength_simple1=strength_simple1, strength_simple2=strength_simple2))
}

#' Main function for GENIUS-MAWII
#'
#' @param z A n*m matrix of SNPs, where n is the sample size, m is the number of SNPs
#' @param a A n-dimentional vector for the exposure
#' @param y A n-dimentional vector for the outcome
#' @param x A n-dimentional vector for the covariate. Default is NULL, when there is no covariates
#' @param alpha Confidence interval has level 1-alpha. Default is 0.05
#' @param diagnostics Should the function returns the residual plot for assumption diagnosis. Default is FALSE
#'
#' @details Computation is fast in the case of no \code{x}. When there are observed covariates \code{x}, there are functions to be estimated by nonparametric kernel, which makes computation slow.
#' @return A list
#' \describe{
#' \item{beta.hat}{Estimated causal effect}
#' \item{beta.se}{Standard error of \code{beta.hat}}
#' \item{ci}{A 1-alpha confidence interval}
#' \item{condition}{A measure that needs to be large for reliable asymptotic approximation based on the GENIUS-MAWII estimator. It is recommended to be greater than 50}
#' }
#'
#' @references Ting Ye, Zhonghua Liu, Baoluo Sun, and Eric Tchetgen Tchetgen (2021). GENIUS-MAWII: For Robust Mendelian Randomization with Many Weak Invalid Instruments.\url{https://arxiv.org/abs/2107.06238}.
#'
#' @import stats ggplot2 np rootSolve
#' @export
#'
#' @examples
#' df<-data_gen(m=20,n=1e5,beta=0.4,gamma=0.1)
#' mr.genius(df$z,df$a,df$y,diagnostics=TRUE)
mr.genius<-function(z,a,y,x=NULL,alpha=0.05, diagnostics=FALSE){ # this function integrates previous functions
  if(is.null(x)){
    df<-list(z=z,a=a,y=y)
    use.x<-FALSE
    omega.hat<-NULL
    theta.hat<-NULL
  }else{
    df<-list(z=z,a=a,y=y,x=x)
    use.x<-TRUE
    a.pred<-lm(a~z+x,data=df)$fitted.values
    z.demean<-apply(df$z,2,function(z) z-lm(z~df$x)$fitted.values) # n*m matrix
    y.pred<-lm(y~z+x,data=df)$fitted.values
    bw.omega <- npregbw((df$a-a.pred)*(df$y-y.pred)~df$x)
    bw.theta <- npregbw((df$a-a.pred)^2~df$x)
    omega.hat<-predict(npreg(bw.omega))
    theta.hat<-predict(npreg(bw.theta))
  }
  interval<-c(-10,10)
  opt.res<-optimize(cue_obj_gf, df=df,x=use.x,omega.hat=omega.hat, theta.hat=theta.hat,interval=interval)$minimum
  if(min(abs(opt.res-interval))<0.01){
    opt.res<-uniroot(deriv_gf,interval=c(-0.5,0.5),extendInt="yes",df=df,x=use.x,omega.hat=omega.hat, theta.hat=theta.hat)$root # need to specify a small interval to make the end points of opposite signs.
  }
  cue.obj.value<-cue_obj_gf(opt.res,df,x=use.x,omega.hat=omega.hat, theta.hat=theta.hat)[1,1]
  tmp<-cue_var_gf(opt.res,df,x=use.x,omega.hat=omega.hat, theta.hat=theta.hat)
  se<-as.numeric(sqrt(tmp$var_est))
  z.alpha<-qnorm(alpha/2,lower.tail = FALSE)
  ci<-c(opt.res[1]-z.alpha*se, opt.res[1]+z.alpha*se)
  #CP<-1*(beta0> ci[1] & beta0<ci[2])
  out<-list(beta.hat=opt.res[1], beta.se=se, ci=ci, #cue.obj.value=cue.obj.value, #CP=CP, #f.statistic=f.statistic,
            condition=tmp$strength_full)
  if(diagnostics){
    if(is.null(x)){
      a.pred<-lm(a~z,data=df)$fitted.values
      t<-(df$a-a.pred)*(df$y-opt.res[1]*df$a)
      t<-t-mean(t)
      p_df<-data.frame(a.pred,t)
      out$diagnostics<-ggplot(p_df)+aes(x=a.pred^2,y=t)+geom_point()+
        geom_smooth(method = "lm",formula = y ~ splines::bs(x, 3))+xlab("A.predicted.squared")+ylab("residual")
    }else{
      tmp<-omega.hat-opt.res[1]*theta.hat
      t<-(df$a-a.pred)*(df$y-opt.res[1]*df$a)
      t<-t-tmp
      p_df<-data.frame(a.pred,t)
      out$diagnostics<-ggplot(p_df)+aes(x=a.pred^2,y=t)+geom_point()+
        geom_smooth(method = "lm",formula = y ~ splines::bs(x, 3))+xlab("A.predicted.squared")+ylab("residual")
    }
  }
  return(out)
}

individ_to_summary<-function(df,out.method=c("lm","logistic")){
  full_df<-data.frame(beta.exposure=numeric(dim(df$z)[2]))
  full_df$beta.exposure<-sapply(1:dim(df$z)[2], function(i) {
    coef <- lm(df$a ~ df$z[,i])$coef[-1]
    return(coef)
  })
  full_df$beta.outcome<-sapply(1:dim(df$z)[2], function(i) {
    coef <- lm(df$y ~ df$z[,i])$coef[-1]
    return(coef)
  })
  full_df$se.exposure<-sapply(1:dim(df$z)[2], function(i) {
    se <- summary( lm(df$a ~ df$z[,i]))$coefficients[2,2]
    return(se)
  })
  full_df$se.outcome<-sapply(1:dim(df$z)[2], function(i) {
    se <- summary( lm(df$y ~ df$z[,i]))$coefficients[2,2]
    return(se)
  })
  full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
  full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  return(full_df)
}

