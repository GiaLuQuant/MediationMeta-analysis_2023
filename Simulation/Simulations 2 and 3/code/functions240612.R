### code for vcov.osmasem provided by Dr. M. Cheung
vcov.osmasem <- function(object, select=c("fixed", "all", "random"),
                         robust=FALSE, ...) {
  if (!is.element("osmasem", class(object)))
    stop("\"object\" must be an object of class \"osmasem\".")
  # labels of the parameters
  ## my.name <- summary(object$mx.fit)$parameters$name
  my.name <- names(omxGetParameters(object$mx.fit))
  my.name <- my.name[!is.na(my.name)]
  
  ## index for untransformed random effects (not the correct ones!)
  index.random <- grep("Tau1_|Cor_", my.name)
  
  ## Fix a bug reported by Zhiming Lu with index.random=integer(0)
  #when there is no random effects with RE.type="Zero"
  select <- match.arg(select)
  
  if (length(index.random) != 0) {
    switch( select,
            ## all = my.name <- my.name,
            fixed =  my.name <- my.name[-index.random],
            random = my.name <- my.name[index.random]
    )
  }
  if (robust) {
    out <- suppressMessages(imxRobustSE(object$mx.fit,
                                        details=TRUE))$cov
  } else {
    out <- vcov(object$mx.fit)
  }
  if (select!="fixed") warning("\"Tau1_xx\" is not the variance
component of the random effects.")
  out[my.name, my.name]
}

# combined variance
c.var = function(var1, var2, n1, n2, m1, m2){
  m = ((n1-1)*m1+(n2-1)*m2)/(n1+n2-2)
  d1 = m1-m; d2=m2-m
  return(((n1-1)*(var1+d1^2)+(n2-1)*(var2+d2^2))/(n1+n2-2))
}

# point biserial correlation 
r_pb = function(MD, sd.cb, n1, n2){
  return(MD*sqrt(((n1-1)*(n2-1))/(n1+n2-2)^2)/sd.cb)
}

# generating a new matrix
Gen.matrix <- function(){
  cm.ps=matrix(NA, 3,3)
  colnames(cm.ps) = rownames(cm.ps) = c('X', 'M', 'Y')
  diag(cm.ps) = 1
  return(cm.ps)
}

# generating sample sizes of individual studies
genN <- function(Nbar, Nsd, K){
  mu = 0.6*Nbar
  size = mu^2/(Nsd^2-mu)
  N = rnbinom(K,size=size, mu=0.6*Nbar) + 0.4*Nbar
  return(N)
}

data_gen <- function(a.s, cp.s, b.s, c.s,Rho,X,n.exp,n.ctl,postVar,p=.5){
  # variances for variables
  sigma.X2 = p*(1-p); sigma.1T2 = sigma.1C2 = sigma.2C2 = 1; sigma.2T2 = postVar
  sigma.Mcs.T2 = sigma.Ycs.T2 = sigma.1T2 - 2*Rho*sqrt(sigma.1T2*sigma.2T2) + sigma.2T2
  sigma.Mcs.C2 = sigma.Ycs.C2 = sigma.1C2 - 2*Rho*sqrt(sigma.1C2*sigma.2C2) + sigma.2C2
  
  # calculating group means 
  MD.M.k = (2*a.s/sqrt(1-a.s^2))*sqrt((sigma.Mcs.T2+sigma.Mcs.C2)/2)
  MD.Y.k = (2*c.s/sqrt(1-c.s^2))*sqrt((sigma.Ycs.T2+sigma.Ycs.C2)/2)
  
  # variances of change scores of M and Y
  sigma.Mcs2 = c.var(var1 = sigma.Mcs.C2, var2 = sigma.Mcs.T2, n1 = n.ctl,n2 = n.exp, m1 = 0, m2 = MD.M.k)
  sigma.Ycs2 = c.var(var1 = sigma.Ycs.C2, var2 = sigma.Ycs.T2, n1 = n.ctl,n2 = n.exp, m1 = 0, m2 = MD.Y.k)
  
  # unstandardized coefficients
  a.u = a.s* sqrt(sigma.Mcs2/sigma.X2)
  b.u = b.s* sqrt(sigma.Ycs2/sigma.Mcs2)
  cp.u = cp.s* sqrt(sigma.Ycs2/sigma.X2)
  c.u = a.u*b.u + cp.u
  
  # generating M
  Rho.pp.T = matrix(c(sigma.1T2,Rho*sqrt(sigma.1T2*sigma.2T2),
                      Rho*sqrt(sigma.1T2*sigma.2T2),sigma.2T2),2,2)
  Rho.pp.C = matrix(c(sigma.1C2,Rho*sqrt(sigma.1C2*sigma.2C2),
                      Rho*sqrt(sigma.1C2*sigma.2C2),sigma.2C2),2,2)
  
  mu.M.exp <- c(0, MD.M.k); mu.M.ctl <- c(0,0)
  
  M.exp <- mvrnorm(n.exp, mu.M.exp, Rho.pp.T)
  M.ctl <- mvrnorm(n.ctl, mu.M.ctl, Rho.pp.C)
  M.T.cs <- M.exp[,2] - M.exp[,1]
  M.C.cs <- M.ctl[,2] - M.ctl[,1]
  M.cs <- c(M.C.cs, M.T.cs)
  M.2 <- c(M.ctl[,2], M.exp[,2])
  M.1 <- c(M.ctl[,1], M.exp[,1])
  
  # generating Y
  sigma.eYcs2C = sigma.Ycs.C2 - b.u^2*sigma.Mcs.C2
  sigma.eYcs2T = sigma.Ycs.T2 - b.u^2*sigma.Mcs.T2
  
  e.Ycs.T = rnorm(n.exp, 0, sqrt(sigma.eYcs2T))
  e.Ycs.C = rnorm(n.ctl, 0, sqrt(sigma.eYcs2C))
  
  Y.T.cs = MD.Y.k - b.u*MD.M.k + b.u*M.T.cs + e.Ycs.T
  Y.C.cs = b.u*M.C.cs + e.Ycs.C
  Y.cs = c(Y.C.cs, Y.T.cs)
  
  beta.ycs1.T = (Rho*sqrt(sigma.1T2*sigma.2T2)-sigma.1T2)/sigma.Ycs.T2
  beta.ycs1.C = (Rho*sqrt(sigma.1C2*sigma.2C2)-sigma.1C2)/sigma.Ycs.C2
  
  i.ycs1.T = -beta.ycs1.T*MD.Y.k
  i.ycs1.C = 0
  
  e.Y1.T = rnorm(n.exp,0, sqrt(sigma.1T2 - beta.ycs1.T^2*sigma.Ycs.T2))
  e.Y1.C = rnorm(n.ctl,0, sqrt(sigma.1C2 - beta.ycs1.C^2*sigma.Ycs.C2))
  
  Y.T.1 = i.ycs1.T + beta.ycs1.T*Y.T.cs + e.Y1.T
  Y.C.1 = i.ycs1.C + beta.ycs1.C*Y.C.cs + e.Y1.C
  
  Y.1 = c(Y.C.1, Y.T.1)
  Y.2 = Y.1+Y.cs 
  
  dat <- list(X=X,Mcs = M.cs, Ycs = Y.cs, M.2 = M.2, Y.2 = Y.2,
              M.T.1 = M.exp[,1], M.T.2 = M.exp[,2], 
              M.C.1 = M.ctl[,1], M.C.2 = M.ctl[,2], 
              Y.T.1 = Y.T.1, Y.T.2 = Y.2[X==1], 
              Y.C.1 = Y.C.1, Y.C.2 = Y.2[X==0])
  return(dat)
}

gen_post_heter <- function(a.s, cp.s, b.s, Rho, beta1, N_l, K_l, postVar, hetInd, p=.5){
  
  a.s = rep(a.s, K_l); b.s = rep(b.s, K_l); cp.s = rep(cp.s, K_l)
  
  n.exp=n.ctl=N_l/2
  X = c(rep(0,n.ctl), rep(1,n.exp))
  
  Z = rnorm(K_l, 0, 1)
  
  if(hetInd){
    a.s = a.s+beta1*Z
    b.s = b.s+beta1*Z
  }else{
    cp.s = cp.s+beta1*Z
  }
  
  c.s = a.s*b.s+cp.s
  
  a.est = c(); cp.est = c(); b.est = c()
  
  for(i in 1:K_l){
    dat = try(data_gen(a.s[i], cp.s[i], b.s[i], c.s[i], Rho,
                       X, n.exp, n.ctl, postVar, p=.5))
    if(is.na(dat$Ycs[1])){
      a.est = c(a.est, NA); b.est = c(b.est, NA); cp.est = c(cp.est, NA)
    }else{
      eq1 = try(lm(dat$M.2~X))
      if(inherits(eq1, 'try-error')){
        a.est = c(a.est, NA)
      }else{
        a.est = c(a.est, coef(eq1)['X']*(sd(X)/sd(dat$M.2)))
      }
      eq2 = try(lm(dat$Y.2~dat$M.2+X))
      if(inherits(eq2, 'try-error')){
        b.est = c(b.est, NA); cp.est = c(cp.est, NA)
      }else{
        cp.est = c(cp.est, coef(eq2)['X']*(sd(X)/sd(dat$Y.2)))
        b.est = c(b.est, coef(eq2)['dat$M.2']*(sd(dat$M.2)/sd(dat$Y.2)))
      }
    }
  }
  if(hetInd){
    beta1.ps = coef(lm(a.est~Z))['Z']
    beta2.ps = coef(lm(b.est~Z))['Z']
  }else{
    beta1.ps = coef(lm(cp.est~Z))['Z']
    beta2.ps = NA}
  
  return(c(mean(a.est), mean(cp.est), mean(b.est),beta1.ps, beta2.ps))
}


dg.ps <- function(StCoef, K, Tau, Nbar, Rho, postVar, beta1, Missing, MissingR, hetInd, p=.5){
  # generate sample size for individual studies
  if(Nbar==80){
    N = genN(Nbar, 23, K)
  }else if(Nbar==120){N = genN(Nbar, 30, K)}
  
  # generating the continuous moderator
  Z = rnorm(K, 0, 1)
  
  # generating true values with heterogeneity for individual studies
  
  if(hetInd){
    cp.s = rep(StCoef['cp'], K)
    a.s = StCoef['a']+beta1*Z; b.s = StCoef['b']+beta1*Z
  }else{
    cp.s = StCoef['cp']+beta1*Z
    a.s = rep(StCoef['a'], K); b.s = rep(StCoef['b'],K)
    }
  c.s = a.s*b.s+cp.s
  r.XY = c.s; r.XM = a.s; r.MY = a.s*cp.s+b.s
  
  para.h = matrix(NA, nrow = K, ncol=3); colnames(para.h) = c('r.XM', 'r.XY', 'r.MY')
  
  for(i in 1:K){
    para.h[i,'r.XM'] = rnorm(1, r.XM[i], Tau)
    para.h[i,'r.XY'] = rnorm(1, r.XY[i], Tau)
    para.h[i,'r.MY'] = rnorm(1, r.MY[i], Tau)
  }
  
  data4MA <- list(CS=list(),PS=list(), N = c())
  data4MA[['N']] <- N
  
  for(ki in 1:K){
    a.s = para.h[ki,'r.XM']; c.s=para.h[ki,'r.XY']
    b.s = (para.h[ki,'r.MY']-a.s*c.s)/(1-a.s^2)
    cp.s = (c.s-a.s*b.s)
    # generating data
    X = sort(rbinom(N[ki], 1, p))
    n.exp=sum(X==1); n.ctl=sum(X==0)
    dat = data_gen(a.s=a.s, cp.s=cp.s, b.s=b.s, c.s=c.s, Rho=Rho, X=X,
                   n.exp=n.exp, n.ctl=n.ctl, postVar=postVar)
    
    for(i in 1:2){
      data4MA[[i]][[ki]] <- Gen.matrix()
    }
    
    ## calculating effect sizes using the 5 methods
    # CS method
    rXM.cs = r_pb(MD = mean(dat$M.T.2-dat$M.T.1)-mean(dat$M.C.2-dat$M.C.1), sd = sd(dat$Mcs), n1 = n.ctl, n2 = n.exp)
    rXY.cs = r_pb(MD = mean(dat$Y.T.2-dat$Y.T.1)-mean(dat$Y.C.2-dat$Y.C.1), sd = sd(dat$Ycs), n1 = n.ctl, n2 = n.exp)
    
    data4MA[['CS']][[ki]][1,2]=data4MA[['CS']][[ki]][2,1]=rXM.cs
    data4MA[['CS']][[ki]][1,3]=data4MA[['CS']][[ki]][3,1]=rXY.cs
    data4MA[['CS']][[ki]][2,3]=data4MA[['CS']][[ki]][3,2]=cor(dat$Mcs, dat$Ycs)
    
    # PS method
    rXM.ps = r_pb(MD = mean(dat$M.T.2)-mean(dat$M.C.2), sd = sd(dat$M.2), n1 = n.ctl, n2 = n.exp)
    rXY.ps = r_pb(MD = mean(dat$Y.T.2)-mean(dat$Y.C.2), sd = sd(dat$Y.2), n1 = n.ctl, n2 = n.exp)
    data4MA[['PS']][[ki]][1,2]=data4MA[['PS']][[ki]][2,1]=rXM.ps
    data4MA[['PS']][[ki]][1,3]=data4MA[['PS']][[ki]][3,1]=rXY.ps
    data4MA[['PS']][[ki]][2,3]=data4MA[['PS']][[ki]][3,2]=cor(dat$M.2, dat$Y.2)
  }
  data4MA$Mod = Z
  
  # Missing completely at random
  if(Missing){
    if(MissingR==0.4){
      MissingR.rXM = 0.4; MissingR.rMY = 0.8
      missing.rXM.index = sort(sample(K, MissingR.rXM*K))
      missing.rMY.index = sort(sample(K, MissingR.rMY*K))
    }else{
      MissingR.M = 0.9
      missing.rXM.index = missing.rMY.index = sort(sample(K, MissingR.M*K))
    }
    
    for(i in missing.rXM.index){
      data4MA[['CS']][[i]][1,2] = data4MA[['CS']][[i]][2,1] = NA
      data4MA[['PS']][[i]][1,2] = data4MA[['PS']][[i]][2,1] = NA
    }
    for(i in missing.rMY.index){
      data4MA[['CS']][[i]][2,3] = data4MA[['CS']][[i]][3,2] = NA
      data4MA[['PS']][[i]][2,3] = data4MA[['PS']][[i]][3,2] = NA
    }
  }
  return(data4MA)
}

extMA <- function(Fit,hetInd){
  s = summary(Fit)
  
  if(hetInd){
    Est = s$parameters[c(1:5),'Estimate']  # a cp b beta1
    
    c.est = Est[1]*Est[3]+Est[2]
    ab.est = Est[1]*Est[3]
    
    se = s$parameters[1:5,'Std.Error']
    pValue = s$parameters[1:5,'Pr(>|z|)']
    
    cov.ab <- vcov.osmasem(Fit)['b','a']
    cov.a <- vcov.osmasem(Fit)['a','a']
    cov.b <- vcov.osmasem(Fit)['b','b']
    se.apb.delta <- sqrt(4*(Est[3]^2*cov.a + Est[1]^2*cov.b + 2*Est[1]*Est[3]*cov.ab))
    
  }else{
    Est = s$parameters[c(1:4, NA),'Estimate']  # a cp b beta1
    
    c.est = Est[1]*Est[3]+Est[2]
    ab.est = Est[1]*Est[3]
    
    se = c(s$parameters[1:4,'Std.Error'], NA)
    pValue = c(s$parameters[1:4,'Pr(>|z|)'], NA)
    
    cov.ab <- vcov.osmasem(Fit)['b','a']
    cov.a <- vcov.osmasem(Fit)['a','a']
    cov.b <- vcov.osmasem(Fit)['b','b']
    se.apb.delta <- sqrt(4*(Est[3]^2*cov.a + Est[1]^2*cov.b + 2*Est[1]*Est[3]*cov.ab))
  }

  infoDef = 1
  
  res = c(Est, c.est, ab.est, se, se.apb.delta, pValue, infoDef)
  return(res)
}

# study 1

ma <- function(StCoef, K, Tau, Nbar, Rho, nrep, postVar, beta1, hetInd, Missing = FALSE, MissingR = 0){
  RES = list()
  res.i = matrix(NA, nrep, 19)
  colnames(res.i) <- c('a.est', 'cp.est', 'b.est', 'beta1.est','beta2.est', 'c.est','ab.est',
                       'a.se','cp.se','b.se','beta1.se','beta2.se','ab.se.delta',
                       'a.p','cp.p','b.p','beta1.p','beta2.p','infoDef')
  for(i in 1:2){ 
    RES[[i]] = res.i
  };names(RES) <- c('CS','PS')
  
  # model
  MedM <- 'M ~ a*X
           Y ~ cp*X + b*M'
  
  RAM1 <- lavaan2RAM(MedM, obs.variables=c("X", "M", "Y"))
  RAM1$S[1,1] <- 1
  
  if(hetInd){
    Ax = matrix(c(0,0,0,
                  '0*data.mod',0,0,
                  0,'0*data.mod',0),3,3, byrow = T)
  }else{
    Ax = matrix(c(0,0,0,
                  0,0,0,
                  '0*data.mod',0,0),3,3, byrow = T)
  }

  for (i in 1:nrep){
    dat = try(dg.ps(StCoef, K, Tau, Nbar, Rho, postVar, beta1, Missing, MissingR, hetInd, p=0.5))
    
    if(inherits(dat, 'try-error')){
      RES[[1]][i,]=RES[[2]][i,]=RES[[3]][i,]=rep(NA, ncol(res.i))
    }else{
      for(j in 1:2){
        df <- Cor2DataFrame(dat[[j]], dat$N)
        df$data = data.frame(df$data, mod=dat$Mod,check.names=FALSE)
        fit0 = try(osmasem(model.name="mediationMA", RAM=RAM1, Ax = Ax, data=df))
        
        if(inherits(fit0, 'try-error')){
          res <- c(rep(NA, ncol(res.i)))
        }else{
          if(summary(fit0)$infoDefinite==FALSE){
            res <- extMA(fit0, hetInd) 
            res[ncol(res.i)] <- -1
          }else{res<- extMA(fit0, hetInd)}
        }
        RES[[j]][i,]=res
      }
    }
  }
  return(RES)
}


# function for deleting negative definite results
select_res <- function(data){
  nrow = nrow(data); ncol = ncol(data)
  for(i in 1:nrow){
    if (data[i,'infoDef']==-1){
      data[i,]=rep(NA,ncol)
    }
  }
  data = remove_empty(data, which='rows')
  return(data)
}

## performance measures
cal.eBias <- function(theta.hat, theta){
  if(theta==0){
    Ebias <- mean(theta.hat, na.rm=T)-theta
  }else{Ebias <- (mean(theta.hat, na.rm=T)-theta)/theta}
  return(Ebias)
}

# coverage rate
cal.CR <- function(SE, theta, theta.hat){
  cv.error = qnorm(1-0.05/2) * SE
  l = theta.hat-cv.error  # lower bound
  u = theta.hat+cv.error  # upper bound
  lcriterion = (theta >= l)
  ucriterion = (theta <= u)
  CR = mean(lcriterion*ucriterion, na.rm=T)
  return(CR)
}

sep.RR <- function(RR, para){
  if(para == 0){
    power = NA
    typeIerror = RR
  }else{
    power = RR
    typeIerror = NA
  }
  return(c(power,typeIerror))
}

# simulation results 
cal_res1 <- function(data,a.s,cp.s,b.s,beta1){
  p=.5
  sigma.X = sqrt(p*(1-p))
  # deleting negative-definite repetitions
  data.s <- select_res(data)
  
  bias = c(mean(data.s[,'a.est']-a.s, na.rm=T)/sigma.X, 
           mean(data.s[,'cp.est']-cp.s, na.rm=T)/sigma.X,
           mean(data.s[,'b.est']-b.s, na.rm=T),
           mean(data.s[,'ab.est']-a.s*b.s, na.rm=T)/sigma.X,
           mean(data.s[,'beta1.est']-beta1, na.rm=T))
  # calculating EBIAS
  Ebias = c(cal.eBias(data.s[,'a.est']/sigma.X, a.s/sigma.X), 
            cal.eBias(data.s[,'cp.est']/sigma.X, cp.s/sigma.X), 
            cal.eBias(data.s[,'b.est'],b.s),
            cal.eBias(data.s[,'ab.est']/sigma.X, a.s*b.s/sigma.X),
            cal.eBias(data.s[,'beta1.est'], beta1))
  
  # calculating empirical SEs of a, b, cp, and beta1
  empSE <- apply(data.s[,1:4],2,sd)
  seM <- apply(data.s[,8:11],2,mean, na.rm=T)
  seRbias <- c(cal.eBias(data.s[,'a.se'],empSE[1]), 
               cal.eBias(data.s[,'cp.se'],empSE[2]),
               cal.eBias(data.s[,'b.se'],empSE[3]),
               cal.eBias(data.s[,'beta1.se'],empSE[4]))
  
  # rejection rate (for a only, b only, and a*b)
  a.RR <- sep.RR(mean(data.s[,'a.p']<=0.05, na.rm=T), a.s)
  cp.RR <- sep.RR(mean(data.s[,'cp.p']<=0.05, na.rm=T), cp.s)
  b.RR <- sep.RR(mean(data.s[,'b.p']<=0.05, na.rm=T), b.s)
  beta1.RR <- sep.RR(mean(data.s[,'beta1.p']<=0.05, na.rm=T), beta1)
  
  Z.delta <- (data.s[,'ab.est']*2) / data.s[,'ab.se.delta']
  ab.p.delta <- 1-pnorm(abs(Z.delta))
  ab.RR.delta <- sep.RR(mean(ab.p.delta<=0.025, na.rm=T),a.s*b.s)
  
  # coverage rate
  a.CR = cal.CR(SE = data.s[,'a.se']*2, theta = a.s*2, theta.hat = data.s[,'a.est']*2)
  cp.CR = cal.CR(SE = data.s[,'cp.se']*2, theta = cp.s*2, theta.hat = data.s[,'cp.est']*2)
  b.CR = cal.CR(SE = data.s[,'b.se'], theta = b.s, theta.hat = data.s[,'b.est'])
  beta1.CR = cal.CR(SE = data.s[,'beta1.se'], theta = (beta1), theta.hat = data.s[,'beta1.est'])
  ab.CR.delta = cal.CR(SE = data.s[,'ab.se.delta']*2, theta = (a.s*b.s*2), theta.hat = data.s[,'ab.est']*2)
  
  # convergence rate & postive definite rate
  ConvRate <- 1 - mean(is.na(data[,'a.est']), na.rm=T)
  posDefRate <- mean(data[,'infoDef']==1, na.rm=T)
  
  # saving results
  res<- matrix(NA, 1, 43)
  res[1,] <- c(a.s,cp.s,b.s,beta1,
               bias, Ebias, empSE, seM, seRbias,
               a.CR, cp.CR, b.CR, beta1.CR, ab.CR.delta,
               a.RR, cp.RR, b.RR, beta1.RR, ab.RR.delta, ConvRate, posDefRate)
  colnames(res) <- c('a.para', 'cp.para', 'b.para', 'beta1.para',
                     'a.bias','cp.bias', 'b.bias', 'ab.bias','beta1.bias',
                     'a.Ebias','cp.Ebias', 'b.Ebias', 'ab.Ebias','beta1.Ebias',
                     'a.empSE', 'cp.empSE','b.empSE', 'beta1.empSE',
                     'a.meanSE', 'cp.meanSE','b.meanSE','beta1.meanSE',
                     'a.seRbias', 'cp.seRbias', 'b.seRbias','beta1.seRbias',
                     'a.CoverageRate','cp.CoverageRate', 'b.CoverageRate','beta1.CoverageRate','ab.CoverageRate', 
                     'a.power', 'a.typeIerror','cp.power', 'cp.typeIerror', 'b.power', 'b.typeIerror', 
                     'beta1.power', 'beta1.typeIerror', 'ab.power', 'ab.typeIerror', 
                     'ConvergenceRate','posDefRate')
  return(res)
}

cal_res2 <- function(data,a.s,cp.s,b.s,beta1,beta2){
  p=.5
  sigma.X = sqrt(p*(1-p))
  # deleting negative-definite repetitions
  data.s <- select_res(data)
  
  bias = c(mean(data.s[,'a.est']-a.s, na.rm=T)/sigma.X, 
           mean(data.s[,'cp.est']-cp.s, na.rm=T)/sigma.X,
           mean(data.s[,'b.est']-b.s, na.rm=T),
           mean(data.s[,'ab.est']-a.s*b.s, na.rm=T)/sigma.X,
           mean(data.s[,'beta1.est']-beta1, na.rm=T),
           mean(data.s[,'beta2.est']-beta2, na.rm=T))
  # calculating EBIAS
  Ebias = c(cal.eBias(data.s[,'a.est']/sigma.X, a.s/sigma.X), 
            cal.eBias(data.s[,'cp.est']/sigma.X, cp.s/sigma.X), 
            cal.eBias(data.s[,'b.est'],b.s),
            cal.eBias(data.s[,'ab.est']/sigma.X, a.s*b.s/sigma.X),
            cal.eBias(data.s[,'beta1.est'], beta1),
            cal.eBias(data.s[,'beta2.est'], beta2))
  
  # calculating empirical SEs of a, b, cp, and beta1
  empSE <- apply(data.s[,1:5],2,sd)
  seM <- apply(data.s[,8:12],2,mean, na.rm=T)
  seRbias <- c(cal.eBias(data.s[,'a.se'],empSE[1]), 
               cal.eBias(data.s[,'cp.se'],empSE[2]),
               cal.eBias(data.s[,'b.se'],empSE[3]),
               cal.eBias(data.s[,'beta1.se'],empSE[4]),
               cal.eBias(data.s[,'beta2.se'],empSE[5]))
  
  # rejection rate (for a only, b only, and a*b)
  a.RR <- sep.RR(mean(data.s[,'a.p']<=0.05, na.rm=T), a.s)
  cp.RR <- sep.RR(mean(data.s[,'cp.p']<=0.05, na.rm=T), cp.s)
  b.RR <- sep.RR(mean(data.s[,'b.p']<=0.05, na.rm=T), b.s)
  beta1.RR <- sep.RR(mean(data.s[,'beta1.p']<=0.05, na.rm=T), beta1)
  beta2.RR <- sep.RR(mean(data.s[,'beta2.p']<=0.05, na.rm=T), beta2)
  
  Z.delta <- (data.s[,'ab.est']*2) / data.s[,'ab.se.delta']
  ab.p.delta <- 1-pnorm(abs(Z.delta))
  ab.RR.delta <- sep.RR(mean(ab.p.delta<=0.025, na.rm=T),a.s*b.s)
  
  # coverage rate
  a.CR = cal.CR(SE = data.s[,'a.se']*2, theta = a.s*2, theta.hat = data.s[,'a.est']*2)
  cp.CR = cal.CR(SE = data.s[,'cp.se']*2, theta = cp.s*2, theta.hat = data.s[,'cp.est']*2)
  b.CR = cal.CR(SE = data.s[,'b.se'], theta = b.s, theta.hat = data.s[,'b.est'])
  beta1.CR = cal.CR(SE = data.s[,'beta1.se'], theta = (beta1), theta.hat = data.s[,'beta1.est'])
  beta2.CR = cal.CR(SE = data.s[,'beta2.se'], theta = (beta2), theta.hat = data.s[,'beta2.est'])
  ab.CR.delta = cal.CR(SE = data.s[,'ab.se.delta']*2, theta = (a.s*b.s*2), theta.hat = data.s[,'ab.est']*2)
  
  # convergence rate & postive definite rate
  ConvRate <- 1 - mean(is.na(data[,'a.est']), na.rm=T)
  posDefRate <- mean(data[,'infoDef']==1, na.rm=T)
  
  # saving results
  res<- matrix(NA, 1, 51)
  res[1,] <- c(a.s,cp.s,b.s,beta1,
               bias, Ebias, empSE, seM, seRbias,
               a.CR, cp.CR, b.CR, beta1.CR, beta2.CR, ab.CR.delta,
               a.RR, cp.RR, b.RR, beta1.RR, beta2.RR, ab.RR.delta, ConvRate, posDefRate)
  colnames(res) <- c('a.para', 'cp.para', 'b.para', 'beta1.para',
                     'a.bias','cp.bias', 'b.bias', 'ab.bias','beta1.bias','beta2.bias',
                     'a.Ebias','cp.Ebias', 'b.Ebias', 'ab.Ebias','beta1.Ebias','beta2.Ebias',
                     'a.empSE', 'cp.empSE','b.empSE', 'beta1.empSE', 'beta2.empSE',
                     'a.meanSE', 'cp.meanSE','b.meanSE','beta1.meanSE', 'beta2.meanSE',
                     'a.seRbias', 'cp.seRbias', 'b.seRbias','beta1.seRbias', 'beta2.seRbias',
                     'a.CoverageRate','cp.CoverageRate', 'b.CoverageRate','beta1.CoverageRate','beta2.CoverageRate','ab.CoverageRate', 
                     'a.power', 'a.typeIerror','cp.power', 'cp.typeIerror', 'b.power', 'b.typeIerror', 
                     'beta1.power', 'beta1.typeIerror', 'beta2.power', 'beta2.typeIerror', 'ab.power', 'ab.typeIerror', 
                     'ConvergenceRate','posDefRate')
  return(res)
}






