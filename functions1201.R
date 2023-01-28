### code for vcov.osmasem provided by Mike Cheung
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

# funtions for converting values
r2d <- function(R){
  return(2*R/sqrt(1-R^2))
}

d2r <- function(D, n1, n2){
  a = (n1+n2)^2 / (n1*n2)
  return(D/sqrt(D^2+a))
}

dcs2dtb <- function(dcs, rho){
  return(dcs * sqrt(2*(1-rho)))
}

dtb2dcs <- function(dtb, rho){
  return(dtb / sqrt(2*(1-rho)))
}

## functions for the 5 methods being compared in study 1
## Methods being compared in Study 1
# the change-score method
cal.cs.d <- function(E.pre, E.post, C.pre, C.post, n.exp, n.ctl){
  cs.ctl <- C.post-C.pre
  cs.exp <- E.post-E.pre
  md <- mean(cs.exp)-mean(cs.ctl)
  s <- sqrt(((n.exp-1)*var(cs.exp)+(n.ctl-1)*var(cs.ctl))/(n.ctl+n.exp-2))
  return(md/s)
}
# the textbook method
cal.tb.d <- function(E.pre, E.post, C.pre, C.post,n.exp, n.ctl){
  md <- mean(E.post-E.pre)-mean(C.post-C.pre)
  s2.exp <- (var(E.pre)+var(E.post))/2
  s2.ctl <- (var(C.pre)+var(C.post))/2
  s <- sqrt(((n.exp-1)*s2.exp + (n.ctl-1)*s2.ctl)/(n.exp+n.ctl-2))
  return(md/s)
}
# the d-change method (Becker, 1988)
cal.dChange <- function(E.pre, E.post, C.pre, C.post){
  d.e = mean(E.post-E.pre) / sd(E.pre)
  d.c = mean(C.post-C.pre) / sd(C.pre)
  return(d.e-d.c)
}
# Cohen's d with pooled pretest sd (recommended by Morris in 2008)
cal.pre.d <- function(E.pre, E.post, C.pre, C.post, n.exp, n.ctl){
  md = mean(E.post-E.pre)- mean(C.post-C.pre)
  s = sqrt(((n.ctl-1)*var(C.pre)+(n.exp-1)*var(E.pre))/(n.ctl+n.exp-2))
  return(md/s)
}

## Methods being compared in Study 2
# the posttest-score method
cal.post.d <- function(E.post, C.post, n.exp, n.ctl){
  md <- mean(E.post)-mean(C.post)
  s <- sqrt(((n.ctl-1)*var(C.post)+(n.exp-1)*var(E.post))/(n.ctl+n.exp-2))
  return(md/s)
}
# the textbook method
cal.tb.d.PS <- function(E.pre, E.post, C.pre, C.post,n.exp, n.ctl){
  md <- mean(E.post)-mean(C.post)
  s2.exp <- (var(E.pre)+var(E.post))/2
  s2.ctl <- (var(C.pre)+var(C.post))/2
  s <- sqrt(((n.exp-1)*s2.exp + (n.ctl-1)*s2.ctl)/(n.exp+n.ctl-2))
  return(md/s)
}
# Cohen's d with pooled pretest sd (recommended by Morris in 2008)
cal.pre.d.PS <- function(E.pre, E.post, C.pre, C.post, n.exp, n.ctl){
  md = mean(E.post)- mean(C.post)
  s = sqrt(((n.ctl-1)*var(C.pre)+(n.exp-1)*var(E.pre))/(n.ctl+n.exp-2))
  return(md/s)
}

# generating a new matrix
Gen.matrix <- function(){
  cm.ps=matrix(NA, 3,3)
  colnames(cm.ps) = rownames(cm.ps) = c('X', 'M', 'Y')
  diag(cm.ps) = 1
  return(cm.ps)
}

# generating parameters with heterogeneity
Gen.para <- function(Parameters, K, Tau){
  # the i_th row of the matrix contains the change-score a/b/c/cp for the Ki study in the MA
  Para.H <- matrix(NA, K, 3)
  colnames(Para.H) <- c('r.XM','r.MY','r.XY')
  for (i in 1:3){
    p <- rnorm(K, Parameters[i], Tau)
    Para.H[, i] <- p
  }
  return(Para.H)
}

genN <- function(Nbar, K){
  if (Nbar==50){
    Nsd=18
    size = (0.2*Nbar^2)/(Nsd^2-.1*Nbar)
    N = rnbinom(K,size=size, mu=0.6*Nbar) + 0.4*Nbar
  }
  if (Nbar==80){
    Nsd=23
    size = (0.1*Nbar^2)/(Nsd^2-.2*Nbar)
    N = rnbinom(K,size=size, mu=0.6*Nbar) + 0.4*Nbar
  }
  if (Nbar==150){
    Nsd=30
    size = (0.1*Nbar^2)/(Nsd^2-.1*Nbar)
    N = rnbinom(K,size=size, mu=0.7*Nbar) + 0.3*Nbar
  }
  return(N)
}

data_gen <- function(a.s,cp.s, b.s, c.s,Rho,X,n.exp,n.ctl,p=.5){
  # value converting
  d.cs.M = r2d(a.s)
  d.cs.Y = r2d(c.s)
  d.tb.M = dcs2dtb(d.cs.M, Rho)
  d.tb.Y = dcs2dtb(d.cs.Y, Rho)
  
  # variances for variables
  sigma.X2 = p*(1-p)
  sigma.Mg2 = sigma.Yg2 = 2*(1-Rho)
  sigma.Mcs2 = sigma.Mg2 + d.tb.M^2/4
  sigma.Ycs2 = sigma.Yg2 + d.tb.Y^2/4
  
  # unstandardized coefficients
  a.u = a.s* sqrt(sigma.Mcs2/sigma.X2)
  b.u = b.s* sqrt(sigma.Ycs2/sigma.Mcs2)
  cp.u = cp.s* sqrt(sigma.Ycs2/sigma.X2)
  c.u <- a.u*b.u + cp.u

  # generating M
  Rho.pp = matrix(c(1,Rho,Rho,1),2,2)
  mu.M.exp <- c(0, d.tb.M); mu.M.ctl <- c(0,0)
  M.exp <- mvrnorm(n.exp, mu.M.exp, Rho.pp)
  M.ctl <- mvrnorm(n.ctl, mu.M.ctl, Rho.pp)
  M.T.cs <- M.exp[,2] - M.exp[,1]
  M.C.cs <- M.ctl[,2] - M.ctl[,1]
  M.cs <- c(M.C.cs, M.T.cs)
  M.post <- c(M.ctl[,2], M.exp[,2])
  M.pre <- c(M.ctl[,1], M.exp[,1])
  
  # Genrating Y
  sigma.eYcs2T =sigma.eYcs2C= sigma.Yg2 - b.u^2*sigma.Mg2 
  e.Ycs.T = rnorm(n.exp, 0, sqrt(sigma.eYcs2T))
  e.Ycs.C = rnorm(n.ctl, 0, sqrt(sigma.eYcs2C))
  e.Ypre.T = rnorm(n.exp,0,sqrt((1+Rho)/2))
  e.Ypre.C = rnorm(n.ctl,0,sqrt((1+Rho)/2))
  Y.T.cs = d.tb.Y-b.u*d.tb.M+b.u*M.T.cs+e.Ycs.T
  Y.C.cs =b.u*M.C.cs+e.Ycs.C
  Y.cs = c(Y.C.cs, Y.T.cs)
  Y.T.pre = d.tb.Y/2 - Y.T.cs/2 + e.Ypre.T
  Y.C.pre = -Y.C.cs/2 + e.Ypre.C
  Y.pre = c(Y.C.pre, Y.T.pre)
  Y.post = Y.pre+Y.cs 
  
  dat <- list(X=X,Mcs = M.cs, Ycs = Y.cs, M.post = M.post, Y.post = Y.post,
              M.T.pre = M.exp[,1], M.T.post = M.exp[,2], 
              M.C.pre = M.ctl[,1], M.C.post = M.ctl[,2], 
              Y.T.pre = Y.T.pre, Y.T.post = Y.post[X==1], 
              Y.C.pre = Y.C.pre, Y.C.post = Y.post[X==0])
  return(dat)
}

data_gen_HV <- function(a.s,cp.s, b.s, c.s,Rho,X,n.exp,n.ctl,p=.5){
  # value converting
  d.cs.M = r2d(a.s)
  d.cs.Y = r2d(c.s)
  d.tb.M = d.cs.M*sqrt(2.25-2*Rho)
  d.tb.Y = d.cs.Y*sqrt(2.25-2*Rho)
  
  # variances for variables
  sigma.X2 = p*(1-p)
  sigma.McsT2 = sigma.eMcsT2 =  2.5-2*Rho
  sigma.McsC2 = sigma.eMcsC2 = 2-2*Rho
  sigma.Mcs2 = (sigma.McsT2+sigma.McsC2)/2 + (1.125*d.tb.M)^2/4
  
  sigma.YcsT2 = 2.5-2*Rho
  sigma.YcsC2 = 2-2*Rho
  sigma.Ycs2 = (sigma.YcsT2+sigma.YcsC2)/2 + (1.125*d.tb.Y)^2/4
  
  # unstandardized coefficients
  a.u = a.s* sqrt(sigma.Mcs2/sigma.X2)
  b.u = b.s* sqrt(sigma.Ycs2/sigma.Mcs2)
  cp.u = cp.s* sqrt(sigma.Ycs2/sigma.X2)
  c.u <- a.u*b.u + cp.u
  
  # generating M
  Rho.pp.C = matrix(c(1,Rho,Rho,1),2,2)
  Rho.pp.T = matrix(c(1,Rho,Rho,1.5),2,2)
  
  mu.M.exp <- c(0, (1.125*d.tb.M)); mu.M.ctl <- c(0,0)
  M.exp <- mvrnorm(n.exp, mu.M.exp, Rho.pp.T)
  M.ctl <- mvrnorm(n.ctl, mu.M.ctl, Rho.pp.C)
  M.T.cs <- M.exp[,2] - M.exp[,1]
  M.C.cs <- M.ctl[,2] - M.ctl[,1]
  M.cs <- c(M.C.cs, M.T.cs)
  M.post <- c(M.ctl[,2], M.exp[,2])
  M.pre <- c(M.ctl[,1], M.exp[,1])
  
  # Genrating Y
  sigma.eYcsT2 = sigma.YcsT2 - b.u^2*sigma.McsT2
  sigma.eYcsC2 = sigma.YcsC2 - b.u^2*sigma.McsC2

  e.Ycs.T = rnorm(n.exp, 0, sqrt(sigma.eYcsT2))
  e.Ycs.C = rnorm(n.ctl, 0, sqrt(sigma.eYcsC2))
  
  sigma.eYpreT2 = 1 - (Rho-1)^2/(2.5-2*Rho)
  e.Ypre.T = rnorm(n.exp,0,sqrt(sigma.eYpreT2))
  e.Ypre.C = rnorm(n.ctl,0,sqrt((1+Rho)/2))
  
  iY.T.cs = 1.125*d.tb.Y - 1.125*b.u*d.tb.M 
  Y.T.cs = iY.T.cs + b.u*M.T.cs + e.Ycs.T
  Y.C.cs = b.u*M.C.cs + e.Ycs.C
  Y.cs = c(Y.C.cs, Y.T.cs)
  
  Y.T.pre = ((Rho-1)/(2.5-2*Rho))*(Y.T.cs - 1.125*d.tb.Y) + e.Ypre.T
  Y.C.pre = -Y.C.cs/2 + e.Ypre.C
  Y.pre = c(Y.C.pre, Y.T.pre)
  Y.post = Y.pre+Y.cs 
  
  dat <- list(X=X,Mcs = M.cs, Ycs = Y.cs, M.post = M.post, Y.post = Y.post,
              M.T.pre = M.exp[,1], M.T.post = M.exp[,2], 
              M.C.pre = M.ctl[,1], M.C.post = M.ctl[,2], 
              Y.T.pre = Y.T.pre, Y.T.post = Y.post[X==1], 
              Y.C.pre = Y.C.pre, Y.C.post = Y.post[X==0])
  return(dat)
}

gen_post <- function(a.s, cp.s, b.s, Rho, N, HV=1,p=.5){
  c.s = a.s*b.s+cp.s
  n.exp=n.ctl=N/2
  X = c(rep(0,n.ctl), rep(1,n.exp))
  if(HV==1){
    dat = data_gen(a.s, cp.s, b.s, c.s, Rho, X, n.exp, n.ctl, p=.5)
  }else{dat = data_gen_HV(a.s=a.s, cp.s=cp.s,b.s=b.s,
                    c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)}
  
  eq1 = lm(dat$M.post~X)
  eq2 = lm(dat$Y.post~dat$M.post+X)
  
  a.ps = summary(eq1)$coefficients['X','Estimate']*(sd(X)/sd(dat$M.post))
  cp.ps = summary(eq2)$coefficients['X','Estimate']*(sd(X)/sd(dat$Y.post))
  b.ps = summary(eq2)$coefficients['dat$M.post','Estimate']*(sd(dat$M.post)/sd(dat$Y.post))
  return(c(a.ps, cp.ps, b.ps))
}

dg.ps <- function(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=.5){
  # generating true values with heterogeneity for individual studies
  c.s = StCoef['a']*StCoef['b']+StCoef['cp']
  r.XY = c.s; r.XM = StCoef['a']; r.MY = StCoef['a']*StCoef['cp']+StCoef['b']
  para.h = Gen.para(c(r.XM,r.MY,r.XY), K, Tau)
  # generate sample size for individual studies
  N = genN(Nbar, K)
  data4MA <- list(TB=list(), CS=list(),Becker=list(), Morris=list(),N = c())
  data4MA[['N']] <- N
  
  for(ki in 1:K){
    # converting r.i into path coefficients
    a.s = para.h[ki,'r.XM'];c.s=para.h[ki,'r.XY']
    b.s = (para.h[ki,'r.MY']-a.s*c.s)/(1-a.s^2)
    cp.s = c.s-a.s*b.s
    
    X = sort(rbinom(N[ki], 1, p))
    n.exp=sum(X==1); n.ctl=sum(X==0)
    if(hetero_post_var==1){
      dat = data_gen(a.s=a.s, cp.s=cp.s,b.s=b.s,
                     c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }else{
      dat = data_gen_HV(a.s=a.s, cp.s=cp.s,b.s=b.s,
                     c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }
    for(i in 1:4){
      data4MA[[i]][[ki]] <- Gen.matrix()
      data4MA[[i]][[ki]][2,3]=data4MA[[i]][[ki]][3,2] = cor(dat$Mcs, dat$Ycs)
    }

    ## calculating effect sizes using the 5 methods
    # TB method 
    d.tb.M <- cal.tb.d(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.tb.Y <- cal.tb.d(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.tb = d2r(d.tb.M, n.exp, n.ctl)
    rXY.tb = d2r(d.tb.Y, n.exp, n.ctl)
    data4MA[['TB']][[ki]][1,2]=data4MA[['TB']][[ki]][2,1]=rXM.tb
    data4MA[['TB']][[ki]][1,3]=data4MA[['TB']][[ki]][3,1]=rXY.tb
    
    # CS method
    d.cs.M = cal.cs.d(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.cs.Y <- cal.cs.d(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.cs = d2r(d.cs.M, n.exp, n.ctl)
    rXY.cs = d2r(d.cs.Y, n.exp, n.ctl)
    data4MA[['CS']][[ki]][1,2]=data4MA[['CS']][[ki]][2,1]=rXM.cs
    data4MA[['CS']][[ki]][1,3]=data4MA[['CS']][[ki]][3,1]=rXY.cs
    
    # Becker method
    d.B.M <- cal.dChange(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post)
    d.B.Y <- cal.dChange(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post)
    rXM.B = d2r(d.B.M, n.exp, n.ctl)
    rXY.B = d2r(d.B.Y, n.exp, n.ctl)
    data4MA[['Becker']][[ki]][1,2]=data4MA[['Becker']][[ki]][2,1]=rXM.B
    data4MA[['Becker']][[ki]][1,3]=data4MA[['Becker']][[ki]][3,1]=rXY.B
    
    # Morris method
    d.mo.M <- cal.pre.d(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.mo.Y <- cal.pre.d(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.mo = d2r(d.mo.M, n.exp, n.ctl)
    rXY.mo = d2r(d.mo.Y, n.exp, n.ctl)
    data4MA[['Morris']][[ki]][1,2]=data4MA[['Morris']][[ki]][2,1]=rXM.mo
    data4MA[['Morris']][[ki]][1,3]=data4MA[['Morris']][[ki]][3,1]=rXY.mo
  }
  return(data4MA)
}

dg.ps2 <- function(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=.5){
  # generating true values with heterogeneity for individual studies
  c.s = StCoef['a']*StCoef['b']+StCoef['cp']
  r.XY = c.s; r.XM = StCoef['a']; r.MY = StCoef['a']*StCoef['cp']+StCoef['b']
  para.h = Gen.para(c(r.XM,r.MY,r.XY), K, Tau)
  # generate sample size for individual studies
  N = genN(Nbar, K)
  data4MA <- list(TB=list(), PS=list(), Morris=list(),N = c())
  data4MA[['N']] <- N
  for(ki in 1:K){
    # converting r.i into path coefficients
    a.s = para.h[ki,'r.XM'];c.s=para.h[ki,'r.XY']
    b.s = (para.h[ki,'r.MY']-a.s*c.s)/(1-a.s^2)
    cp.s = c.s-a.s*b.s
    
    X = sort(rbinom(N[ki], 1, p))
    n.exp=sum(X==1); n.ctl=sum(X==0)
    if(hetero_post_var==1){
      dat = data_gen(a.s=a.s, cp.s=cp.s,b.s=b.s,
                     c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }else{
      dat = data_gen_HV(a.s=a.s, cp.s=cp.s,b.s=b.s,
                        c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }
    for(i in 1:3){
      data4MA[[i]][[ki]] <- Gen.matrix()
      data4MA[[i]][[ki]][2,3]=data4MA[[i]][[ki]][3,2] = cor(dat$M.post, dat$Y.post)
      
    }

    ## calculating effect sizes using the 3 methods
    # TB method 
    d.tb.M <- cal.tb.d.PS(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.tb.Y <- cal.tb.d.PS(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.tb = d2r(d.tb.M, n.exp, n.ctl)
    rXY.tb = d2r(d.tb.Y, n.exp, n.ctl)
    data4MA[['TB']][[ki]][1,2]=data4MA[['TB']][[ki]][2,1]=rXM.tb
    data4MA[['TB']][[ki]][1,3]=data4MA[['TB']][[ki]][3,1]=rXY.tb
    
    # PS method
    d.ps.M <- cal.post.d(dat$M.T.post, dat$M.C.post, n.exp, n.ctl)
    d.ps.Y <- cal.post.d(dat$Y.T.post, dat$Y.C.post, n.exp, n.ctl)
    rXM.ps = d2r(d.ps.M, n.exp, n.ctl)
    rXY.ps = d2r(d.ps.Y, n.exp, n.ctl)
    data4MA[['PS']][[ki]][1,2]=data4MA[['PS']][[ki]][2,1]=rXM.ps
    data4MA[['PS']][[ki]][1,3]=data4MA[['PS']][[ki]][3,1]=rXY.ps
    
    # Morris method
    d.mo.M <- cal.pre.d.PS(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.mo.Y <- cal.pre.d.PS(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.mo = d2r(d.mo.M, n.exp, n.ctl)
    rXY.mo = d2r(d.mo.Y, n.exp, n.ctl)
    data4MA[['Morris']][[ki]][1,2]=data4MA[['Morris']][[ki]][2,1]=rXM.mo
    data4MA[['Morris']][[ki]][1,3]=data4MA[['Morris']][[ki]][3,1]=rXY.mo
  }
  return(data4MA)
}

dg.ps3 <- function(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=.5){
  # generating true values with heterogeneity for individual studies
  c.s = StCoef['a']*StCoef['b']+StCoef['cp']
  r.XY = c.s; r.XM = StCoef['a']; r.MY = StCoef['a']*StCoef['cp']+StCoef['b']
  para.h = Gen.para(c(r.XM,r.MY,r.XY), K, Tau)
  # generate sample size for individual studies
  N = genN(Nbar, K)
  data4MA <- list(CS=list(),PS=list(), Morris=list(),N = c())
  data4MA[['N']] <- N
  
  for(ki in 1:K){
    # converting r.i into path coefficients
    a.s = para.h[ki,'r.XM'];c.s=para.h[ki,'r.XY']
    b.s = (para.h[ki,'r.MY']-a.s*c.s)/(1-a.s^2)
    cp.s = c.s-a.s*b.s
    
    X = sort(rbinom(N[ki], 1, p))
    n.exp=sum(X==1); n.ctl=sum(X==0)
    if(hetero_post_var==1){
      dat = data_gen(a.s=a.s, cp.s=cp.s,b.s=b.s,
                     c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }else{
      dat = data_gen_HV(a.s=a.s, cp.s=cp.s,b.s=b.s,
                        c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }
    
    for(i in 1:3){
      data4MA[[i]][[ki]] <- Gen.matrix()
    }

    ## calculating effect sizes using the 5 methods
    # CS method
    d.cs.M = cal.cs.d(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.cs.Y <- cal.cs.d(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.cs = d2r(d.cs.M, n.exp, n.ctl)
    rXY.cs = d2r(d.cs.Y, n.exp, n.ctl)
    data4MA[['CS']][[ki]][1,2]=data4MA[['CS']][[ki]][2,1]=rXM.cs
    data4MA[['CS']][[ki]][1,3]=data4MA[['CS']][[ki]][3,1]=rXY.cs
    data4MA[['CS']][[ki]][2,3]=data4MA[['CS']][[ki]][3,2]=cor(dat$Mcs, dat$Ycs)
    # PS method
    d.ps.M <- cal.post.d(dat$M.T.post, dat$M.C.post, n.exp, n.ctl)
    d.ps.Y <- cal.post.d(dat$Y.T.post, dat$Y.C.post, n.exp, n.ctl)
    rXM.ps = d2r(d.ps.M, n.exp, n.ctl)
    rXY.ps = d2r(d.ps.Y, n.exp, n.ctl)
    data4MA[['PS']][[ki]][1,2]=data4MA[['PS']][[ki]][2,1]=rXM.ps
    data4MA[['PS']][[ki]][1,3]=data4MA[['PS']][[ki]][3,1]=rXY.ps
    data4MA[['PS']][[ki]][2,3]=data4MA[['PS']][[ki]][3,2]=cor(dat$M.post, dat$Y.post)

    # Morris method
    d.mo.M <- cal.pre.d.PS(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.mo.Y <- cal.pre.d.PS(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.mo = d2r(d.mo.M, n.exp, n.ctl)
    rXY.mo = d2r(d.mo.Y, n.exp, n.ctl)
    data4MA[['Morris']][[ki]][1,2]=data4MA[['Morris']][[ki]][2,1]=rXM.mo
    data4MA[['Morris']][[ki]][1,3]=data4MA[['Morris']][[ki]][3,1]=rXY.mo
    data4MA[['Morris']][[ki]][2,3]=data4MA[['Morris']][[ki]][3,2]=cor(dat$M.post, dat$Y.post)
  }
  return(data4MA)
}


dg.ps4 <- function(StCoef, K, Tau, Nbar, Rho, hetero_post_var, MCrate,missingPos, p=.5){
  # generating true values with heterogeneity for individual studies
  c.s = StCoef['a']*StCoef['b']+StCoef['cp']
  r.XY = c.s; r.XM = StCoef['a']; r.MY = StCoef['a']*StCoef['cp']+StCoef['b']
  para.h = Gen.para(c(r.XM,r.MY,r.XY), K, Tau)
  # generate sample size for individual studies
  N = genN(Nbar, K)
  data4MA <- list(CS=list(), PS=list(),Morris=list(),N = c())
  data4MA[['N']] <- N
  # generate index for studies WITH missing correlation. 
  # 0 represents NO missing correlation introduced
  num.missing=MCrate*K; num.NM = K-num.missing
  missingCor_index=rep(c(0,1),times=c(K-MCrate*K, MCrate*K))
  missingCor_index=sample(missingCor_index)
  for(ki in 1:K){
    # converting r.i into path coefficients
    a.s = para.h[ki,'r.XM'];c.s=para.h[ki,'r.XY']
    b.s = (para.h[ki,'r.MY']-a.s*c.s)/(1-a.s^2)
    cp.s = c.s-a.s*b.s
    
    X = sort(rbinom(N[ki], 1, p))
    n.exp=sum(X==1); n.ctl=sum(X==0)
    
    if(hetero_post_var==1){
      dat = data_gen(a.s=a.s, cp.s=cp.s,b.s=b.s,
                     c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }else{
      dat = data_gen_HV(a.s=a.s, cp.s=cp.s,b.s=b.s,
                        c.s=c.s,Rho=Rho, X,n.exp=n.exp, n.ctl=n.ctl)
    }
    
    for(i in 1:3){
      data4MA[[i]][[ki]] <- Gen.matrix()
    }
    
    # CS method
    d.cs.M = cal.cs.d(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.cs.Y <- cal.cs.d(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.cs = d2r(d.cs.M, n.exp, n.ctl)
    rXY.cs = d2r(d.cs.Y, n.exp, n.ctl)
    data4MA[['CS']][[ki]][1,2]=data4MA[['CS']][[ki]][2,1]=rXM.cs
    data4MA[['CS']][[ki]][1,3]=data4MA[['CS']][[ki]][3,1]=rXY.cs
    data4MA[['CS']][[ki]][2,3]=data4MA[['CS']][[ki]][3,2] = cor(dat$Mcs, dat$Ycs)
    
    # PS method
    d.ps.M <- cal.post.d(dat$M.T.post, dat$M.C.post, n.exp, n.ctl)
    d.ps.Y <- cal.post.d(dat$Y.T.post, dat$Y.C.post, n.exp, n.ctl)
    rXM.ps = d2r(d.ps.M, n.exp, n.ctl)
    rXY.ps = d2r(d.ps.Y, n.exp, n.ctl)
    data4MA[['PS']][[ki]][1,2]=data4MA[['PS']][[ki]][2,1]=rXM.ps
    data4MA[['PS']][[ki]][1,3]=data4MA[['PS']][[ki]][3,1]=rXY.ps
    data4MA[['PS']][[ki]][2,3]=data4MA[['PS']][[ki]][3,2] = cor(dat$M.post, dat$Y.post)
    
    # Morris method
    d.mo.M <- cal.pre.d.PS(dat$M.T.pre,dat$M.T.post,dat$M.C.pre, dat$M.C.post, n.exp, n.ctl)
    d.mo.Y <- cal.pre.d.PS(dat$Y.T.pre,dat$Y.T.post,dat$Y.C.pre, dat$Y.C.post, n.exp, n.ctl)
    rXM.mo = d2r(d.mo.M, n.exp, n.ctl)
    rXY.mo = d2r(d.mo.Y, n.exp, n.ctl)
    data4MA[['Morris']][[ki]][1,2]=data4MA[['Morris']][[ki]][2,1]=rXM.mo
    data4MA[['Morris']][[ki]][1,3]=data4MA[['Morris']][[ki]][3,1]=rXY.mo
    data4MA[['Morris']][[ki]][2,3]=data4MA[['Morris']][[ki]][3,2] = cor(dat$M.post, dat$Y.post)
    
    if(missingPos == '1'){
      if(missingCor_index[ki]==1){
        for(i in 1:3){
          data4MA[[i]][[ki]][2,3]=data4MA[[i]][[ki]][3,2] = NA
        }
      }
    }
    if(missingPos=='0'){
      if(missingCor_index[ki]==1){
        for(i in 1:3){
          data4MA[[i]][[ki]][1,2]=data4MA[[i]][[ki]][2,1] = NA
          data4MA[[i]][[ki]][1,3]=data4MA[[i]][[ki]][3,1] = NA
        }
      }
    }
  }
  return(data4MA)
}


extMA <- function(Fit){
  s = summary(Fit)
  Est = s$parameters[1:3,'Estimate']  # a cp b
  
  c.est = Est[1]*Est[3]+Est[2]
  ab.est = Est[1]*Est[3]
  
  se = s$parameters[1:3,'Std.Error'] 
  pValue = s$parameters[1:3,'Pr(>|z|)']
  
  cov.ab <- vcov.osmasem(Fit)['b','a']
  cov.a <- vcov.osmasem(Fit)['a','a']
  cov.b <- vcov.osmasem(Fit)['b','b']
  se.ab.delta <- sqrt(Est[3]^2*cov.a + Est[1]^2*cov.b + 2*Est[1]*Est[3]*cov.ab)
  
  infoDef = 1
  
  res = c(Est, c.est, ab.est, se, se.ab.delta, 
          cov.a, cov.b, cov.ab, pValue, infoDef)
  return(res)
}


ma1 <- function(StCoef, K, Tau, Nbar, Rho, nrep, hetero_post_var=1){
  RES = list()
  res.i = matrix(NA, nrep, 16)
  colnames(res.i) <- c('a.est', 'cp.est', 'b.est', 'c.est','ab.est',
                       'a.se','cp.se','b.se','ab.se.delta',
                       'cov.a','cov.b','cov.ab','a.p','cp.p','b.p','infoDef')
  for(i in 1:4){ 
    RES[[i]] = res.i
  };names(RES) <- c('TB','CS','Becker','Morris')

  # model
  MedM <- 'M ~ a*X
           Y ~ cp*X + b*M'
  
  RAM1 <- lavaan2RAM(MedM, obs.variables=c("X", "M", "Y"))
  RAM1$S[1,1] <- 1
  M0 <- create.vechsR(A0=RAM1$A, S0=RAM1$S)
  T0 <- create.Tau2(RAM=RAM1, RE.type='Diag', Transform="expLog", RE.startvalues=0.05)

  for (i in 1:nrep){
    dat = try(dg.ps(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=.5))
    if(inherits(dat, 'try-error')){
      RES[[1]][i,]=RES[[2]][i,]=RES[[3]][i,]=RES[[4]][i,]=rep(NA, 16)
    }else{
      for(j in 1:4){
        df <- Cor2DataFrame(dat[[j]], dat$N)
        fit0 = try(osmasem(model.name="mediationMA", Mmatrix=M0, Tmatrix=T0, data=df))
        if(inherits(fit0, 'try-error')){
          res <- c(rep(NA, 16))
        }else{
          if(summary(fit0)$infoDefinite==FALSE){
            res <- extMA(fit0) 
            res[16] <- -1
          }else{res<- extMA(fit0)}
        }
        RES[[j]][i,]=res
      }
    }
  }
  return(RES)
}

ma2 <- function(StCoef, K, Tau, Nbar, Rho, nrep, hetero_post_var=1){
  RES = list()
  res.i = matrix(NA, nrep, 16)
  colnames(res.i) <- c('a.est', 'cp.est', 'b.est', 'c.est','ab.est',
                       'a.se','cp.se','b.se','ab.se.delta',
                       'cov.a','cov.b','cov.ab','a.p','cp.p','b.p','infoDef')
  for(i in 1:3){ 
    RES[[i]] = res.i
  };names(RES) <- c('TB','PS','Morris')
  
  # model
  MedM <- 'M ~ a*X
           Y ~ cp*X + b*M'
  
  RAM1 <- lavaan2RAM(MedM, obs.variables=c("X", "M", "Y"))
  RAM1$S[1,1] <- 1
  M0 <- create.vechsR(A0=RAM1$A, S0=RAM1$S)
  T0 <- create.Tau2(RAM=RAM1, RE.type='Diag', Transform="expLog", RE.startvalues=0.05)
  
  for (i in 1:nrep){
    dat = try(dg.ps2(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=0.5))
    if(inherits(dat, 'try-error')){
      RES[[1]][i,]=RES[[2]][i,]=RES[[3]][i,]=rep(NA, 16)
    }else{
      for(j in 1:3){
        df <- Cor2DataFrame(dat[[j]], dat$N)
        fit0 = try(osmasem(model.name="mediationMA", Mmatrix=M0, Tmatrix=T0, data=df))
        if(inherits(fit0, 'try-error')){
          res <- c(rep(NA, 16))
        }else{
          if(summary(fit0)$infoDefinite==FALSE){
            res <- extMA(fit0) 
            res[16] <- -1
          }else{res<- extMA(fit0)}
        }
        RES[[j]][i,]=res
      }
    }
  }
  return(RES)
}

ma34 <- function(StCoef, K, Tau, Nbar, Rho, nrep,MCrate=0,missingPos=NULL, hetero_post_var=1){
  RES = list()
  res.i = matrix(NA, nrep, 16)
  colnames(res.i) <- c('a.est', 'cp.est', 'b.est', 'c.est','ab.est',
                       'a.se','cp.se','b.se','ab.se.delta',
                       'cov.a','cov.b','cov.ab','a.p','cp.p','b.p','infoDef')
  for(i in 1:3){ 
    RES[[i]] = res.i
  };names(RES) <- c('CS','PS','Morris')
  
  # model
  MedM <- 'M ~ a*X
           Y ~ cp*X + b*M'
  
  RAM1 <- lavaan2RAM(MedM, obs.variables=c("X", "M", "Y"))
  RAM1$S[1,1] <- 1
  M0 <- create.vechsR(A0=RAM1$A, S0=RAM1$S)
  T0 <- create.Tau2(RAM=RAM1, RE.type='Diag', Transform="expLog", RE.startvalues=0.05)
  
  for (i in 1:nrep){
    if(MCrate==0){
      dat = try(dg.ps3(StCoef, K, Tau, Nbar, Rho, hetero_post_var, p=0.5))
    }else{dat = try(dg.ps4(StCoef, K, Tau, Nbar, Rho, hetero_post_var,MCrate,missingPos, p=0.5))}

    if(inherits(dat, 'try-error')){
      RES[[1]][i,]=RES[[2]][i,]=RES[[3]][i,]=rep(NA, 16)
    }else{
      for(j in 1:3){
        df <- Cor2DataFrame(dat[[j]], dat$N)
        fit0 = try(osmasem(model.name="mediationMA", Mmatrix=M0, Tmatrix=T0, data=df))
        if(inherits(fit0, 'try-error')){
          res <- c(rep(NA, 16))
        }else{
          if(summary(fit0)$infoDefinite==FALSE){
            res <- extMA(fit0) 
            res[16] <- -1
          }else{res<- extMA(fit0)}
        }
        RES[[j]][i,]=res
      }
    }
  }
  return(RES)
}



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
    Ebias <- mean(theta.hat)-theta
  }else{Ebias <- (mean(theta.hat)-theta)/theta}
  return(Ebias)
}

# coverage rate
cal.CR <- function(SE, theta, theta.hat){
  cv.error = qnorm(1-0.05/2) * SE
  l = theta.hat-cv.error  # lower bound
  u = theta.hat+cv.error  # upper bound
  lcriterion = (theta >= l)
  ucriterion = (theta <= u)
  CR = mean(lcriterion*ucriterion)
  return(CR)
}

cal_res <- function(data,a.s,cp.s,b.s){
  p=.5
  sigma.X = sqrt(p*(1-p))
  # deleting negative-definite repetitions
  data.s <- select_res(data)
  
  # calculating EBIAS
  bias = c(mean(data.s[,'a.est']-a.s)/sigma.X, 
           mean(data.s[,'b.est']-b.s),
           mean(data.s[,'ab.est']-a.s*b.s)/sigma.X)
  Ebias = c(cal.eBias(data.s[,'a.est']/sigma.X, a.s/sigma.X), 
            cal.eBias(data.s[,'b.est'],b.s),
            cal.eBias(data.s[,'ab.est']/sigma.X, a.s*b.s/sigma.X))
  
  # calculating empirical SEs of a, b, and cp
  empSE <- apply(data.s[,1:3],2,sd); names(empSE)=c('a.eSE','cp.eSE','b.eSE')
  seM <- apply(data.s[,6:8],2,mean)
  seRbias <- c(cal.eBias(data.s[,'a.se'],empSE['a.eSE']), 
              cal.eBias(data.s[,'cp.se'],empSE['cp.eSE']),
              cal.eBias(data.s[,'b.se'],empSE['b.eSE']))

  # rejection rate (for a only, b only, and a*b)
  a.RR <- mean(data.s[,'a.p']<=0.05)
  b.RR <- mean(data.s[,'b.p']<=0.05)
  
  Z.delta <- data.s[,'ab.est'] / data.s[,'ab.se.delta']
  ab.p.delta <- 1-pnorm(abs(Z.delta))
  ab.RR.delta <- mean(ab.p.delta<=0.025)
  
  # coverage rate
  a.CR = cal.CR(SE = sqrt(data.s[,'cov.a']), theta = a.s, theta.hat = data.s[,'a.est'])
  b.CR = cal.CR(SE = sqrt(data.s[,'cov.b']), theta = b.s, theta.hat = data.s[,'b.est'])
  ab.CR.delta = cal.CR(SE = data.s[,'ab.se.delta'], theta = (a.s*b.s), theta.hat = data.s[,'ab.est'])

  # convergence rate & postive definite rate
  ConvRate <- 1 - mean(is.na(data[,'a.est']))
  posDefRate <- mean(data[,'infoDef']==1)
  
  # saving results
  res<- matrix(NA, 1, 23)
  res[1,] <- c(Ebias, bias, empSE, seM, seRbias,
               a.CR, b.CR, ab.CR.delta,
               a.RR, b.RR, ab.RR.delta, ConvRate, posDefRate)
  colnames(res) <- c('a.Ebias', 'b.Ebias', 'ab.Ebias',
                     'a.bias','b.bias','ab.bias',
                     'a.empSE', 'cp.empSE','b.empSE',
                     'a.meanSE', 'cp.meanSE','b.meanSE',
                     'a.seRbias', 'cp.seRbias', 'b.seRbias',
                     'a.CoverageRate', 'b.CoverageRate','ab.CoverageRate', 
                     'a.RejectionRate', 'b.RejectionRate','ab.RejectionRate',
                     'ConvergenceRate','posDefRate')
  return(res)
}







