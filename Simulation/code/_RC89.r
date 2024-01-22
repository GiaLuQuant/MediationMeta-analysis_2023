library('MASS')
library('metaSEM')
library('fs')
library('xlsx')
library('janitor')

wd.win = 'C:/Users/officeuser/Desktop/MMA0609/'
wd.mac = '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/code/'
setwd(wd.mac)


source('functions0105.R')

# No. of simulations
nrep=5000

# setting parameters

Tau <- .1
p <- 0.5
Nbar <- 80

K <- 30
Rho <- 0.45
beta1 <- 0.1
postVar <- 1.5

a.s <- 0
b.s<-0
cp.s<-0.3

ab.s <- a.s*b.s
c.s <- cp.s+ab.s

StCoef <- c(a.s, cp.s, b.s)
names(StCoef)=c('a','cp','b')

Missing = TRUE
MissingR = 0.6


# generating data
wd.d.win <- 'C:/Users/officeuser/Desktop/MMA0609/MAres/'
wd.d.mac <- '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/MAres/'
setwd(wd.d.mac)

data <- ma(StCoef, K, Tau, Nbar, Rho, nrep, postVar, beta1, Missing, MissingR)

for(i in 1:2){
  write.csv(data[[i]], file=paste0(wd.d.mac,'MACond', 89,'_M',i, '.csv'))
}


# simulation results
wd.r.win <- 'C:/Users/officeuser/Desktop/MMA0609/SIMres/'
wd.r.mac <- '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/SIMres/'
setwd(wd.r.mac)

Res1 <- cal_res(data[[1]], a.s,cp.s,b.s,beta1)
write.csv(Res1, file=paste0(wd.r.mac,'MMAresCond', 89,'_M',1, '.csv'))
  
N_l = 1000; K_l = 10000; 
postpara = gen_post_heter(a.s, cp.s, b.s, Rho, beta1, N_l, K_l, postVar, p=.5)

Res2 <- cal_res(data[[2]], postpara[1] ,postpara[2],postpara[3],postpara[4])

write.csv(Res2, file=paste0(wd.r.mac,'MMAresCond', 89,'_M',2, '.csv'))
  











