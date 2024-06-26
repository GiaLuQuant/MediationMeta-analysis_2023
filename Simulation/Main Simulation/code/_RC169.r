library('MASS')
library('metaSEM')
library('fs')
library('xlsx')
library('janitor')

wd.win = 'C:/Users/officeuser/Desktop/MMA0623/'
#wd.mac = '/Users/luzhiming/Desktop/MMA240623/code/'
setwd(wd.win)


source('functions240612.R')

# No. of simulations
nrep=5000

# setting parameters

Tau <- .05
p <- 0.5
Nbar <- 120

K <- 10
Rho <- 0.45
beta1 <- 0.1
hetInd <- 0
postVar <- 1.5

a.s <- 0
b.s<-0
cp.s<-0.3

ab.s <- a.s*b.s
c.s <- cp.s+ab.s

StCoef <- c(a.s, cp.s, b.s)
names(StCoef)=c('a','cp','b')



# generating data
wd.d.win <- 'C:/Users/officeuser/Desktop/MMA0623/MAres/'
#wd.d.mac <- '/Users/luzhiming/Desktop/MMA240623/MAres/'
setwd(wd.d.win)

data <- ma(StCoef = StCoef, K = K, Tau = Tau, Nbar = Nbar, Rho = Rho, 
           nrep = nrep, postVar = postVar, beta1 = beta1, hetInd = hetInd)

for(i in 1:2){
  write.csv(data[[i]], file=paste0(wd.d.win,'MACond', 169,'_M',i, '.csv'))
}


# simulation results
wd.r.win <- 'C:/Users/officeuser/Desktop/MMA0623/SIMres/'
#wd.r.mac <- '/Users/luzhiming/Desktop/MMA240623/SIMres/'
setwd(wd.r.win)

if(hetInd){
   Res1 <- cal_res2(data[[1]], a.s,cp.s,b.s,beta1,beta1)
}else{Res1 <- cal_res1(data[[1]], a.s,cp.s,b.s,beta1)}

write.csv(Res1, file=paste0(wd.r.win,'MMAresCond', 169,'_M',1, '.csv'))
  
N_l = 1000; K_l = 10000; 
postpara = gen_post_heter(a.s, cp.s, b.s, Rho, beta1, N_l, K_l, postVar, hetInd, p=.5)

if(hetInd){
Res2 <- cal_res2(data[[2]], postpara[1] ,postpara[2],postpara[3],postpara[4],postpara[5])
    }else{Res2 <- cal_res1(data[[2]], postpara[1] ,postpara[2],postpara[3],postpara[4])}

write.csv(Res2, file=paste0(wd.r.win,'MMAresCond', 169,'_M',2, '.csv'))
  











