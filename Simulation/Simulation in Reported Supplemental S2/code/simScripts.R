wd.win = 'C:/Users/officeuser/Desktop/MMA0609/code/'
wd.mac = '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/code/'
setwd(wd.mac)

a.para.all <- c(0,0.3)
b.para.all <- c(0,0.3)
#cp.para.all <- c(0,0.3)
rho.all <- c(0.45, 0.9)
K.all <- c(5, 10, 30)
#beta1.all <- c(0, 0.1)
postV.all <- c(1,1.5)
missingR.all <- c(0.4, 0.6)

total <- length(a.para.all)*length(b.para.all)*length(rho.all)*length(K.all)*length(missingR.all)*length(postV.all)
Cond.list<- matrix(NA, total,7)
colnames(Cond.list) <- c("Cond#",'Rho','a.ps', 'b.s','K','postVar', 'missingRate')

Condi = 0
for(m in missingR.all){
for(p in postV.all){
for(k in K.all){
for(r in rho.all){
for(i in a.para.all){
  for(j in b.para.all){
      Condi = Condi+1
      newname = paste0('_RC',Condi,'.r')
      tx=readLines('TemplateS1S1.R')
      tx0=gsub('PVARii', replace=p, x = tx)
      tx1=gsub('a.COEFi', replace=i, x = tx0)
      tx2=gsub('b.COEFi', replace=j, x = tx1)
      tx4=gsub('RHOii', replace=r, x = tx2)
      tx5=gsub('Kii', replace=k, x = tx4)
      tx6=gsub('MissingRateii', replace=m, x = tx5)
      tx7=gsub('CONDi', replace=Condi, x=tx6)
      writeLines(tx7, con=newname)
      Cond.list[Condi,] = c(paste0('Cond', Condi),r,i/.5,j,k,p,m)
    }}}}}}
win = 'C:/Users/officeuser/Desktop/MMA0609/'
mac = '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/SIMres/'
setwd(mac)
write.csv(Cond.list,'CondList.csv', row.names = F)






















