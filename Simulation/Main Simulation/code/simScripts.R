wd.win = 'C:/Users/officeuser/Desktop/MMA0623/code/'
# wd.mac = '/Users/luzhiming/Desktop/MMA240623/code/'
setwd(wd.win)

a.para.all <- c(0,0.3)
b.para.all <- c(0,0.3)
HetInd <- c(1, 0)
rho.all <- c(0.45, 0.9)
K.all <- c(5, 10, 30)
Nbar.all <- c(80, 120)
postV.all <- c(1,1.5)

total <- length(a.para.all)*length(b.para.all)*length(HetInd)*length(rho.all)*length(K.all)*length(Nbar.all)*length(postV.all)
Cond.list<- matrix(NA, total,8)
colnames(Cond.list) <- c("Cond#",'Rho','HetInd','a.ps', 'b.s','K','Nbar','postVar')

Condi = 0

for(h in HetInd){
for(p in postV.all){
for(k in K.all){
for(n in Nbar.all){
for(r in rho.all){
for(i in a.para.all){
for(j in b.para.all){
  Condi = Condi+1
  newname = paste0('_RC',Condi,'.r')
  tx=readLines('TemplateMMA.R')
  tx0=gsub('PVARii', replace=p, x = tx)
  tx1=gsub('a.COEFi', replace=i, x = tx0)
  tx2=gsub('b.COEFi', replace=j, x = tx1)
  tx4=gsub('RHOii', replace=r, x = tx2)
  tx5=gsub('Kii', replace=k, x = tx4)
  tx7=gsub('Nbarii', replace=n, x = tx5)
  tx8=gsub('HetIndii', replace=h, x = tx7)
  tx9=gsub('CONDi', replace=Condi, x=tx8)
  writeLines(tx9, con=newname)
  Cond.list[Condi,] = c(paste0('Cond', Condi),r,h,i/.5,j,k,n,p)
  }}}}}}}

win = 'C:/Users/officeuser/Desktop/MMA0623/'
#mac = '/Users/luzhiming/Desktop/MMA240623/'
setwd(win)
write.csv(Cond.list,'CondList.csv', row.names = F)






















