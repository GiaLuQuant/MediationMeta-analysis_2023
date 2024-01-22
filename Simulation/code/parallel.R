library('doParallel')
wd.mac = '/Users/luzhiming/000我/sysu/自己/主要/MMA/MMA_revision/sim/MMA240105/code/'

setwd(wd.mac)

myfun <- function(i){source(paste0('_RC',i,'.r'))}

cores = detectCores()-1

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (92:96)) %dopar% myfun(i)
stopCluster(clust)


clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*7) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*8) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*9) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*10) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*11) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*12) %dopar% myfun(i)
stopCluster(clust)




clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*0) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*1) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*2) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*3) %dopar% myfun(i)
stopCluster(clust)


clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*4) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*5) %dopar% myfun(i)
stopCluster(clust)

clust <- makeCluster(cores)
registerDoParallel(clust)
result = foreach(i = (1:cores)+cores*6) %dopar% myfun(i)
stopCluster(clust)

