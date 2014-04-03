
library(snow)

cl <- c('haast','leon','tieke','taiko','tui')

#cores per workstation
cpc <- 3

cl.act <- makeSOCKcluster(rep(cl,cpc))

#########################################
### test conversion coeff sensitivity ###
#########################################

source('CC_var_test.R')

clusterExport(cl.act, "test.cc.dep")
cc.test <- clusterApply(cl.act,rep(rep(seq(0.1,0.5,0.1),5),length(cl)*cpc),test.cc.dep,n.fats=12,n.preys=3,nsamples=20,sep=3,n.preds=1)

save(cc.test,file='cc.test.Rdata')

###############################
### test variable selection ###
###############################

source('var.select.tests.R')

arrg <- as.list(as.data.frame(matrix(rep(c(0.75,0.85,0.9,0.95,0.98,0.99),5*length(cl)*cpc),ncol=5*length(cl)*cpc),byrow=F))

clusterExport(cl.act, "var.select.tests")
var.select.test <- clusterApply(cl.act,arrg,var.select.tests,n.fats=20,n.preys=3,nsamples=20,sep=4,n.preds=1)

save(var.select.test,file='var.select.test.Rdata')

######################################
### collinearity (c), eveness (e), ###
### source separation (s)test  #######
######################################

source('ces.test.R')

arrg <- as.list(as.data.frame(matrix(rep(seq(2,4,0.25),5*length(cl)*cpc),ncol=5*length(cl)*cpc),byrow=F))

clusterExport(cl.act, "ces.tests")
ces.test <- clusterApply(cl.act,rep(rep(seq(2,4,0.25),5),length(cl)*cpc),ces.tests,n.fats=15,n.preys=3,nsamples=40,n.preds=1)

save(ces.test,file='ces.test.Rdata')

stopCluster(cl.act)

# killing in case something crashed

for (i in 1:length(cl)){
        
    system(paste('ssh',cl[i],'pkill -u philipp'))

}
