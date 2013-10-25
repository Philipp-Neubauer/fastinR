library(snow)

cl <- c('haast','leon','robin','tieke','taiko','tui')

#cores per workstation
cpc <- 3

cl.act <- makeSOCKcluster(rep(cl,cpc))

#########################################
### test conversion coeff sensitivity ###
#########################################

clusterExport(cl.act, "test.cc.dep")
cc.test <- clusterApply(cl.act,rep(seq(0,0.5,0.1),length(cl)*cpc),test.cc.dep,n.fats=12,n.preys=3,nsamples=20,sep=3,n.preds=5)

save(cc.test,file='cc.test.Rdata')

###############################
### test variable selection ###
###############################

clusterExport(cl.act, "var.select.tests")
var.select.test <- clusterApply(cl.act,rep(c(0.75,0.9,0.95,0.98,0.99),length(cl)*cpc),var.select.tests,n.fats=15,n.preys=3,nsamples=20,sep=4,n.preds=6)

save(var.select.test,file='var.select.test.Rdata')
######################################
### collinearity (c), eveness (e), ###
### source separation (s)test  #######
######################################

clusterExport(cl.act, "ces.tests")
ces.test <- clusterApply(cl.act,seq(0.5,5,0.25),length(cl)*cpc),ces.tests,n.fats=12,n.preys=5,nsamples=20,sep=3,n.preds=10)

stopCluster(cl.act)
