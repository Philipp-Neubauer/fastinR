load('var.select.test.Rdata')

vstest <- matrix(unlist(var.select.test),ncol=2,byrow=T)
vstest[,2] <- vstest[,2]/rep(vstest[vstest[,1]==1,2],each=7)
boxplot(vstest[vstest[,1]!=0.99,2] ~ vstest[vstest[,1]!=0.99,1])
abline(h=1)

