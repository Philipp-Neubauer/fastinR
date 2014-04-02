load(file='ces.test.Rdata')

ces <- matrix(unlist(ces.test),ncol=4,byrow=T)

colnames(ces) <- c('Eveness','Condition','Sep','Error')

ces[,1] <- log(ces[,1]/(1-ces[,1]))/(2*sd(ces[,1]))
ces[,1:3] <- apply(ces[,1:3],2,function(x){(x-mean(x))/(2*sd(x))})


ces <- as.data.frame(ces)
attach(ces)

mod <- lm(log(Error) ~ Eveness*Sep*Condition,data=ces)
plot(mod)
summary(step(mod))

require(ggplot2)

qplot(Eveness,log(Error),color=log(Sep),data=ces,geom=c('point','smooth'),method='lm')

require(effects)

eff <- effect('Eveness*Sep',mod)

plot(eff,main='')
