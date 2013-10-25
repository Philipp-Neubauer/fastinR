load(file='ces.test.Rdata')

ces <- matrix(unlist(ces.test),ncol=4,byrow=T)

colnames(ces) <- c('Eveness','Condition','Sep','Error')

ces <- as.data.frame(ces)

mod <- lm(log(Error) ~ 0+Eveness+Condition*Sep,data=ces)
plot(mod)
summary(mod)

require(ggplot2)

qplot(Eveness,log(Error),color=log(Sep),data=ces,geom=c('point','smooth'),method='lm')

require(effects)

eff <- effect('Eveness',mod)
plot(eff)
