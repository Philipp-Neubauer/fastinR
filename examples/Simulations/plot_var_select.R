load('./examples//Simulations/var.select.test.Rdata')

require(reshape2)

vstest <- do.call(rbind,var.select.test)
vstest2 <- melt(var.select.test)
vstest[,2] <- vstest[,2]/rep(vstest[vstest[,1]==0.99,2],each=6)
boxplot(vstest[vstest[,1]!=0.99,2] ~ vstest[vstest[,1]!=0.99,1])
abline(h=1)

vst <- cbind(vstest,subset(vstest2,Var2==2))
names(vst)<-c('sim','err','val','ob','err2','L1')

model <- lm((err) ~  (val), data=data.frame(vst))
grid <- with(vst, expand.grid(
  val = unique((val))))

grid$err<- stats::predict(model, newdata=grid)

errs <- stats::predict(model, newdata=grid, se = TRUE)
grid$ucl <- errs$fit + 1.96 * errs$se.fit
grid$lcl <- errs$fit - 1.96 * errs$se.fit

require(ggplot2)
m <- ggplot(data=vst,aes(y=err,x=(val),col=rep(tapply(log(err),factor(L1),mean),each=6),group=factor(L1))) + geom_point() + scale_y_log10()+scale_x_discrete(labels=unique(vst$sim)) +xlab('variation retained') + ylab('log(error)')
m <- m+ geom_line(linetype=3)+  theme(legend.position="none")
m <- m+ scale_colour_gradient(low='green',high='red')
m + stat_smooth(aes(y=err,x=(val),group=NULL),data=vst,alpha=0.2,col=1) 

ggsave('./examples//Simulations/varselectplot.pdf',scale=0.6)
