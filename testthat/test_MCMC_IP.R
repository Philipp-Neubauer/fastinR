context('MCMC Ind.Prop testing')

test_that('Ind Props with FA give sensible answers',{
  
  data('Sim',envir = environment())
  
  dats <- select_vars(datas,ix=c(2,7,3,6,8))
  
  dat <- run_MCMC(dats,nIter=10000,nBurnin=5000,nChains=3,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Individual.proportions',even=0.2,plott=T)
  
  
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat))[1:3],matrix(colMeans(do.call('rbind',dat))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_false(max(abs(props-mmat))>0.2)
  
})

test_that('Pop Props with SI give sensible answers',{
  
  data('Sim',envir = environment())
  
  dat <- run_MCMC(datas,nIter=10000,nBurnin=1000,nChains=3,nThin=10,Data.Type='Stable.Isotopes',Analysis.Type='Individual.proportions',even=0.2,plott=F)
  
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat))[1:3],matrix(colMeans(do.call('rbind',dat))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_true(max(abs(props-mmat))>0.2)
  
})

test_that('Pop Props with combined give sensible answers',{
  
  data('Sim',envir = environment())
  
  dats <- select_vars(datas,ix=c(2,7,3,6,8))
  
  dat <- run_MCMC(dats,nIter=10000,nBurnin=1000,nChains=3,nThin=10,Data.Type='Combined.Analysis',Analysis.Type='Individual.proportions',even=0.2,plott=F)
  
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat))[1:3],matrix(colMeans(do.call('rbind',dat))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_false(max(abs(props-mmat))>0.2)
  
})