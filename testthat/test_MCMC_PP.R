context('MCMC pop.prop testing')

test_that('Pop Props with FA give sensible answers',{
  
  data('Sim',envir = environment())
  
  dats <- selectvars(datas,ix=c(2,7,3,6,8))
  
  dat <- run_MCMC(dats,nIter=10000,nBurnin=1000,nChains=3,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.1,plott=F)
  
  expect_is(dat,"pop_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  expect_false(sum(abs(colMeans(props)-colMeans(do.call('rbind',dat))))>0.2)
  
})

test_that('Pop Props with SI give sensible answers',{
  
  data('Sim',envir = environment())
  
  dat <- run_MCMC(datas,nIter=10000,nBurnin=1000,nChains=3,nThin=10,Data.Type='Stable.Isotopes',Analysis.Type='Population.proportions',even=0.1,plott=F)
  
  expect_is(dat,"pop_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  expect_true(sum(abs(colMeans(props)-colMeans(do.call('rbind',dat))))>0.2)
  
})

test_that('Pop Props with combined give sensible answers',{
  
  data('Sim',envir = environment())
  
  dats <- selectvars(datas,ix=c(2,7,3,6,8))
  
  dat <- run_MCMC(dats,nIter=10000,nBurnin=1000,nChains=3,nThin=10,Data.Type='Combined.Analysis',Analysis.Type='Population.proportions',even=0.1,plott=F)
  
  expect_is(dat,"pop_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="FASTIN")
  props <- read.csv(prop,header=F,row.names=1)
  expect_false(sum(abs(colMeans(props)-colMeans(do.call('rbind',dat))))>0.2)
  
})