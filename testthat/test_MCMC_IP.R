context('MCMC Ind.Prop testing')

test_that('Ind Props with FA give sensible answers',{
  
  FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="fastinR")
  FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="fastinR")
  Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="fastinR")
  Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_var.csv", package="fastinR")
  fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package="fastinR")
  SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="fastinR")
  SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="fastinR")
  Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="fastinR")
  Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package="fastinR")
  
  dats <- add_FA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
  datas <- add_SI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var,datas=dats)
    
  dats <- select_vars(datas,ix=c(2,7,3,6,8))
  
  dat <- run_MCMC(dats,nIter=1000,nBurnin=500,nChains=3,nThin=1,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Individual.proportions',even=0.2,plott=F)
    
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="fastinR")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat[1:3]))[1:3],matrix(colMeans(do.call('rbind',dat[1:3]))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_false(max(abs(props-mmat))>0.2)
  
})

test_that('Pop Props with SI give sensible answers',{
  
  dat <- run_MCMC(datas,nIter=1000,nBurnin=500,nChains=3,nThin=1,Data.Type='Stable.Isotopes',Analysis.Type='Individual.proportions',even=0.2,plott=F)
  
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="fastinR")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat[1:3]))[1:3],matrix(colMeans(do.call('rbind',dat[1:3]))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_true(max(abs(props-mmat))>0.2)
  
})

test_that('Pop Props with combined give sensible answers',{
  
  
  dat <- run_MCMC(dats,nIter=1000,nBurnin=500,nChains=3,nThin=1,Data.Type='Combined.Analysis',Analysis.Type='Individual.proportions',even=0.2,plott=F)
  
  expect_is(dat,"ind_props")
  
  prop <- system.file("extdata", "Simdata_props.csv", package="fastinR")
  props <- read.csv(prop,header=F,row.names=1)
  mmat <- rbind(colMeans(do.call('rbind',dat[1:3]))[1:3],matrix(colMeans(do.call('rbind',dat[1:3]))[4:33],10,3))
  props  <- rbind(colMeans(props),props)
  expect_false(max(abs(props-mmat))>0.2)
  
})