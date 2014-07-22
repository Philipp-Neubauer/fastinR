 context('Variable selection test')

test_that('Variable selection works',{
  
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
  
  dats <- select_vars(datas,ix=1:3)
  expect_is(dats,'Combined_Markers')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))  
  expect_false(any(is.na(dats$datas.FA$preym)))  
  expect_false(any(is.na(dats$datas.FA$preds)))  
  expect_true(ncol(dats$datas.FA$preys)==3)  
  expect_true(ncol(dats$datas.FA$preym)==2)  
  expect_true(ncol(dats$datas.FA$preds)==2)  
  expect_true(ncol(dats$datas.FA$preds.FA)==3)  
  
})