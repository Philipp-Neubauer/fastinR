context('IO test')

test_that('SI import works correctly',{
  
  SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="FASTIN")
  SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="FASTIN")
  Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="FASTIN")
  Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_sd.csv", package="FASTIN")
  
  dats <- addSI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
  expect_is(dats,'list')
  expect_is(dats$datas.SI,'list')
  expect_false(any(is.na(dats$datas.SI$preys.SI)))

})

test_that('FA import works correctly',{
  
  FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
  FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
  Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
  Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_sd.csv", package="FASTIN")
  fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package="FASTIN")
  
  dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
  expect_is(dats,'list')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))
  
})