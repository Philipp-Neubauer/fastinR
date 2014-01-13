context('IO test')

test_that('SI import works correctly with files',{
  
  SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="FASTIN")
  SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="FASTIN")
  Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="FASTIN")
  Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_sd.csv", package="FASTIN")
  
  dats <- addSI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
  expect_is(dats,'list')
  expect_is(dats$datas.SI,'list')
  expect_false(any(is.na(dats$datas.SI$preys.SI)))
  expect_false(any(is.na(dats$datas.SI$preds.SI)))
  expect_false(any(is.na(dats$datas.SI$preym.SI)))
})

test_that('SI import works correctly FC supplied directly',{
  
  SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="FASTIN")
  SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="FASTIN")
  Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="FASTIN")
  Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_sd.csv", package="FASTIN")
  
  dats <- addSI(SI.predators=SI.predators,SI.preys=SI.preys,FC.mean=c(1,1),FC.var=c(2,2))
  expect_is(dats,'list')
  expect_is(dats$datas.SI,'list')
  expect_false(any(is.na(dats$datas.SI$preys.SI)))
  expect_false(any(is.na(dats$datas.SI$preds.SI)))
  expect_false(any(is.na(dats$datas.SI$preym.SI)))
  expect_true(all(dats$datas.SI$mean_cs==1))
  expect_true(all(dats$datas.SI$tau_cs==0.5))
  expect_true(all(dim(dats$datas.SI$mean_cs)[1]==2))
})

test_that('FA import works correctly with all files',{
  
  FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
  FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
  Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
  Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_sd.csv", package="FASTIN")
  fat.conts <- system.file("extdata", "Simdata_fat.cont.csv", package="FASTIN")
  
  dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
  expect_is(dats,'list')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))
  expect_false(any(is.na(dats$datas.FA$preym)))  
  expect_false(any(is.na(dats$datas.FA$preds)))  
})

test_that('FA import works correctly CC not from file',{
  
  FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
  FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
  Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
  Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_sd.csv", package="FASTIN")
  fat.conts <- system.file("extdata", "Simdata_fat.cont.csv", package="FASTIN")
  
  dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,CC.mean=1,CC.var=2)
  expect_is(dats,'list')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))
  expect_false(any(is.na(dats$datas.FA$preym)))  
  expect_false(any(is.na(dats$datas.FA$preds))) 
  expect_true(all(dats$datas.FA$mean_c==1))
  expect_true(all(dats$datas.FA$tau_c==0.5))
  expect_true(dim(dats$datas.FA$mean_c)[2]==10 & dim(dats$datas.FA$mean_c)[1]==3)
})
  
test_that('FA import works correctly CC and FC not from file',{
  
  FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
  FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
  Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
  Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_sd.csv", package="FASTIN")
  fat.conts <- system.file("extdata", "Simdata_fat.cont.csv", package="FASTIN")
  
  dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,CC.mean=1,CC.var=2,FC.mean=2,FC.var=2)
  expect_is(dats,'list')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))
  expect_false(any(is.na(dats$datas.FA$preym)))  
  expect_false(any(is.na(dats$datas.FA$preds))) 
  expect_true(all(dats$datas.FA$mean_c==1))
  expect_true(all(dats$datas.FA$tau_c==0.5))
  expect_true(dim(dats$datas.FA$mean_c)[2]==10 & dim(dats$datas.FA$mean_c)[1]==3)
  expect_true(length(dats$datas.FA$fc_mean)==3)
  
})