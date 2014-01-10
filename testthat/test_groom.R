context('Variable selection test')

test_that('Variable selection works',{
  
  data('SimEx')
  
  dats <- selectvars(datas,ix=1:3)
  expect_is(dats,'list')
  expect_is(dats$datas.FA,'list')
  expect_false(any(is.na(dats$datas.FA$preys)))  
  expect_false(any(is.na(dats$datas.FA$preym)))  
  expect_false(any(is.na(dats$datas.FA$preds)))  
  expect_true(ncol(dats$datas.FA$preys)==3)  
  expect_true(ncol(dats$datas.FA$preym)==2)  
  expect_true(ncol(dats$datas.FA$preds)==2)  
  expect_true(ncol(dats$datas.FA$preds.FA)==3)  
  
})