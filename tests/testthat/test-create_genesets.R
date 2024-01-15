
# simulate coordiantes for genes
set.seed(15)
trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
                                  y = runif(10, min=1, max=10), feature_name="A"),
                            cbind(x = runif(5, min=10, max=24),
                                  y = runif(5, min=1, max=10), feature_name="B"),
                            cbind(x = runif(10, min=10, max=24),
                                  y = runif(10, min=10, max=24), feature_name="C")))
trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)

data=list(trans_info=trans)

geneset_res = create_genesets(data_lst=list("rep1"= data),
                                name_lst=list(dummy_A=c("A","C"),
                                              dummy_B=c("A","B","C")),
                                bin_type="square",
                                bin_param=c(2,2),
                                w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for gene set- output mathces", {
  expect_equal(colnames(geneset_res), c("dummy_A","dummy_B"))
  expect_equal(as.vector(geneset_res$dummy_A), c(2,7,10,1))
  expect_equal(as.vector(geneset_res$dummy_B), c(2,7,11,5))
})
