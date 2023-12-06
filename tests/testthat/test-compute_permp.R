set.seed(100)
# simulate coordiantes for clusters
df_clA = data.frame(x_centroid = rnorm(n=100, mean=20, sd=5),
                    y_centroid = rnorm(n=100, mean=20, sd=5), cluster="A")
df_clB = data.frame(x_centroid = rnorm(n=100, mean=100, sd=5),
                    y_centroid = rnorm(n=100, mean=100, sd=5), cluster="B")

clusters = rbind(df_clA, df_clB)
clusters$rep="rep1"
# simulate coordiantes for genes
trans_info = data.frame(rbind(cbind(x_location = rnorm(n=100, mean=20, sd=5),
                                    y_location = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A1"),
                              cbind(x_location = rnorm(n=100, mean=20, sd=5),
                                    y_location = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A2"),
                              cbind(x_location = rnorm(n=100, mean=100, sd=5),
                                    y_location = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B1"),
                              cbind(x_location = rnorm(n=100, mean=100, sd=5),
                                    y_location = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B2")))
trans_info$x_location=as.numeric(trans_info$x_location)
trans_info$y_location=as.numeric(trans_info$y_location)

w_x =  c(min(floor(min(trans_info$x_location)),
             floor(min(clusters$x_centroid))),
         max(ceiling(max(trans_info$x_location)),
             ceiling(max(clusters$x_centroid))))
w_y =  c(min(floor(min(trans_info$y_location)),
             floor(min(clusters$y_centroid))),
         max(ceiling(max(trans_info$y_location)),
             ceiling(max(clusters$y_centroid))))

rep1 = list(trans_info = trans_info)

set.seed(100)
perm_p_lst = compute_permp(data=rep1,
                       cluster_info=clusters,
                       perm.size=100,
                       bin_type="square",
                       bin_param=c(2,2),
                       all_genes=unique(trans_info$feature_name),
                       correlation_method = "pearson",
                       n.cores=2,
                       correction_method="BH",
                       w_x=w_x ,
                       w_y=w_y)
test_that("Test permutation result - output dimension matches", {
  expect_equal(length(perm_p_lst), 4)
  expect_equal(dim(perm_p_lst$perm.arrays), c(4,2,100))
  expect_equal(dim(perm_p_lst$obs.stat), c(4,2))
  expect_equal(dim(perm_p_lst$perm.pval.adj), c(4,2))
  expect_equal(dim(perm_p_lst$perm.pval), c(4,2))
  expect_equal(names(perm_p_lst),
               c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))

})

test_that("Test permutation result - observed stat matches", {
  expect_equal(as.vector(perm_p_lst$obs.stat),
               c(1,1, -1/3, -1/3,-1/3,-1/3,1,1))
})
