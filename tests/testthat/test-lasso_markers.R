
set.seed(12)
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

data = list(trans_info = trans_info)
vecs_lst = get_vectors(data_lst=list(rep1=data), cluster_info = clusters,
                       bin_type = "square",
                       bin_param = c(20,20),
                       all_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
                       w_x = w_x, w_y=w_y)
set.seed(100)
lasso_res1 = lasso_markers(gene_mt=vecs_lst$gene_mt,
                                cluster_mt = vecs_lst$cluster_mt,
                                sample_names=c("rep1"),
                                keep_positive=TRUE,
                                coef_cutoff=0.05,
                                background=NULL)
########################################################################
test_that("lasso- one sample case", {
  expect_equal(lasso_res1$lasso_top_result$top_cluster, c("A","A","B","B"))
})

########################################################################
dummy_cl = vecs_lst$cluster_mt
colnames(dummy_cl) = c("1","B")
set.seed(100)
test_that("Return error when column names contain integers", {
  expect_error(
    lasso_markers(gene_mt=vecs_lst$gene_mt,
                  cluster_mt = dummy_cl,
                  sample_names=c("rep1"),
                  keep_positive=TRUE,
                  coef_cutoff=0.05,
                  background=NULL))
})

########################################################################
set.seed(100)
test_that("Return error when some genes non-zero entry < n_fold", {
  expect_error(
    lasso_markers(gene_mt= vecs_lst$gene_mt,
                  cluster_mt =vecs_lst$cluster_mt,
                  sample_names=c("rep1"),
                  keep_positive=TRUE,
                  coef_cutoff=0.05,n_fold = 20,
                  background=NULL))
})
########################################################################
set.seed(100)
dummy_cl = vecs_lst$cluster_mt
dummy_cl[,1] = 0
test_that("Return error when some clusters have zero variance", {
  expect_error(
    lasso_markers(gene_mt=vecs_lst$gene_mt,
                  cluster_mt = dummy_cl,
                  sample_names=c("rep1"),
                  keep_positive=TRUE,
                  coef_cutoff=0.05,
                  background=NULL))
})

dummy_ge = vecs_lst$gene_mt
dummy_ge[,1] = 0
test_that("Return error when some genes have zero variance", {
  expect_error(
    lasso_markers(gene_mt=dummy_ge,
                  cluster_mt =vecs_lst$cluster_mt,
                  sample_names=c("rep1"),
                  keep_positive=TRUE,
                  coef_cutoff=0.05,
                  background=NULL))
})
########################################################################


