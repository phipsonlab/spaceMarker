
# simulate coordiantes for clusters
clusters = data.frame(x = c(1,2,20,21,22,23,24),
                    y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"

w_x=c(0,25)
w_y=c(0,25)
vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       bin_type = "square",
                       bin_param = c(2,2),
                       all_genes = NULL,
                       w_x = w_x, w_y=w_y)
test_that("Test can only vectorise clusters - output length matches", {
  expect_equal(length(vecs_lst), 1)
})


test_that("Test can only vectorise clusters - output vector matches", {
  expect_equal(as.vector(vecs_lst$cluster_mt), c(2,0,0,5))
})

#############################################################################
# simulate coordiantes for genes
trans = data.frame(x = c(1,2,20,21,22,23,24),
                   y = c(23, 24, 1,2,3,4,5),
                   feature_name="A")
data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"=data),
                            cluster_info = NULL,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A"),
                            w_x = w_x, w_y=w_y)
test_that("Test can only vectorise genes - output length mathces", {
  expect_equal(length(vecs_lst_gene), 1)
})

#############################################################################
# simulate coordiantes for genes
trans = as.data.frame(rbind(cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="A"),
                            cbind(x = c(1,20),
                                  y = c(15, 10),
                                  feature_name="B"),
                            cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="C")))

trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)
clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="rep1"

data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"= data),
                            cluster_info = clusters,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A","B","C"),
                            w_x = w_x, w_y=w_y)


test_that("Test can vectorise genes and clusters - output mathces", {
  expect_equal(as.vector(vecs_lst_gene$cluster_mt), c(2,0,2,3))
  expect_equal(as.vector(vecs_lst_gene$gene_mt), c(2, 0, 0 ,5 ,
                                                 1,0, 0, 1,
                                                 2, 0, 0, 5))
})
#############################################################################
# generate gene vector from count matrix
cm <- data.frame(rbind("gene_A"=c(0,0,2,0,0,0,2),
                 "gene_B"=c(5,3,3,13,0,1,14),
                 "gene_C"=c(5,0,1,5,1,0,7),
                 "gene_D"=c(0,1,1,2,0,0,2)))
colnames(cm)= paste("cell_", 1:7, sep="")

# simulate coordiantes for clusters
clusters = data.frame(x = c(1, 2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"
clusters$cell_id= colnames(cm)
# simulate coordiantes for genes
w_x=c(0,25)
w_y=c(0,25)
# cell_1 = (1,0,0,0)
# cell_2 = (1,0,0,0)
# cell_3 = (0,0,0,1)
# cell_4 = (0,0,0,1)
# cell_5 = (0,0,0,1)
# cell_6 = (0,0,0,1)
# cell_7 = (0,0,0,1)

vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       cm_lst=list(rep1=cm),
                       bin_type = "square",
                       bin_param = c(2,2),
                       all_genes = row.names(cm),
                       w_x = w_x, w_y=w_y)

test_that("Test can use count matrix (square) - output vector matches", {
  expect_equal(as.vector(vecs_lst$gene_mt),
               c(0,0,0,4,
                 8,0,0,31,
                 5,0,0,14,
                 1,0,0,5))

})
# cell_1 = c()
vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       cm_lst=list(rep1=cm),
                       bin_type = "hexagon",
                       bin_param = c(10),
                       all_genes = row.names(cm),
                       w_x = w_x, w_y=w_y)

test_that("Test can use count matrix (hex) - output vector matches", {
  expect_equal(as.vector(vecs_lst$gene_mt),
               c(0,0,0,0,4,0,0,
                 0,0,0,0,31,8,0,
                 0,0,0,0,14,5,0,
                 0,0,0,0,5,1,0))

})
