test_that("`plot_vectors()` plots as expected", {
    # simulate coordiantes for clusters
    clusters = data.frame(x = c(1,2,20,21,22,23,24),
                          y = c(23, 24, 1,2,3,4,5), cluster="A")

    w_x=c(0,25)
    w_y=c(0,25)
    vdiffr::expect_doppelganger(title = "plot_vectors",
                                fig = plot_vectors(loc_mt=clusters, to_plot=c("A"),
                                to_plot_column_name="cluster",
                                bin_type="square",bin_param=c(2,2),
                                w_x=w_x, w_y=w_y)

  )
})
