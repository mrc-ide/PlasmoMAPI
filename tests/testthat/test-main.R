
#------------------------------------------------
test_that("exits gracefuly when not enough values (or SD=0) in spatial distance groups", {
  
  set.seed(1)
  
  # define node distribution
  L <- 5
  n_deme <- 1e2
  node_df <- data.frame(long = runif(n_deme, -L, L),
                        lat = runif(n_deme, -L, L))
  
  # define statistical distance equal to euclidian distance for convenience
  stat_distance <- as.matrix(dist(node_df))
  
  # create new PlasmoMAPI project
  p <- pm_project()
  
  # load node coordinates
  p <- load_coords(p, node_df$long, node_df$lat)
  
  # set up map
  p <- create_map(p)
  
  # assign edges to hexes
  p <- assign_map(p, eccentricity = 0.9)
  
  # load data
  p <- load_data(p, stat_distance, check_delete_output = FALSE)
  
  # run analysis with min_group_size=10
  p <- pm_analysis(p, n_perms = 1e1, n_breaks = 50,
                   min_dist = 0, max_dist = Inf,
                   min_group_size = 10)
  
  # should not see any distance group with more than 10 edges in output
  n_edges <- p$output$spatial_group_num$n_edges
  n_edges[n_edges == 0] <- NA
  expect_true(min(n_edges, na.rm = TRUE) >= 10)
  
  # should fail when running analysis with min_group_size=1000
  expect_error(pm_analysis(p, n_perms = 1e1, n_breaks = 50,
                           min_dist = 0, max_dist = Inf,
                           min_group_size = 1000))
  
  # load new stat data that has zero variance
  stat_distance <- matrix(1, n_deme, n_deme)
  p <- load_data(p, stat_distance, check_delete_output = FALSE)
  
  # should fail due to lack of variance in and spatial group
  expect_error(pm_analysis(p, n_perms = 1e1, n_breaks = 50,
                           min_dist = 0, max_dist = Inf,
                           min_group_size = 2))
  
})
