
#------------------------------------------------
test_that("get_barrier_intersect() works for complex barrier shapes", {
  
  # define node distribution
  node_df <- expand.grid(long = seq(-5, 5, 5),
                         lat = seq(-5, 5, 5))
  
  # define barrier shape
  barrier_list <- list(data.frame(long = c(-4,-4,4,4,3,3,-3,-3,-4), lat = c(-4,4,4,-4,-4,3,3,-4,-4)))
  
  # expect no warnings when drawing pairwise distances
  expect_silent(get_barrier_intersect(node_df$long, node_df$lat,
                                      barrier_list = barrier_list,
                                      barrier_penalty = 1,
                                      barrier_method = 1))
  
})
