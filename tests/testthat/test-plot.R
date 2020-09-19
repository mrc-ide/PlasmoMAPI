
#------------------------------------------------
test_that("plot_network produces error with bad input", {
  
  set.seed(1)
  
  # create new PlasmoMAPI project
  p <- pm_project()
  
  # expect error when run with no coords loaded
  expect_error(plot_network(p))
  
})
