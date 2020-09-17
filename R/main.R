
#------------------------------------------------
#' @title Check that PlasmoMAPI package has loaded successfully
#' 
#' @description Simple function to check that PlasmoMAPI package has loaded
#'   successfully. Prints "PlasmoMAPI loaded successfully!" if so.
#' 
#' @export

check_PlasmoMAPI_loaded <- function() {
  message("PlasmoMAPI loaded successfully!")
}

#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file.
#'
#' @export

pm_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'PlasmoMAPI', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

#' #------------------------------------------------
#' @title Define empty PlasmoMAPI project
#'
#' @description Define empty PlasmoMAPI project.
#'
#' @export

pm_project <- function() {
  
  # initialise project with default values
  ret <- list(data = list(coords = NULL,
                          spatial_dist = NULL,
                          stat_dist = NULL),
              map = NULL,
              output = NULL)
  
  # create custom class and return
  class(ret) <- "pm_project"
  return(ret)
}

#------------------------------------------------
#' @title Load coordinates of sampling locations into PlasmoMAPI project
#'
#' @description Given an existing PlasmoMAPI project, load the spatial coordinates of
#'   points at which samples were obtained.
#'
#' @param proj object of class \code{pm_project}.
#' @param long,lat vectors of sampling longitudes and latitudes.
#' @param check_delete_output if \code{TRUE} (the default) then check before
#'   overwriting any existing data loaded into a project.
#'
#' @export

load_coords <- function(proj, long, lat, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  assert_vector_numeric(long)
  assert_vector_numeric(lat)
  assert_same_length(long, lat)
  assert_single_logical(check_delete_output)
  
  # check whether there is data loaded into project already
  if (!is.null(proj$data$coords)) {
    
    # continuing will overwrite existing values. Return existing project if user
    # not happy to continue
    if (check_delete_output) {
      if (!user_yes_no("All existing output for this project will be lost. Continue? (Y/N): ")) {
        message("returning original project")
        return(proj)
      }
    }
    
    # inform that deleting old output
    message("overwriting data")
  }
  
  # calculate pairwise spatial distance between coords
  spatial_dist <- get_spatial_distance(long, lat)
  
  # update project with new data
  proj$data$coords <- data.frame(long = long, lat = lat)
  proj$data$spatial_dist <- spatial_dist
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# check that project has coordinates loaded
#' @noRd
pm_proj.check_coords_loaded <- function(proj) {
  assert_custom_class(proj, "pm_project")
  assert_non_null(proj$data$coords)
}

#------------------------------------------------
#' @title Create map composed of hex tiles
#'
#' @description Given an PlasmoMAPI project with sampling coordinates already loaded,
#'   create a hex grid around these points.
#'
#' @param proj object of class \code{pm_project}.
#' @param hex_width width of hexagons.
#' @param border_coords dataframe giving coordinates (long, lat) of a polygon
#'   within which the map is defined. If null then this is generated
#'   automatically from the convex hull of the sampling locations.
#'
#' @import sf
#' @importFrom grDevices chull
#' @export

create_map <- function(proj, hex_width = NULL, border_coords = NULL) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  if (!is.null(hex_width)) {
    assert_single_pos(hex_width, zero_allowed = FALSE)
  }
  if (!is.null(border_coords)) {
    assert_dataframe(border_coords)
    assert_in(c("long", "lat"), names(border_coords))
    assert_vector_numeric(border_coords$long)
    assert_vector_numeric(border_coords$lat)
  }
  
  message("Creating hex map")
  
  # calculate default hex size from data
  if (is.null(hex_width)) {
    min_range <- min(apply(proj$data$coords, 2, function(x) diff(range(x))))
    hex_width <- min_range/20
    message(sprintf("hex size chosen automatically: %s", signif(hex_width, 3)))
  }
  
  # get border_coords from convex hull of data
  if (is.null(border_coords)) {
    ch_data <- chull(proj$data$coords)
    border_coords <- proj$data$coords[c(ch_data, ch_data[1]),]
  }
  
  # get convex hull into sf polygon format
  bounding_poly <- sf::st_sfc(st_polygon(list(as.matrix(border_coords))))
  
  # make sf hex grid from poly
  hex_polys <- sf::st_make_grid(bounding_poly, cellsize = hex_width, square = FALSE)
  nhex <- length(hex_polys)
  
  # get hex centroid points
  hex_pts <- sf::st_centroid(hex_polys)
  hex_pts_df <- as.data.frame(t(mapply(as.vector, hex_pts)))
  names(hex_pts_df) <- c("long", "lat")
  
  message(sprintf("%s hexagons created", nhex))
  
  # add to project
  proj$map$hex <- hex_polys
  proj$map$hex_centroid <- hex_pts_df
  proj$map$hex_width <- hex_width
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# check that project has map loaded
#' @noRd
pm_proj.check_map_loaded <- function(proj) {
  assert_custom_class(proj, "pm_project")
  assert_non_null(proj$map)
}

#------------------------------------------------
#' @title Assign edges to the hex grid
#'
#' @description Given an PlasmoMAPI project with a hex map already created, determine
#'   which edges intersect each hex. Assumes an elliptical projection with the
#'   start- and end-points of an edge becoming the two foci of an ellipse.
#'
#' @param proj object of class \code{pm_project}.
#' @param eccentricity eccentricity of ellipses, defined as half the distance
#'   between foci divided by the semi-major axis. \eqn{e = sqrt{1 - b^2/a^2}},
#'   where \eqn{e} is the eccentricity, \eqn{a} is the length of the semi-major
#'   axis, and \eqn{b} is the length of the semi-minor axis. Eccentricity ranges
#'   between 0 (perfect circle) and 1 (straight line between foci).
#' @param assign_type Assignment method to use (1=standard, 2=new method with 
#'   function for dealing with duplicate ellipses)
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console during the permutation testing procedure.
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#'
#' @importFrom stats as.dist
#' @importFrom utils txtProgressBar
#' @export

assign_map <- function(proj, eccentricity = 0.9, assign_type = 1, 
                       report_progress = TRUE, pb_markdown = FALSE) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  pm_proj.check_map_loaded(proj)
  assert_in(assign_type,c(1,2))
  assert_single_bounded(eccentricity, inclusive_left = FALSE, inclusive_right = TRUE)
  assert_single_logical(report_progress)
  assert_single_logical(pb_markdown)
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  args_progress <- list()
  if (report_progress) {
    pb <- txtProgressBar(0, length(proj$map$hex), initial = NA, style = 3)
    args_progress <- list(pb = pb)
  }
  
  # create argument list
  args <- list(node_long = proj$data$coords$long,
               node_lat = proj$data$coords$lat,
               hex_long = proj$map$hex_centroid$long,
               hex_lat = proj$map$hex_centroid$lat,
               hex_width = proj$map$hex_width,
               eccentricity = eccentricity,
               report_progress = report_progress,
               pb_markdown = pb_markdown)
  
  
  # ---------------------------------------------
  # Run efficient C++ code, process results to project
  
  if(assign_type==1){
    output_raw <- assign_map_cpp(args, args_functions, args_progress)
    proj$map$hex_edges <- output_raw$hex_edges
  } else {
    output_raw <- assign_map2_cpp(args, args_functions, args_progress)
    proj$map$hex_edges <- output_raw$hex_edges
    proj$map$duplicate_labels <- data.frame(x=output_raw$loc_long,y=output_raw$loc_lat,label=output_raw$nDuplicates)
    # nhex=length(proj$map$hex_edges)
    # for(i in 1:nhex){
    #   proj$map$hex_edges[i][[1]]<-sort(proj$map$hex_edges[i][[1]])
    # }
  }
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# check that project has map loaded and assigned
#' @noRd
pm_proj.check_map_assigned <- function(proj) {
  assert_custom_class(proj, "pm_project")
  pm_proj.check_map_loaded(proj)
  assert_non_null(proj$map$hex_edges)
}

#------------------------------------------------
#' @title Load pairwise data into PlasmoMAPI project
#'
#' @description Load pairwise statistical data into an existing PlasmoMAPI project.
#'
#' @param proj object of class \code{pm_project}.
#' @param pairwise_data matrix of pairwise statistics between nodes.
#' @param check_delete_output if \code{TRUE} (the default) then check before
#'   overwriting any existing data loaded into a project.
#'
#' @importFrom stats as.dist
#' @export

load_data <- function(proj, pairwise_data, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  assert_matrix_numeric(pairwise_data)
  assert_symmetric_matrix(pairwise_data)
  assert_single_logical(check_delete_output)
  
  # check whether there is data loaded into project already
  if (!is.null(proj$data$stat_dist)) {
    
    # continuing will overwrite existing values. Return existing project if user
    # not happy to continue
    if (check_delete_output) {
      if (!user_yes_no("All existing output for this project will be lost. Continue? (Y/N): ")) {
        message("returning original project")
        return(proj)
      }
    }
    
    # inform that deleting old output
    message("overwriting data")
  }
  
  # update project with new data
  proj$data$stat_dist <- stats::as.dist(pairwise_data, upper = TRUE)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# check that project has pairwise data loaded
#' @noRd
pm_proj.check_data_loaded <- function(proj) {
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  assert_non_null(proj$data$stat_dist)
}

#------------------------------------------------
#' @title Calculate map hex values without permutation
#' 
#' @description TODO
#' 
#' @param proj object of class \code{pm_project}.
#' @param min_dist,max_dist minimum and maximum edge lengths to be included in
#'   the analysis. Anything outside this range is ignored.
#' 
#' @export

calc_simple_hex_values <- function(proj, min_dist=0.0,max_dist=Inf){
  
  assert_custom_class(proj, "pm_project")
  
  assert_single_pos(min_dist, zero_allowed = TRUE)
  assert_single_pos(max_dist, zero_allowed = FALSE)
  assert_gr(max_dist, min_dist)
  pm_proj.check_coords_loaded(proj)
  pm_proj.check_map_assigned(proj)
  pm_proj.check_data_loaded(proj)
  
  # Create subset of edges for which distances lie in required range
  x <- proj$data$spatial_dist
  y <- proj$data$stat_dist
  w <- which(x > min_dist & x < max_dist & !is.na(y))
  
  nHexes=length(proj$map$hex_edges)
  hex_values=rep(0,nHexes)
  for(i in 1:nHexes){
    edge_list=proj$map$hex_edges[i][[1]]
    if(is.na(edge_list[1])==FALSE){
      edge_list2=edge_list[which(edge_list %in% w)]
      hex_values[i]=mean(proj$data$stat_dist[edge_list2])
    } else {
      hex_values[i]=NA
    }
  }
  proj$output$hex_values2=hex_values
  
  invisible(proj)
  
}

#------------------------------------------------
#' @title Perform PlasmoMAPI analysis
#'
#' @description TODO.
#'
#' @param proj object of class \code{pm_project}.
#' @param n_perms number of permutations in test.
#' @param n_breaks values are broken into this many groups based on spatial
#'   distances. Permutation testing then occurs within groups. Hence, a larger
#'   value of \code{n_breaks} does a better job of conditioning on spatial
#'   distance, but comes at the cost of statistical power in permutation
#'   testing. A warning is printed if any group ends up containing fewer than 10
#'   edges.
#' @param min_dist,max_dist minimum and maximum edge lengths to be included in
#'   the analysis. Anything outside this range is ignored.
#' @param min_group_size minimum number of edges within a spatial permutation
#'   group, otherwise edges within this group are replaced with \code{NA}.
#' @param check_internal if \code{TRUE} then spatial distance is used as
#'   statistical distance in place of the loaded matrix. This option acts as an
#'   internal check, because statistical significance should not be seen when
#'   spatial distance is tested against itself. If you do see areas of
#'   significance then it is likely you are not using enough \code{n_breaks} to
#'   account for the trend with distance in the data.
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console during the permutation testing procedure.
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#'
#' @importFrom stats sd
#' @export

pm_analysis <- function(proj,
                        n_perms = 1e3,
                        n_breaks = 50, min_dist = 0, max_dist = Inf,
                        min_group_size = 5,
                        check_internal = FALSE,
                        report_progress = TRUE,
                        pb_markdown = FALSE) {
  
  # check project
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  pm_proj.check_map_assigned(proj)
  pm_proj.check_data_loaded(proj)
  
  # check other inputs
  assert_single_pos_int(n_perms, zero_allowed = FALSE)
  assert_greq(n_perms, 10, message = "minimum 10 permutations")
  assert_single_pos_int(n_breaks, zero_allowed = FALSE)
  assert_single_pos(min_dist, zero_allowed = TRUE)
  assert_single_pos(max_dist, zero_allowed = FALSE)
  assert_gr(max_dist, min_dist)
  assert_single_pos_int(min_group_size)
  assert_greq(min_group_size, 2)
  assert_single_logical(check_internal)
  assert_single_logical(report_progress)
  assert_single_logical(pb_markdown)
  
  
  # ---------------------------------------------
  # Pre-process data
  
  message("Pre-processing");
  
  # convert spatial and statistical values to x and y for convenience
  x <- as.vector(proj$data$spatial_dist)
  if (check_internal) {
    y <- as.vector(proj$data$spatial_dist)
  } else {
    y <- as.vector(proj$data$stat_dist)
  }
  
  # keep track of original index of these values to be used later in subsetting
  index <- 1:length(x)
  
  # subset to values that are within specified distance range, and are not NA
  w <- which(x > min_dist & x < max_dist & !is.na(y))
  x <- x[w]
  y <- y[w]
  index <- index[w]
  
  # break y into groups based on spatial distance
  cut_breaks <- seq(min(x), max(x), l = n_breaks + 1)
  x_cut <- cut(x, breaks = cut_breaks, include.lowest = TRUE)
  perm_group <- as.numeric(x_cut)
  y_perm <- mapply(function(i) {
    ret <- y[perm_group == i]
    if (length(ret) >= min_group_size) {
      return(ret)
    }
  }, 1:n_breaks, SIMPLIFY = FALSE)
  
  # store number of values in each bin
  df_group_num <- data.frame(dist_min = cut_breaks[-(n_breaks+1)],
                             dist_max = cut_breaks[-1],
                             n_edges = mapply(length, y_perm))
  
  # exit if all bins empty
  if (all(df_group_num$n_edges == 0)) {
    stop("not enough values within each spatial bin. Decrease value of n_breaks or min_group_size and try again")
  }
  
  # for each group, calculate the mean and sd
  y_perm_mean <- mapply(function(x) ifelse(is.null(x), NA, mean(x)), y_perm)
  y_perm_sd <- mapply(function(x) ifelse(is.null(x), NA, sd(x)), y_perm)
  
  # drop y values in groups that have zero variance
  bad_groups <- which(is.na(y_perm_sd) | (y_perm_sd == 0))
  w <- which(perm_group %in% bad_groups)
  if (any(w)) {
    y <- y[-w]
    perm_group <- perm_group[-w]
    index <- index[-w]
  }
  
  # exit if no variance in any bin
  if (all(is.na(y_perm_sd) | (y_perm_sd == 0))) {
    stop("no variance in any spatial bin")
  }
  
  # use mean and sd to normalise y values
  y_norm <- (y - y_perm_mean[perm_group])/y_perm_sd[perm_group]
  y_perm_norm <- mapply(function(i) {
    (y_perm[[i]] - y_perm_mean[i])/y_perm_sd[i]
    }, 1:n_breaks, SIMPLIFY = FALSE)
  
  # indices of edges may have changed. Update hex_edges to account for this
  hex_edges <- mapply(function(z) {
    ret <- z[z %in% index]
    ret <- match(ret, index)
    return(ret)
  }, proj$map$hex_edges, SIMPLIFY = FALSE)
  
  # calculate hex coverage
  hex_coverage <- mapply(length, hex_edges)
  
  # calculate final "observed" y values. There is one value per hex
  y_obs <- mapply(function(x) mean(y_norm[x]), hex_edges)
  
  # check that no NAs in final y_norm vector
  if (any(is.na(y_norm))) {
    stop("bug: y_norm still contains NA or NaN values")
  }
  
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  args_progress <- list()
  if (report_progress) {
    pb <- txtProgressBar(0, n_perms, initial = NA, style = 3)
    args_progress <- list(pb = pb)
  }
  
  # create argument list
  args <- list(perm_group = perm_group,
               perm_list = y_perm_norm,
               hex_edges = hex_edges,
               n_perms = n_perms,
               report_progress = report_progress,
               pb_markdown = pb_markdown)
  
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- pm_analysis_cpp(args, args_functions, args_progress)
  
  
  # ---------------------------------------------
  # Process raw output
  
  # get mean and covariance matrix of null distribution
  null_mean <- output_raw$ret_sum/n_perms
  ret_sum_sq <- rcpp_to_matrix(output_raw$ret_sum_sq)
  ret_sum_sq[upper.tri(ret_sum_sq)] <- t(ret_sum_sq)[upper.tri(ret_sum_sq)]
  cov_mat <- ret_sum_sq/n_perms - outer(null_mean, null_mean)
  
  # use null distribution to convert y_obs into a z-score
  z_score <- (y_obs - null_mean)/sqrt(diag(cov_mat))
  
  # get correlation matrix
  v <- diag(cov_mat)
  cor_mat <- cov_mat/outer(sqrt(v), sqrt(v))
  
  # get effective number of independent samples from Cheverud (2001) method
  w <- which(hex_coverage > 0)
  eg <- eigen(cor_mat[w,w])$values
  n_raw <- length(w)
  n_eff <- n_raw - (n_raw - 1)*var(eg) / n_raw
  
  # ---------------------------------------------
  # Save output as list
  
  proj$output <- list(hex_values = z_score,
                      hex_coverage = hex_coverage,
                      spatial_group_num = df_group_num,
                      n_eff = n_eff)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# check that project has output
#' @noRd
pm_proj.check_output_exists <- function(proj) {
  assert_custom_class(proj, "pm_project")
  assert_non_null(proj$output)
}

#------------------------------------------------
#' @title Get significant hexes
#'
#' @description Given a completed PlasmoMAPI analysis, return the number of hexes
#'   that pass a stated significance threshold. Calculation takes account of the
#'   inherent correlation between hexes through the effective sample size
#'   calculation.
#'
#' @param proj object of class \code{pm_project}.
#' @param empirical_tail whether to calculate empirical p-values using a
#'   one-sided test (\code{empirical_tail = "left"} or \code{empirical_tail =
#'   "right"}) or a two-sided test (\code{empirical_tail = "both"}).
#' @param alpha_raw the significance threshold used to determine significantly
#'   high/low values. This raw value is Bonferroni corrected based on the
#'   effective number of independent samples, and hence applies to the whole map
#'   and not just a single hex.
#' @param min_hex_coverage minimum coverage (number of edges assigned to a hex)
#'   for it to be included in the final result.
#'
#' @importFrom stats qnorm
#' @export

get_significant_hexes <- function(proj,
                                  empirical_tail = "both",
                                  alpha_raw = 0.05,
                                  min_hex_coverage = 10) {
  
  # check project
  assert_custom_class(proj, "pm_project")
  pm_proj.check_output_exists(proj)
  
  assert_single_string(empirical_tail)
  assert_in(empirical_tail, c("left", "right", "both"))
  assert_single_bounded(alpha_raw)
  
  # get Bonferroni-corrected upper and lower thresholds on significant z-score
  #alpha_new <- 1 - (1 - alpha_raw)^(1/proj$output$n_eff)  # (exact Bonferroni correction)
  alpha_new <- alpha_raw/proj$output$n_eff  # (approximate Bonferroni correction)
  p_bounds <- switch (empirical_tail,
                      "lower" = c(alpha_new, 1.0),
                      "upper" = c(0, 1.0 - alpha_new),
                      "both" = c(alpha_new/2, 1.0 - alpha_new/2)
  )
  z_thresh <- qnorm(p_bounds)
  
  # get which hexes (if any) are significant, and satisfy minimum hex coverage
  which_lower <- which(proj$output$hex_values < z_thresh[1] & proj$output$hex_coverage >= min_hex_coverage)
  which_upper <- which(proj$output$hex_values > z_thresh[2] & proj$output$hex_coverage >= min_hex_coverage)
  
  # return list
  ret <- list(which_lower = which_lower,
              which_upper = which_upper,
              z_thresh = z_thresh)
  return(ret)
}

#------------------------------------------------
#' @title Calculate ellipse polygon coordinates from foci and eccentricity
#'
#' @description Calculate ellipse polygon coordinates from foci and eccentricity.
#'
#' @param f1 x- and y-coordinates of the first focus.
#' @param f2 x- and y-coordinates of the first focus.
#' @param ecc eccentricity of the ellipse, defined as half the distance between
#'   foci divided by the semi-major axis. We can say \eqn{e = sqrt{1 -
#'   b^2/a^2}}, where \eqn{e} is the eccentricity, \eqn{a} is the length of the
#'   semi-major axis, and \eqn{b} is the length of the semi-minor axis.
#'   Eccentricity ranges between 0 (perfect circle) and 1 (straight line between
#'   foci).
#' @param n number of points in polygon.
#'
#' @export

get_ellipse <- function(f1 = c(-3,-2), f2 = c(3,2), ecc = 0.8, n = 100) {
  
  # check inputs
  assert_vector_numeric(f1)
  assert_length(f1, 2)
  assert_vector_numeric(f2)
  assert_length(f2, 2)
  assert_single_pos(ecc)
  assert_bounded(ecc, inclusive_left = FALSE)
  assert_single_pos_int(n)
  
  # define half-distance between foci (c), semi-major axis (a) and semi-minor
  # axis(b)
  c <- 0.5*sqrt(sum((f2-f1)^2))
  a <- c/ecc
  b <- sqrt(a^2-c^2)
  
  # define slope of ellipse (alpha) and angle of points from centre (theta)
  alpha <- atan2(f2[2]-f1[2], f2[1]-f1[1])
  theta <- seq(0, 2*pi, l = n+1)
  
  # define x and y coordinates
  x <- (f1[1]+f2[1])/2 + a*cos(theta)*cos(alpha) - b*sin(theta)*sin(alpha)
  y <- (f1[2]+f2[2])/2 + a*cos(theta)*sin(alpha) + b*sin(theta)*cos(alpha)
  
  # ensure ellipse closes perfectly
  x[n+1] <- x[1]
  y[n+1] <- y[1]
  
  # return as dataframe
  return(data.frame(x = x, y = y))
}

