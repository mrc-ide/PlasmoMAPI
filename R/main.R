
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
  bounding_poly <- sf::st_sfc(sf::st_polygon(list(as.matrix(border_coords))))
  
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
  
  if (assign_type == 1) {
    output_raw <- assign_map_cpp(args, args_functions, args_progress)
    proj$map$hex_edges <- output_raw$hex_edges
  } else {
    output_raw <- assign_map2_cpp(args, args_functions, args_progress)
    proj$map$hex_edges <- output_raw$hex_edges
    proj$map$duplicate_labels <- data.frame(x=output_raw$loc_long,y=output_raw$loc_lat,label=output_raw$nDuplicates)
  }
  proj$map$eccentricity = eccentricity
  
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
#' @description Runs the main PlasmoMAPI analysis. Statistical distances are
#'   binned by edge lengths and normalised within each bin. Then raw scores are
#'   calculated for eac hex and compared against their sampling distribution,
#'   obtained via permutation testing, to arrive at z-scores. A range-limited
#'   version of this analysis can be specified by setting the minimum and
#'   maximum allowed distances.
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
  
  if (report_progress) {
    message("Pre-processing");
  }
  
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
  }, seq_len(n_breaks), SIMPLIFY = FALSE)
  
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
  y_norm <- (y - y_perm_mean[perm_group]) / y_perm_sd[perm_group]
  
  # TODO - remove? Maybe fixes false positives problem
  #tmp <- df_group_num$n_edges[perm_group]
  #y_norm <- y_norm / sqrt((tmp - 1)/tmp)
  
  # apply same normalisation to values in the perm list
  y_perm_norm <- mapply(function(i) {
    (y_perm[[i]] - y_perm_mean[i]) / y_perm_sd[i]
  }, seq_len(n_breaks), SIMPLIFY = FALSE)
  
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
               y_norm = y_norm,
               report_progress = report_progress,
               pb_markdown = pb_markdown)
  
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- pm_analysis_cpp(args, args_functions, args_progress)
  
  
  # ---------------------------------------------
  # Process raw output
  
  # get mean and variance of null distribution
  null_mean <- output_raw$ret_sum/n_perms
  null_var <- (output_raw$ret_sum_sq - output_raw$ret_sum^2/n_perms) / (n_perms - 1)
  
  # use null distribution to convert y_obs into a z-score
  z_score <- (y_obs - null_mean)/sqrt(null_var)
  
  # TODO - delete. Calculating empirical p-values
  if (TRUE) {
    
    #browser()
    
    z <- do.call(rbind, output_raw$ret_all)
    
    empirical_p <- colSums(sweep(z, 2, y_obs, "<"))
    empirical_p <- (empirical_p + 1) / (nrow(z) + 2)
    z_score2 <- qnorm(empirical_p)
    
  }
  
  # ---------------------------------------------
  # Save output as list
  
  proj$output <- list(hex_values = z_score,
                      z_score = z_score,
                      z_score2 = z_score2,
                      hex_coverage = hex_coverage,
                      spatial_group_num = df_group_num)
  
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
#' @description Given a completed PlasmoMAPI analysis, return a list of which
#'   hexes are significant outliers.
#'
#' @param proj object of class \code{pm_project}.
#' @param empirical_tail whether to calculate empirical p-values using a
#'   one-sided test (\code{empirical_tail = "left"} or \code{empirical_tail =
#'   "right"}) or a two-sided test (\code{empirical_tail = "both"}).
#' @param FDR the false discovery rate, i.e. the probability that a hex
#'   identified as significant is actually a false positive.
#' @param min_hex_coverage minimum coverage (number of edges assigned to a hex)
#'   for it to be included in the final result.
#'
#' @importFrom stats pnorm
#' @export

get_significant_hexes <- function(proj,
                                  empirical_tail = "both",
                                  FDR = 0.05,
                                  min_hex_coverage = 10) {
  
  # avoid no visible binding note
  hex_coverage <- NULL
  
  # check project
  assert_custom_class(proj, "pm_project")
  pm_proj.check_output_exists(proj)
  
  assert_single_string(empirical_tail)
  assert_in(empirical_tail, c("left", "right", "both"))
  assert_single_bounded(FDR)
  assert_single_pos_int(min_hex_coverage, zero_allowed = TRUE)
  
  # get results into dataframe
  df <- data.frame(hex = seq_along(proj$output$hex_values),
                   hex_values = proj$output$hex_values,
                   hex_coverage = proj$output$hex_coverage)
  
  # subset based on coverage
  if (!any(df$hex_coverage >= min_hex_coverage)) {
    stop("no hexes wth sufficient coverage when calculating significance")
  }
  df <- subset(df, hex_coverage >= min_hex_coverage)
  
  # calculate p-values
  if (empirical_tail == "left") {
    df$p <- pnorm(df$hex_values) 
  } else if (empirical_tail == "right") {
    df$p <- pnorm(df$hex_values, lower.tail = FALSE) 
  } else if (empirical_tail == "both") {
    df$p <- 2*pnorm(-abs(df$hex_values))
  }
  
  # sort in order of increasing p
  df <- df[order(df$p),]
  
  # Bejamini and Yekutieli (2001) method for identifying significant results
  # while fixing the false descovery rate
  df$BY <- FDR * seq_along(df$p) / nrow(df)
  which_lower <- which_upper <- integer()
  if (any(df$p <= df$BY, na.rm = TRUE)) {
    
    w <- which(df$p <= df$BY)
    which_upper <- df$hex[w][df$hex_values[w] > 0]
    which_lower <- df$hex[w][df$hex_values[w] <= 0]
    
  }
  
  # return list
  ret <- list(which_lower = which_lower,
              which_upper = which_upper)
  return(ret)
}


