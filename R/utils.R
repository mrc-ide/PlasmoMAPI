
#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
#' @title Print unclassed object
#'
#' @description Print object after unclassing, thereby removing any custom print
#'   method.
#'
#' @param x object to print in full.
#'
#' @export

print_full <- function(x) {
  print(unclass(x))
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

#------------------------------------------------
#' @title Get great circle distance between spatial points
#'
#' @description Get great circle distance between spatial points, defined by a
#'   vector of longitudes and latitudes. Distances are returned in a pairwise
#'   distance matrix.
#' 
#' @param long,lat vector of longitudes and latitudes.
#'
#' @export

get_spatial_distance <- function(long, lat) {
  
  # check inputs
  assert_vector(long)
  assert_numeric(long)
  assert_vector(lat)
  assert_numeric(lat)
  assert_same_length(long, lat)
  
  # calculate distance matrix
  ret <- apply(cbind(long, lat), 1, function(y) {lonlat_to_bearing(long, lat, y[1], y[2])$gc_dist})
  ret <- as.dist(ret, upper = TRUE)
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates, defined by longitude and latitude of both origin and
#'   destination points.
#'
#' @param origin_lon,origin_lat the origin longitude and latitude.
#' @param dest_lon,dest_lat the destination longitude and latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # check inputs
  assert_vector(origin_lon)
  assert_numeric(origin_lon)
  assert_vector(origin_lat)
  assert_numeric(origin_lat)
  assert_vector(dest_lon)
  assert_numeric(dest_lon)
  assert_vector(dest_lat)
  assert_numeric(dest_lat)
  
  # convert input arguments to radians
  origin_lon <- origin_lon*2*pi/360
  origin_lat <- origin_lat*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  
  # get change in lon
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat) - sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
}

#------------------------------------------------
#' @title Get distance between points taking into account barriers
#'
#' @description Given a set of lat/lon coordinates and a list of barriers in the
#'   form of polygons, returns the "distance" between points where distance is
#'   equal to the great-circle distance with a penalty applied if the line
#'   intersects a barrier. The exact way in which barriers modify distances can
#'   be varied (see \code{barrier_method} argument).
#'
#' @param node_long longitudes of nodes.
#' @param node_lat latitudes of nodes.
#' @param barrier_list list of polygons representing barriers. Each element of
#'   the list must be a dataframe with columns \code{long} and \code{lat}
#'   specifying the coordinates of points that make up the polygon. Polygons
#'   must be complete rings, meaning the final row of the dataframe must equal
#'   the first row.
#' @param barrier_penalty penalty values of each barrier. If a single value is
#'   provided then this value will be used for all barriers.
#' @param barrier_method the method by which penalties are applied:
#'   \enumerate{
#'     \item{compare pairwise lines to barriers. If the line intersects then add
#'     a fixed \code{barrier_penalty} to the spatial distance.}
#'     \item{compare pairwise lines to barriers. Calculate the intersection of
#'     the two, multiply this by the \code{barrier_penalty} and add to the
#'     spatial distance. For example, a \code{barrier_penalty} of 1 would mean
#'     there is double the "friction" when moving through a barrier.}
#'     \item{compare pairwise ellipses to barriers. Calculate the intersection
#'     area of the two, multiply this by the \code{barrier_penalty} and add to
#'     the spatial distance.}
#'   }
#' @param max_barrier_range edges that are longer than this distance are
#'   unaffected by any barriers. Makes it possible to model barriers that only
#'   apply locally.
#' @param eccentricity eccentricity of ellipses (only used under
#'   \code{barrier_method = 3}).
#' @param noise_sd standard deviation of Gaussian noise added to all distances
#'   (after the application of barriers).
#' @param n_ell number of points that make up an ellipse (only used under
#'   \code{barrier_method = 3}).
#'
#' @import sf
#' @importFrom stats dist rnorm
#' @export

get_barrier_intersect <- function(node_long,
                                  node_lat,
                                  barrier_list = list(),
                                  barrier_penalty = numeric(),
                                  barrier_method = 1,
                                  max_barrier_range = Inf,
                                  eccentricity = 0.9,
                                  noise_sd = 0,
                                  n_ell = 20) {
  
  # check inputs
  assert_vector_numeric(node_long)
  assert_vector_numeric(node_lat)
  assert_same_length(node_long, node_lat)
  assert_list(barrier_list)
  nb <- length(barrier_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(barrier_list[[i]])
      assert_in(c("long", "lat"), names(barrier_list[[i]]))
      assert_eq(barrier_list[[i]][1,], barrier_list[[i]][nrow(barrier_list[[i]]),], 
                message = "barrier polygons must be closed, i.e. the last node coordinate equals the first")
    }
  }
  assert_vector_numeric(barrier_penalty)
  assert_single_pos_int(barrier_method)
  assert_in(barrier_method, 1:3)
  assert_single_pos(max_barrier_range, zero_allowed = TRUE)
  assert_single_bounded(eccentricity, inclusive_left = FALSE)
  assert_single_pos(noise_sd, zero_allowed = TRUE)
  assert_single_pos_int(n_ell, zero_allowed = FALSE)
  
  # force barrier_penalty to vector
  barrier_penalty <- force_vector(barrier_penalty, length(barrier_list))
  assert_same_length(barrier_penalty, barrier_list)
  
  # create mask for ignoring edges greater then
  distance_mask <- 1
  if (is.finite(max_barrier_range)) {
    d <- as.vector(get_spatial_distance(node_long, node_lat))
    distance_mask <- (d < max_barrier_range)
  }
  
  # apply barrier penalties
  intersect_penalty <- 0
  if (nb > 0 & any(barrier_penalty != 0)) {
    
    # convert barrier list to st_polygon
    poly_list <- list()
    for (i in 1:length(barrier_list)) {
      poly_list[[i]] <- sf::st_polygon(list(as.matrix(barrier_list[[i]])))
    }
    
    # get node coordinates in matrix
    node_mat <- cbind(node_long, node_lat)
    
    # if comparing lines
    if (barrier_method %in% c(1,2)) {
      
      # create all pairwise sf_linestring between nodes
      line_list <- list()
      n_node <- length(node_long)
      i2 <- 0
      for (i in 1:(n_node-1)) {
        for (j in (i+1):n_node) {
          i2 <- i2 + 1
          line_list[[i2]] <- sf::st_linestring(node_mat[c(i,j),])
        }
      }
      
      # convert lines and polys to st_sfc
      line_sfc <- sf::st_sfc(line_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(line_sfc, poly_sfc))
      
      # convert to (great circle) length of intersection if using method 2
      if (barrier_method == 2) {
        intersect_mat[intersect_mat == TRUE] <- mapply(function(x) {
          get_spatial_distance(x[c(1,3)], x[c(2,4)])[1]
        }, sf::st_intersection(line_sfc, poly_sfc))
      }
    }
    
    # mask out edges that are beyond limit distance
    intersect_mat <- sweep(intersect_mat, 1, distance_mask, "*")
    
    # if comparing ellipse
    if (barrier_method == 3) {
      
      # create all pairwise ellipses between nodes
      ell_list <- list()
      n_node <- length(node_long)
      i2 <- 0
      for (i in 1:(n_node-1)) {
        for (j in (i+1):n_node) {
          i2 <- i2 + 1
          ell_df <- get_ellipse(f1 = node_mat[i,], f2 = node_mat[j,], ecc = eccentricity, n = n_ell)
          ell_list[[i2]] <- sf::st_polygon(list(as.matrix(ell_df)))
        }
      }
      
      # convert ellipses and polys to st_sfc
      ell_sfc <- sf::st_sfc(ell_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(ell_sfc, poly_sfc))
      
      # mask out ellipses that are beyond limit distance
      intersect_mat <- sweep(intersect_mat, 1, distance_mask, "*")
      
      # convert to area of intersection
      intersect_areas <- mapply(function(x) {
        sf::st_area(x)
      }, sf::st_intersection(ell_sfc, poly_sfc))
      intersect_mat[intersect_mat == TRUE] <- intersect_areas[intersect_mat == TRUE]
      
    }
    
    # apply penalty
    intersect_penalty <- rowSums(sweep(intersect_mat, 2, barrier_penalty, '*'))
    
  }  # end apply barrier penalties
  
  # get pairwise distance plus penalty
  d <- get_spatial_distance(node_long, node_lat) + intersect_penalty
  d <- d + rnorm(length(d), sd = noise_sd)
  
  # return matrix
  return(as.matrix(d))
}

#------------------------------------------------
# pass in a series of polygons (class sfc_POLYGON). Expand by a buffer distance
# d and merge polygons
#' @noRd
get_merged_poly <- function(hex_polys, d = 0.1) {
  
  # expand polygons
  ret <- st_buffer(hex_polys, d)
  
  # merge polygons
  ret <- sf::st_union(sf::st_sf(ret))
  
  # undo expansion
  ret <- st_buffer(ret, -d)
  
  return(ret)
}
