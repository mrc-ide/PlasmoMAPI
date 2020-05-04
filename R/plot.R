
#------------------------------------------------
#' @title Red to blue colours
#'
#' @description Simple sequence of red-to-blue colours.
#'
#' @param n the number of colours.
#'
#' @importFrom grDevices colorRampPalette
#' @export

col_hotcold <- function(n = 6) {
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_pal <- colorRampPalette(raw_cols)
  return(my_pal(n))
}

#------------------------------------------------
# a series of internally-used colours
#' @noRd
daily_cols <- function() {
  c("firebrick1", "chartreuse3", "dodgerblue", "dodgerblue4", "purple", "darkorange", "firebrick4")
}

#------------------------------------------------
#' @title Plot pairwise spatial distance against loaded statistic
#'
#' @description Plot pairwise spatial distance against loaded statistic.
#'
#' @param proj object of class \code{pm_project}.
#' @param col the colour of points.
#' @param overlay_model whether to overlay the model fit, if \code{fit_model()}
#'   has been called.
#' 
#' @import ggplot2
#' @export

plot_dist <- function(proj, col = "#00000050", overlay_model = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  assert_length(col, 1)
  assert_single_logical(overlay_model)
  
  # create basic plot
  df_plot <- data.frame(spatial <- as.vector(proj$data$spatial_dist), stat <- as.vector(proj$data$stat_dist))
  plot1 <- ggplot(df_plot) + theme_bw()
  
  # add points
  plot1 <- plot1 + geom_point(aes_(x = ~spatial, y = ~stat), data=df_plot, col = col)
  
  # titles etc
  plot1 <- plot1 + xlab("spatial distance") + ylab("statistical distance")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot hex coverage
#'
#' @description Given an PlasmoMAPI project with a map already assigned, plot
#'   the coverage (number of intersecting ellipses) of each hex. This plot acts
#'   as a diagnostic, as hexes with low coverage will conform less well to the
#'   assumptions of the permutation test (see \emph{details}).
#'
#' @details Good coverage is needed to ensure the validity of the statistical
#'   procedure. For each hex, the observed value is the mean of the normalised
#'   values of the edges that intersect it. Under the null model that these
#'   normalised values have no systematic bias, the distribution of this test
#'   statistic is approximately normal as a consequence of the central limit
#'   theorem. Hence, we can use a permutation test to characterise the mean and
#'   standard deviation of this null distribution, then we can quantify how
#'   extreme the observed value is in terms of its z-score. However, if coverage
#'   is too low then the null distribution will not be normally distributed, and
#'   hence the z-score will not be an accurate description of how extreme the
#'   observed data really is.
#'
#' @param proj object of class \code{pm_project}.
#' @param breaks the sequence of coverage breaks used.
#' 
#' @import sf
#' @import ggplot2
#' @export

plot_coverage <- function(proj, breaks = c(0,10,20,30,40,50,100,Inf)) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  pm_proj.check_coords_loaded(proj)
  pm_proj.check_map_assigned(proj)
  pm_proj.check_output_exists(proj)
  assert_vector_pos(breaks)
  assert_increasing(breaks)
  assert_greq(length(breaks), 2)
  
  # bin hex coverage
  hex_coverage <- proj$output$hex_coverage
  intersect_bin <- cut(hex_coverage, breaks = breaks)
  
  # -------- Map --------
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add hexs
  plot1 <- plot1 + geom_sf(aes_(fill = ~intersect_bin), color = NA, data = proj$map$hex)
  
  # add points
  coords <- data.frame(long <- proj$coords$long, lat <- proj$coords$lat)
  plot1 <- plot1 + geom_point(aes_(x = ~long, y = ~lat), data = coords, size = 0.5)
  
  # titles and legends
  plot1 <- plot1 + scale_fill_manual(values = col_hotcold(length(breaks)-1),
                                     limits = levels(intersect_bin),
                                     name = "hex coverage")
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  plot1 <- plot1 + guides(fill = guide_legend(reverse = TRUE))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot hex map of PlasmoMAPI output
#'
#' @description Plot hex map of PlasmoMAPI output.
#'
#' @param proj object of class \code{pm_project}.
#' @param col_scale the colour scale to use.
#' @param plot_sampling_points whether to overlay sampling locations.
#' @param plot_hex_grid whether to plot the hex map. If the project contains
#'   output then hexes will be coloured based on z-scores, otherwise a flat
#'   colour scheme will be used.
#' @param min_hex_coverage minimum coverage (number of edges assigned to a hex)
#'   for it to be included in the final result, otherwise these hexes are given
#'   the value \code{NA}.
#' @param plot_significance whether to outline areas that were identified as
#'   significant outliers.
#' @param empirical_tail if plotting significance, whether to calculate
#'   empirical p-values using a one-sided test (\code{empirical_tail = "left"}
#'   or \code{empirical_tail = "right"}) or a two-sided test
#'   (\code{empirical_tail = "both"}).
#' @param alpha_raw the significance threshold used to determine significantly
#'   high/low values. This raw value is Bonferroni corrected based on the
#'   effective number of independent samples, and hence applies to the whole map
#'   and not just a single hex.
#' @param base_plot optional base plot (object of class \code{ggplot}) on which
#'   this function builds. If \code{NULL} then a simple empty plot is used.
#' @param poly_list optional list of polygon coordinates that are added to plot.
#' @param labeled_points optional data frame of labeled points to add to graph
#' @param plot_data_values whether to plot aggregated data values instead of 
#'                         z-score, if available
#' 
#' @import ggplot2
#' @importFrom viridisLite magma
#' @export

plot_map <- function(proj,
                     col_scale = viridisLite::magma(100),
                     plot_sampling_points = TRUE,
                     plot_hex_grid = TRUE,
                     min_hex_coverage = 10,
                     plot_significance = TRUE,
                     empirical_tail = "both",
                     alpha_raw = 0.05,
                     base_plot = NULL,
                     poly_list = list(),
                     labeled_points=NULL,
                     plot_data_values = FALSE) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  assert_single_logical(plot_sampling_points)
  assert_single_logical(plot_hex_grid)
  assert_single_pos_int(min_hex_coverage, zero_allowed = TRUE)
  assert_single_logical(plot_significance)
  assert_single_string(empirical_tail)
  assert_in(empirical_tail, c("left", "right", "both"))
  assert_single_bounded(alpha_raw)
  if (!is.null(base_plot)) {
    assert_custom_class(base_plot, "ggplot")
  }
  assert_list(poly_list)
  nb <- length(poly_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(poly_list[[i]])
      assert_in("long", names(poly_list[[i]]))
      assert_in("lat", names(poly_list[[i]]))
      assert_eq(poly_list[[i]][1,], poly_list[[i]][nrow(poly_list[[i]]),], 
                message = "barrier polygons must be closed, i.e. the last node coordinate equals the first")
    }
  }
  if(is.null(labeled_points)==FALSE){
    assert_dataframe(labeled_points)
    assert_vector_numeric(labeled_points$x)
    assert_vector_numeric(labeled_points$y)
  }
  assert_single_logical(plot_data_values)
  
  # determine which aspects can/should be plotted
  if (plot_sampling_points) {
    plot_sampling_points <- !is.null(proj$data$coords)
  }
  if (plot_hex_grid) {
    plot_hex_grid <- !is.null(proj$map$hex)
  }
  if(plot_data_values){
    plot_hex_values <- !is.null(proj$output$hex_values2) 
  } else {
    plot_hex_values <- !is.null(proj$output$hex_values) 
  }
  
  # determine plotting values
  if (plot_hex_values) {
    add_legend <- TRUE
    if(plot_data_values==FALSE){
      y <- proj$output$hex_values
      legend_name="z-score"
    } else {
      y <- proj$output$hex_values2
      legend_name="Value"
    }
  } else {
    add_legend <- FALSE
    y <- rep(0, length(proj$map$hex))
    plot_significance <- FALSE
  }
  
  # replace plotting values with NA if below minimum hex coverage
  y[proj$output$hex_coverage < min_hex_coverage] <- NA
  
  # produce basic plot
  if (is.null(base_plot)) {
    plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
  }
  
  # add hexs
  if (plot_hex_grid) {
    plot1 <- plot1 + geom_sf(aes_(fill = ~y), color = NA, data = proj$map$hex)
  }
  
  # outline significance
  if (plot_significance) {
    
    # get significant hexes
    hex_signif <- get_significant_hexes(proj,
                                        empirical_tail = empirical_tail,
                                        alpha_raw = alpha_raw,
                                        min_hex_coverage = min_hex_coverage)
    
    # outline low values
    w <- hex_signif$which_lower
    if (length(w) != 0) {
      merged_poly_lower <- get_merged_poly(proj$map$hex[w], d = proj$map$hex_width/10)
      plot1 <- plot1 + geom_sf(color = "white", fill = NA, data = merged_poly_lower)
    }
    
    w <- hex_signif$which_upper
    if (length(w) != 0) {
      merged_poly_upper <- get_merged_poly(proj$map$hex[w], d = proj$map$hex_width/10)
      plot1 <- plot1 + geom_sf(color = "black", fill = NA, data = merged_poly_upper)
    }
    
  }
  
  # add points
  if (plot_sampling_points) {
    plot1 <- plot1 + geom_point(aes_(x = ~long, y = ~lat),
                                shape = 21, color = "white", fill = "black", size = 1,
                                data = proj$data$coords)
  }
  
  # add labeled points
  if (is.null(labeled_points)==FALSE){
    plot1 <- plot1 + geom_point(aes_(x=~x, y=~y),data=labeled_points,
                                shape=21, color="black", fill="white", size=1)
    plot1 <- plot1 + geom_text(aes(x=~x, y=~y,label=~label),data=labeled_points,
                               color="white",size=3,hjust=0, vjust=0) 
  }
  
  # titles and legends
  if (add_legend) {
    plot1 <- plot1 + scale_fill_gradientn(colours = col_scale, name = legend_name)
  } else {
    plot1 <- plot1 + scale_fill_gradientn(colours = NA)
    plot1 <- plot1 + guides(fill = FALSE)
  }
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
  # add barrier polygons
  if (nb > 0) {
    for (i in 1:nb) {
      plot1 <- plot1 + geom_polygon(aes_(x = ~long, y = ~lat),
                                    col = "white", fill = NA, linetype = "dashed",
                                    data = as.data.frame(poly_list[[i]]))
    }
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Interactive hex map of PlasmoMAPI output
#'
#' @description Interactive hex map of PlasmoMAPI output.
#'
#' @param proj object of class \code{pm_project}.
#' @param fill_opacity opacity of fill; 0 = fully transparent, 1 = fully opaque.
#' @param legend_opacity opacity of legend; 0 = fully transparent, 1 = fully
#'   opaque.
#' @param map_type an index from 1 to 137 indicating the type of base map. The
#'   map types are taken from \code{leaflet::providers}, see
#'   \href{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}{here}
#'   for an interactive gallary.
#' @param col_scale the colour scale to use.
#' 
#' @import leaflet
#' @importFrom viridisLite magma
#' @importFrom grDevices colorRamp
#' @importFrom methods slot
#' @export

plot_leaflet <- function(proj, fill_opacity = 0.8, legend_opacity = 1,
                         map_type = 110, col_scale = viridisLite::magma(100)) {
  
  # check inputs
  assert_custom_class(proj, "pm_project")
  assert_single_bounded(fill_opacity)
  assert_single_bounded(legend_opacity)
  assert_single_pos_int(map_type)
  
  # produce default colours
  add_legend <- TRUE
  #if ("hex_values" %in% names(proj$output)) {
  #  x <- proj$output$hex_values
  #} else {
  #  add_legend <- FALSE
  #  x <- rep(0, length(proj$map$hex))
  #}
  
  # add colour variable to hex map
  hex_values <- proj$output$hex_values
  hex_sf <- st_sf(geom = proj$map$hex, col = hex_values)
  
  # define colour ramp and palette
  col_ramp <- colorRamp(col_scale)
  if (all(hex_values == 0)) {
    col_ramp <- colorRamp(grey(0.5))
  }
  pal <- colorNumeric(col_ramp, domain = hex_values)
  
  # produce basic leaflet plot
  plot1 <- leaflet(hex_sf)
  plot1 <- addProviderTiles(plot1, leaflet::providers[[map_type]])
  
  # add hex polygons
  plot1 <- addPolygons(plot1, color = NA, fillColor = ~pal(col), fillOpacity = fill_opacity)
  
  # add legend
  if (add_legend) {
    plot1 <- addLegend(plot1, position = "bottomright", pal = pal, values = ~col,
                       title = "z-score", opacity = legend_opacity)
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Add points to dynamic map
#'
#' @description Add points to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_leaflet()} function.
#' @param lon,lat longitude and latitude of points.
#' @param col colour of points.
#' @param size size of points.
#' @param opacity opacity of points.
#'
#' @import leaflet
#' @export

overlay_points <- function(myplot, lon, lat, col = "black", size = 2, opacity = 1.0) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_vector_numeric(lon)
  assert_vector_numeric(lat)
  assert_same_length(lon, lat)
  assert_single_string(col)
  assert_single_pos(size, zero_allowed = FALSE)
  assert_single_pos(opacity, zero_allowed = TRUE)
  assert_bounded(opacity)
  
  # add circle markers
  myplot <- addCircleMarkers(myplot, lng = lon, lat = lat, radius = size,
                             fillColor = col, stroke = FALSE, fillOpacity = opacity)
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Plot daily counts of each host state from simulation
#'
#' @description For a set of simulation output produced by the function
#'   \code{sim_falciparum()}, plots daily simulated values in a given deme.
#'
#' @param x object of class \code{pm_sim}.
#' @param deme which deme to plot.
#' @param states which states to plot. Can be any subset of \code{c("S", "E",
#'   "I", "Sv", "Ev", "Iv")}.
#'
#' @importFrom grDevices grey
#' @import tidyr
#' @export

plot_daily_states <- function(x, deme = 1, states = c("S", "E", "I")) {
  
  # check inputs
  assert_custom_class(x, "pm_sim")
  assert_leq(deme, length(x$daily_values), message = "deme not found within simulation output")
  assert_vector(states)
  assert_in(states, c("S", "E", "I", "Sv", "Ev", "Iv"))
  
  # subset to desired rows and columns
  df_wide <- x$daily_values[[deme]][, c("time", states), drop = FALSE]
  
  # get to long format
  df_long <- tidyr::gather(df_wide, states, factor_key = TRUE)
  
  # choose plotting colours
  raw_cols <- daily_cols()
  plot_cols <- c(S = grey(0.5), E = raw_cols[2], I = raw_cols[1],
                 Sv = grey(0.8), Ev = raw_cols[6], Iv = raw_cols[7],
                 EIR = grey(0.0))
  
  # produce plot
  ggplot2::ggplot(df_long) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes_(x = ~time, y = ~count, color = ~state)) +
    ggplot2::scale_color_manual(values = plot_cols) +
    ggplot2::xlab("time (days)") + ggplot2::ggtitle(sprintf("Deme %s", deme))
}
