#################################################################################################################
# 
# Orientation analysis - written by Shannon P. McPherron
#
# Please cite McPherron (2018) from PLOS One.
#
#################################################################################################################
#
# The following is intended to take a data frame where each line is a shot, and convert it into
# a dataframe where each line is both shots.  It uses a counter field (starting with 0 and normally
# called suffix) to re-arrange the data.
#
generalized_prepare_dataset = function(xyz, ID = NA, SUFFIX = NA) {

  if (sum(c('X','Y','Z') %in% names(xyz)) != 3) {
    stop('There must be a field X, Y, and Z.')
  }
  
  xyz = unique(xyz)
  
  if (is.na(ID)) {
    if (!('UNIT' %in% names(xyz)) | !('ID' %in% names(xyz))) {
      stop('If ID data are not provided, then UNIT and ID columns must be present.')
    }
    ID = paste(xyz$UNIT,xyz$ID)
  }
  
  twoshots = table(ID)==2
  
  twoshots = twoshots[twoshots==TRUE]
  
  twoshots = rownames(twoshots)
  
  if (is.na(SUFFIX)) {
    if (!('SUFFIX') %in% names(xyz)) {
      stop('If SUFFIX data are not provided, then a column named SUFFIX must be present.')
    }
    SUFFIX = xyz$SUFFIX
  }

  xyz_twoshots = xyz[ID %in% twoshots & SUFFIX==0,]
  
  shot_two = xyz[ID %in% twoshots & SUFFIX==1, c('X','Y','Z')]
  
  colnames(shot_two) = c('X2','Y2','Z2')
  
  colnames(xyz_twoshots)[colnames(xyz_twoshots)=='X'] = 'X1'
  colnames(xyz_twoshots)[colnames(xyz_twoshots)=='Y'] = 'Y1'
  colnames(xyz_twoshots)[colnames(xyz_twoshots)=='Z'] = 'Z1'
  
  return(cbind(xyz_twoshots,shot_two))
  
}

#################################################################################################################
# Convert radians to degrees
#
deg = function (radian) { (radian * 180) / pi }

#################################################################################################################
#  Add a list containing XY pairs of line segments to a plot
#
add_overlays = function(overlays, col = "black", lty = 1) { for (overlay in overlays) lines(overlay$x, overlay$y, col = col, lty = lty) }

#################################################################################################################
#  Compute limits based on either the data or the data+overlays
#
set_limits = function(xyz, overlays = NULL, limits = NULL) {
	if (!is.null(limits)) {
		if (limits[1] == "data") {
			limits = c(floor(min(xyz$X1,xyz$X2)), floor(max(xyz$X1,xyz$X2) + 1),
			           floor(min(xyz$Y1,xyz$Y2)), floor(max(xyz$Y1,xyz$Y2) + 1))
		} else {
			if (limits[1] == "dataplus") {
				limits = c(min(xyz$X1,xyz$X2), max(xyz$X1,xyz$X2),
				           min(xyz$Y1,xyz$Y2), max(xyz$Y1,xyz$Y2))
				for (overlay in overlays) {
					limits[1] = floor(min(limits[1], overlay$x))
					limits[2] = floor(max(limits[2], overlay$x) + 1)
					limits[3] = floor(min(limits[3], overlay$y))
					limits[4] = floor(max(limits[4], overlay$y) + 1) } } } } 
	return(limits) }

#################################################################################################################
#  Basic bearing and plunge analysis
#
orientations = function(xyz, level = "All points", min_sample = 50, 
                        main = "", limits = NULL, overlays = NULL, spatial_benn = FALSE, 
                        spatial_bearing = FALSE, p = .05, background_image = NULL,...) {
 
	limits = set_limits(xyz, overlays, limits)
	
	results = circular_statistics(xyz, level, min_sample)
	
	if (main=="") {
	  titles = rownames(results)
	} else {
	  titles = rep(main, nrow(results))
	}
	
	# Go through now and print figures for each level
	for (i in 1:nrow(results)) {

		if (nrow(results)==1) xyz_level = xyz else xyz_level = subset(xyz, level==rownames(results)[i])
			
		# Build a graphic that consists of the two shots plotted in the site
		# plus the plunge rose diagram, and the bearing schmidt diagram
		titletext  = titles[i]
		
		# If there are at least min_sample objects, then color code by Benn values
		if ((nrow(xyz_level) >= min_sample) & (spatial_benn)) {
		  spatial_benn_plots(xyz_level, titletext, limits, overlays,
		                    background_image = background_image,...) 
		  t_mai = par(c("mai"))
      par(mai = c(.05, .05, .05, .05))
		  rose_diagram(plunge_and_bearing(xyz_level)$bearing,
		               main = "Bearing", cex = .8)
		  schmidt_diagram(angles = plunge_and_bearing(xyz_level), redraw = FALSE)
		  rose_diagram_plunge(plunge_and_bearing(xyz_level)$plunge,
		                      main = "Plunge", cex = .8)
		  par(mai = t_mai) }
		else {
      plunge_and_bearing_plots(xyz_level, titletext, limits,
                               overlays, background_image,
                               spatial_benn = spatial_benn)
    }
	}

	return(results)
}

#################################################################################################################
# Create basic figures for bearing and plunge information
#
plunge_and_bearing_plots = function(xyz, main = "", limits = NULL,
                                    overlays = NULL, background_image = NULL,
                                    spatial_benn = FALSE) {

	t_mai = par(c("mai"))
	par(mai = c(.5, .5, .5, .05))
	plot_2shot(xyz, main = main, limits = limits, background_image = background_image) 
	if (!is.null(overlays)) add_overlays(overlays)
	if (spatial_benn) frame()
	par(mai = c(.05, .05, .05, .05))
	rose_diagram(plunge_and_bearing(xyz)$bearing, main = "Bearing", cex = .8)
	schmidt_diagram(angles = plunge_and_bearing(xyz), redraw=FALSE)
	rose_diagram_plunge(plunge_and_bearing(xyz)$plunge, main = "Plunge", cex = .8)
	par(mai = t_mai)
		
}

#################################################################################################################
# Rose.diag from CircStats package adapted to do only plunge angle from 0 to 90,
# to deal with angles that turn clockwise, and some other customizations.
#
rose_diagram_plunge = function (x, bins = 10, main = "", prop = 1, 
                                cex = 1, pch = 19, pts_on_edge = FALSE,
                                color_codes = NULL, color_filled = NULL,
                                pnt_col = 'black', bar_col = 'white', bg = 'white',
                                dotsep = 40, shrink = 1,...) {
	x = rad(x)
  x <- x%%(2 * pi)
	plot(cos(seq(3/4 * 2 * pi, 2 * pi, length = 1000)),
	     sin(seq(3/4 * 2 * pi, 2 * pi, length = 1000)),
	     axes = FALSE, xlab = "", ylab = "", main = "", type = "l",
	     xlim = shrink * c(-.15, 1.15), ylim = shrink * c(-1.15, 0.15), asp = 1)
  polygon(c(0,cos(seq(3/4 * 2 * pi, 2 * pi, length = 1000))),
          c(0,sin(seq(3/4 * 2 * pi, 2 * pi, length = 1000))), col = bg)
  #text(-.88,1.05,main,cex=1.1)
	text(-.03, .1, main, cex = cex * 1.1)
	lines(c(0, 0), c(-0.9, -1))
	text(0.005, -1.075, "90", cex = cex * 1.1)
	lines(c(0.9, 1), c(0, 0))
	text(1.05, 0, "0", cex = cex * 1.1)
	lines(c(0,0),c(-1,0))
	n <- length(x)
	freq <- c(1:bins)
	arc <- (1/2 * pi)/bins
	for (i in 1:bins) {
		newi = bins - i + 1 									# This turns the angles clockwise
		freq[i] <- sum(x <= newi * arc & x > (newi - 1) * arc)
	}
	rel.freq <- freq/n
	radius <- sqrt(rel.freq) * prop
	#radius <- freq/max(freq) * prop     					  # This will bring bars to circle but with flat proportions
	radius <- radius/max(radius) * prop  						# This will bring bars to circle but with exponential proportions
	sector <- seq(0, 1/2 * pi - (1/2 * pi)/bins, length = bins)
	sector <- sector - (pi/2) 									# This rotates them to the right spot
	mids <- seq(arc/2, 1/2 * pi - (pi/4)/bins, length = bins)
	index <- cex/dotsep
	for (i in 1:bins) {
		if (rel.freq[i] != 0) {
			xp = c(0, radius[i] * cos(sector[i]), radius[i] * cos(sector[i] + (1/2 * pi)/bins))
			yp = c(0, radius[i] * sin(sector[i]), radius[i] * sin(sector[i] +  (1/2 * pi)/bins))
			polygon(xp, yp, col = bar_col,...)
			
			if (pts_on_edge) {
			  xp = cos(x)
			  yp = -sin(x)
			  if (!is.null(color_codes)) {
			    points(xp[color_filled], yp[color_filled], cex = (cex * .8),
			           col=color_codes[color_filled], pch = 19)	
			    points(xp[!color_filled], yp[!color_filled], cex = (cex * .8),
			           col=color_codes[!color_filled], pch = 21)	
			  } else {	
			    points(xp, yp, cex= (cex * .8), pch = pch, col = pnt_col) }
			}	} }
}

#################################################################################################################
# Modifed version of the Watson.Two from CircStats package to do two things:
#	1) Return the t-statistic
#	2) Return False or True 
#
watson_two = function (x, y, alpha = 0, plot = FALSE) {
	n1 <- length(x)
	n2 <- length(y)
	n <- n1 + n2
	if (n < 18) 
		cat("Total Sample Size < 18:  Consult tabulated critical values", "\n", "\n")
	if (plot == TRUE) {
		x <- sort(x%%(2 * pi))
		y <- sort(y%%(2 * pi))
		plot.edf(x, main = "Comparison of Empirical CDFs", xlab = "", ylab = "")
		par(new = TRUE)
		plot.edf(y, xlab = "", ylab = "", axes = FALSE, lty = 2)
	}
	#cat("\n", "      Watson's Two-Sample Test of Homogeneity", "\n", "\n")
	x <- cbind(sort(x%%(2 * pi)), rep(1, n1))
	y <- cbind(sort(y%%(2 * pi)), rep(2, n2))
	xx <- rbind(x, y)
	rank <- order(xx[, 1])
	xx <- cbind(xx[rank, ], seq(1:n))
	a <- c(1:n)
	b <- c(1:n)
	for (i in 1:n) {
		a[i] <- sum(xx[1:i, 2] == 1)
		b[i] <- sum(xx[1:i, 2] == 2)
	}
	d <- b/n2 - a/n1
	dbar <- mean(d)
	u2 <- (n1 * n2)/n^2 * sum((d - dbar)^2)
	crits <- c(99, 0.385, 0.268, 0.187, 0.152)
	#cat("Test Statistic:", round(u2, 4), "\n")
	if (sum(alpha == c(0, 0.001, 0.01, 0.05, 0.1)) == 0) 
		stop("Invalid input for alpha")
	else if (alpha == 0) {
		if (u2 > 0.385) 
			cat("P-value < 0.001", "\n", "\n")
		else if (u2 > 0.268) 
			cat("0.001 < P-value < 0.01", "\n", "\n")
		else if (u2 > 0.187) 
			cat("0.01 < P-value < 0.05", "\n", "\n")
		else if (u2 > 0.152) 
			cat("0.05 < P-value < 0.10", "\n", "\n")
		else cat("P-value > 0.10", "\n", "\n")
	}
	else {
		index <- (1:5)[alpha == c(0, 0.001, 0.01, 0.05, 0.1)]
		Critical <- crits[index]
		if (u2 > Critical) {
			Reject <- "Reject Null Hypothesis"
			reject_truefalse=TRUE }
		else {
			Reject <- "Do Not Reject Null Hypothesis"
			reject_truefalse=FALSE }
		#cat("Level", alpha, "Critical Value:", round(Critical, 4), "\n")
		#cat(Reject, "\n", "\n")
	}
	return(list(t.statistic=round(u2,4),reject=reject_truefalse))
}

#################################################################################################################
# Function to plot artifacts used in orientation analysis
#
plot_2shot = function(xyz, level = "All Points", lty = 3, lwd = .5, fg = "black",
                      bg = "white", cex = 1, limits = NULL, draw_grid = TRUE,
                      color_codes = NULL, color_filled = NULL, col = 'black',
                      background_image = NULL,...) {
	
	#par(mai=c(.5,.5,.5,.5))

  level = factor(level)
  
  for (l in levels(level)) {

    xyz_level = subset(xyz, level == l)

    # Either prepare a plot to fit the points or to the specified limits
    if (length(levels(level))>1 & !exists('main')) {
      if (!is.null(limits)) {
        plot(xyz_level$X1, xyz_level$Y1, asp=1, type="n", xlab="X", ylab="Y",
             xlim=c(limits[1],limits[2]), ylim=c(limits[3],limits[4]), main = l,...) 
      } else {
        plot(xyz_level$X1, xyz_level$Y1, asp=1, type="n", xlab="X", ylab="Y", main = l,...) }
    } else {
    	if (!is.null(limits)) {
    		plot(xyz_level$X1, xyz_level$Y1, asp=1, type="n", xlab="X", ylab="Y",
    		     xlim=c(limits[1],limits[2]), ylim=c(limits[3],limits[4]),...) 
    	} else {
    		plot(xyz_level$X1, xyz_level$Y1, asp=1, type="n", xlab="X", ylab="Y",...) } }
  
  	# Make the background black if color coding is being used and there is no background image
  	if (!is.null(color_codes) & is.null(background_image)) rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 'grey50')
  	
    # If a georeferenced background image (geotiff) is specified, get the limits and then display it
    if (!is.null(background_image)) {
      library(raster)
      library(tiff)
      library(rgdal)
      r_image = raster(background_image)
      rasterImage(suppressWarnings(readTIFF(background_image,convert=TRUE)), r_image@extent@xmin, r_image@extent@ymin, r_image@extent@xmax, r_image@extent@ymax)
    }
    
  	# Draw a meter grid
  	if (draw_grid) {
  		for (x in floor(par("usr")[1]):floor(par("usr")[2])) {
  			segments(x,par("usr")[3],x, par("usr")[4], col = fg, lty = lty, lwd = lwd) }
  			
  		for (y in floor(par("usr")[3]):floor(par("usr")[4])) {
  			segments(par("usr")[1],y, par("usr")[2],y, col = fg, lty = lty, lwd = lwd) } }
  			
  	# Plot the points color coded
  	if (!is.null(color_codes)) {
  		#points(xyz$X1,xyz$Y1,pch=19,cex=.5,col=xyz$point_color) 
  		if (is.null(color_filled)) {
  			points(xyz_level$X1, xyz_level$Y1, cex = (cex * .8), col = color_codes, pch = 19)	
  			#points(xyz_level$X1[which(xyz_level$Z1>xyz_level$Z2)], xyz_level$Y1[which(xyz_level$Z1>xyz_level$Z2)], pch=19, cex = (cex * .8),col = color_codes)		
  			#points(xyz_level$X2[which(xyz_level$Z2>xyz_level$Z1)], xyz_level$Y2[which(xyz_level$Z2>xyz_level$Z1)], pch=19, cex = (cex * .8),col = color_codes)		
  			#segments(xyz_level$X1, xyz_level$Y1, xyz_level$X2, xyz_level$Y2, xaxt="n")
  		} else {
  			points(xyz_level$X1[color_filled], xyz_level$Y1[color_filled], cex = (cex * .8),
  			       col = color_codes[color_filled], pch = 19)	
  			points(xyz_level$X1[!color_filled], xyz_level$Y1[!color_filled], cex = (cex * .8),
  			       col = color_codes[!color_filled], pch = 21)	
  		}
  	} else {
  		segments(xyz_level$X1, xyz_level$Y1, xyz_level$X2, xyz_level$Y2,
  		         xaxt="n", col = col)
  		points(xyz_level$X1[which(xyz_level$Z1>xyz_level$Z2)],
  		       xyz_level$Y1[which(xyz_level$Z1>xyz_level$Z2)],
  		       pch = 19, cex = (cex * .25), col = col)		
  		points(xyz_level$X2[which(xyz_level$Z2>xyz_level$Z1)],
  		       xyz_level$Y2[which(xyz_level$Z2>xyz_level$Z1)],
  		       pch = 19, cex = (cex * .25), col = col)		
  	}
  }		
	#par(mar=c(5, 4, 4, 2) + 0.5)
	
}


#################################################################################################################
# Internal function to setup a circular graph for rose or schmidt diagrams
#
draw_circle_diagram = function(bg = "white", cex = 1, main = "",...) {

  plot(cos(seq(0, 2 * pi, length = 1000)), sin(seq(0, 2 * pi, length = 1000)), axes = FALSE,
       xlab = "", ylab = "", main = "", type = "n", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)  
  polygon(cos(seq(0, 2 * pi, length = 1000)), sin(seq(0, 2 * pi, length = 1000)), col = bg)
  
  text(-.88, 1.05, main, cex = (cex * 1.1))
  lines(c(0, 0), c(0.9, 1))
  text(0.005, 1.10, "0", cex = (cex * 1.1))
  lines(c(0, 0), c(-0.9, -1))
  text(0.005, -1.11, "180", cex = (cex * 1.1))
  lines(c(-1, -0.9), c(0, 0))
  text(-1.16, 0, "270", cex = (cex * 1.1))
  lines(c(0.9, 1), c(0, 0))
  text(1.13, 0, "90", cex = (cex * 1.1))
  lines(c(-.05, .05), c(0, 0))
  lines(c(0, 0), c(-0.05, .05))
}

#################################################################################################################
# Function to plot points on circle in Schmidt equal area space
# Code to draw circle taken from rose.diag in CircStats package
# Angles should be in decimal degrees
#
schmidt_diagram = function(bearing = NULL, plunge = NULL, angles = NULL, 
                           level = "All Points", color_codes = NULL,
                           color_filled = FALSE, redraw = TRUE, col = "black",
                           pch = 19, cex = 1, main = "",...) {

  library(CircStats)
  
  if (is.null(angles) & (is.null(bearing) | is.null(plunge))) stop("Not enough data passed to schmidt.diagram.2shot. Bearing and plunge angles are required.")

  if (!is.null(angles)) {
    plunge = angles[,1]
    bearing = angles[,2]
  }
  
  level = factor(level)
  
  for (l in levels(level)) {

    bearing_angle_level = subset(bearing, level == l)
    plunge_angle_level = subset(plunge, level == l)
    
    # If the circle doesn't already exist (from doing a Rose diagram), then make it
  	if (length(levels(level))>1) {
      if (!exists('main')) draw_circle_diagram(main = l,...) else draw_circle_diagram(main = main[which(levels(level)==l)],...)
  	} else {  
      if (redraw) draw_circle_diagram(...)
  	}
    
  	# Shift points into Schmidt space 
  	d = sin(rad((90 - plunge_angle_level) / 2)) / sin(rad(45))
  	x = d * sin(rad(bearing_angle_level))
  	y = d * cos(rad(bearing_angle_level))
  
  	# Plot them using color coding or not
  	if (!is.null(color_codes)) {
  		points(x[color_filled],y[color_filled], cex = (cex * .8),
  		       col = color_codes[color_filled], pch = 19)	
  		points(x[!color_filled],y[!color_filled], cex = (cex * .8),
  		       col = color_codes[!color_filled], pch = 21)	
  	} else {	
  		points(x, y, cex = (cex * .8), pch = pch, col = col) }
  }	
}

#################################################################################################################
# Rose.diag from CircStats package adapted for orientation analysis.
#
rose_diagram = function(bearing = NULL, plunge = NULL, angles = NULL,
                        bins = 36, level = "All points", main = "",
                        color_codes = NULL, color_filled = NULL, prop = 1,
                        pts_on_edge = FALSE, pts_schmidt = FALSE,
                        pch = 16, dotsep = 40, cex = 1,
                        bg = "white", bar_col = "white", pnt_col = "black",...) {
  
  if (is.null(angles) & is.null(bearing)) stop("No data passed to rose_diagram.  Need bearing angles.")

  if (!is.null(angles)) {
    plunge = angles[,1]
    bearing = angles[,2]
  }
  
  #bearing_angle <- bearing_angle%%(2 * pi)
  bearing <- bearing%%(360)
  
  level = factor(level)
  
  for (l in levels(level)) {
    
    draw_circle_diagram(bg = bg, cex = cex, main = main[which(levels(level)==l)],...)

    bearing_angle_level = rad(subset(bearing, level == l))
    
    n <- length(bearing_angle_level)
		freq <- c(1:bins)
		arc <- (2 * pi) / bins
		for (i in 1:bins) {
			newi = bins - i + 1 									# This turns the angles clockwise
			freq[i] <- sum(bearing_angle_level <= newi * arc & bearing_angle_level > (newi - 1) * arc)
		}
		rel.freq <- freq / n
		radius <- sqrt(rel.freq) * prop
		#radius <- freq/max(freq) * prop     					# This will bring bars to circle but with flat proportions
		radius <- radius/max(radius) * prop  						# This will bring bars to circle but with exponential proportions
		sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
		sector <- sector + (pi/2) 									# This shifts them into 0 degrees up
		mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
		index <- cex/dotsep
		for (i in 1:bins) {
			if (rel.freq[i] != 0) {
			  x = c(0, radius[i] * cos(sector[i]), radius[i] *  cos(sector[i] + (2 * pi)/bins))
			  y = c(0, radius[i] * sin(sector[i]), radius[i] * sin(sector[i] + (2 * pi)/bins))
			  polygon(x, y, col = bar_col,...) } }
    if (pts_schmidt) schmidt_diagram(bearing = subset(bearing, level == l),
                                     plunge = subset(plunge, level == l),
                                     color_codes = color_codes,
                                     color_filled = color_filled,
                                     col = pnt_col,
                                     redraw = FALSE, cex = cex,...) 
		if (pts_on_edge) {
		  x = sin(bearing_angle_level)
		  y = cos(bearing_angle_level)
		  if (!is.null(color_codes)) {
		    points(x[color_filled], y[color_filled], cex = (cex * .8),
		           col=color_codes[color_filled], pch = 19)	
		    points(x[!color_filled], y[!color_filled], cex = (cex * .8),
		           col=color_codes[!color_filled], pch = 21)	
		  } else {	
		    points(x, y, cex= (cex * .8), pch = pch, col = pnt_col) }
		}
  }
}

#################################################################################################################
# Compute plunge and bearing angles from XYZ 2-shot data
#
plunge_and_bearing = function(xyz) {

	# Compute the plunge angle
	# Note: Angles are converted to degrees
	run = sqrt((xyz$X2 - xyz$X1) ^ 2 + (xyz$Y2 - xyz$Y1) ^ 2)
	rise = abs(xyz$Z2 - xyz$Z1)
	plunge_angle = ifelse(run == 0, 90, deg(atan(rise / run)))

	# Compute Schmidt bearing (lower hemisphere)
	# Note: Angle are adjusted to be clockwise with north (0) = positive y-axis
	# Note: Angles are converted to degrees.
	run = ifelse(xyz$Z1 >= xyz$Z2, xyz$X2 - xyz$X1, xyz$X1 - xyz$X2)
	rise = ifelse(xyz$Z1 >= xyz$Z2, xyz$Y2 - xyz$Y1, xyz$Y1 - xyz$Y2)
	slope = ifelse(run != 0, rise / run, 1000000)
	bearing_angle = 90 - deg(atan(slope))
	bearing_angle[run <= 0] = 180 + bearing_angle[run <= 0]
	
	return(data.frame(plunge = plunge_angle, bearing = bearing_angle))   }

#################################################################################################################
# Return circular stats for XYZ object
#
circular_statistics = function(data, level = "All observations", min_sample = 50) {

  library(CircStats)
  
	level = factor(level)

	# Create a data frames for the results
	if (sum(c('X1','X2','Y1','Y2','Z1','Z2') %in% names(data)) == 6) {
    xyz = TRUE
	  results = matrix(NA, nrow=length(levels(level)), ncol=11, 
	                   dimnames = list(levels(level),c("N","Length_Mean","Length_SD","Bearing_L","Bearing_Mean",
	                                                   "Bearing_Var","Bearing_P","Plunge_L","Plunge_Mean","Plunge_Var","Plunge_P"))) }
	else {
	  xyz = FALSE
	  results = matrix(NA, nrow=length(levels(level)), ncol=9, 
	                   dimnames = list(levels(level),c("N","Bearing_L","Bearing_Mean","Bearing_Var","Bearing_P",
	                                                   "Plunge_L","Plunge_Mean","Plunge_Var","Plunge_P"))) }

	for (a_level in levels(level)) {

		data_level = subset(data, level==a_level)
		
		# If there are to few cases, don't try to do orientations
		if (nrow(data_level) < min_sample) {
			results[a_level,"N"] = nrow(data_level)
			results[a_level,2:ncol(results)] = 0 
		} else {
			
		  # If XYZ data, then compute average observation length and SD
		  if (xyz) {
			  l = compute_lengths(data_level)
			  vector_length = round(mean(l), 3) 
			  vector_length_sd = round(sd(l), 3) 
  		  pb = plunge_and_bearing(data_level) }
		  else {
        vector_length = NA
        vector_length_sd = NA
		    pb = data_level
		  }

			# Circle stats on bearing
			Bearing_L = round(r.test(pb$bearing,degree=TRUE)$r.bar, 2)
			Bearing_Mean = round(deg(circ.mean(rad(pb$bearing))), 1) %% 360
			Bearing_Var = round(deg(circ.disp(rad(pb$bearing))$var), 1)
			Bearing_P = round(r.test(pb$bearing,degree=TRUE)$p.value, 2)

			# Stats on plunge (circle stats are not used because plunge is not circular - it varies from 0 to 90 and 0<>90)
			# Kolmogorov-Smirnov test is used to test for non-uniform distribution
			Plunge_L = round(r.test(pb$plunge,degree=TRUE)$r.bar,2)
			Plunge_Mean = round(mean(pb$plunge),1)
			Plunge_Var = round(sd(pb$plunge),1)
			while (any(duplicated(pb$plunge))) {
			  # This code is necessary because ks.test gives a warning when there are ties
				pb$plunge[duplicated(pb$plunge)]=pb$plunge[duplicated(pb$plunge)] + .000001 }
			Plunge_P = round(ks.test(pb$plunge,punif,min=0,max=90)$p.value, 2)  
			
      if (xyz) {
			  results[a_level,] = c(nrow(data_level), vector_length, vector_length_sd, Bearing_L,
			                        Bearing_Mean, Bearing_Var, Bearing_P, Plunge_L, Plunge_Mean,
			                        Plunge_Var, Plunge_P) }
		  else {
		    results[a_level,] = c(nrow(data_level), Bearing_L, Bearing_Mean, Bearing_Var, Bearing_P,
		                          Plunge_L, Plunge_Mean, Plunge_Var, Plunge_P) }
			} }
	
	return(results) }
			
#################################################################################################################
# Calc. color codes based on average bearing of nearest objects
#
spatial_bearing = function(xyz, nearest = 40) {
					
	# Make a place to hold the computed mean bearings of nearest neighbors
	near_avg_bearing = vector(mode = "numeric", length = nrow(xyz))
	near_avg_plunge = vector(mode = "numeric", length = nrow(xyz))
	near_avg_bearing_p = vector(mode = "numeric", length = nrow(xyz))
	
	# Go through each artifact, get the nearest (default is 40), and compute mean bearing angle
	for (k in 1:nrow(xyz)) {
		centerx = xyz$X1[k]
		centery = xyz$Y1[k]
		d = sqrt((centerx-xyz$X1)^2 + (centery-xyz$Y1)^2)
		xyz_subsample = xyz[order(d)[1:nearest],]

		# Get the mean bearing angle, test significance, and mean plunge angle of this subset of artifacts
		near_avg_bearing[k] = round(deg(circ.mean(rad(plunge_and_bearing(xyz_subsample)$bearing))), 1) 
		near_avg_bearing_p[k] = round(r.test(plunge_and_bearing(xyz_subsample)$bearing,degree=TRUE)$p.value, 2)
		near_avg_plunge[k] = round(mean(plunge_and_bearing(xyz_subsample)$plunge), 1) }

	# Color coding based on average bearing and on average plunge (higher plunge angles are more less saturated - i.e. more white)
	point_color = hsv((near_avg_bearing %% 360) / 360, 1-(near_avg_plunge %% 90) / 90, 1)
	return(data.frame(point_color = point_color, p = near_avg_bearing_p)) }

#################################################################################################################
# Create figures with color coding based on bearing of artifacts in immediate vicinity
# Significant probabilities are plotted as filled circles
#
spatial_bearing_plots = function(xyz, main = "", limits = NULL, overlays = NULL, p = .05, background_image = NULL) {
	t_mai = par(c("mai"))
	par(mai = c(.5, .5, .5, .05))
	coloring = spatial_bearing(xyz)
	plot_2shot(xyz, main = paste(main," - Color by Avg. Bearing"),
	           limits = limits,
	           color_codes = coloring$point_color,
	           color_filled = (coloring$p <= .05),
	           background_image = background_image) 
	if (!is.null(overlays)) add_overlays(overlays, col = "white")
	par(mai = c(.05,.05,.05,.05))
	rose_diagram(plunge_and_bearing(xyz)$bearing, main="Bearing", colorcoded = TRUE)
	schmidt_diagram(plunge_and_bearing(xyz)$bearing, plunge_and_bearing(xyz)$plunge,
	                color_codes = rep("black",nrow(xyz)), color_filled = (coloring$p <= .05), redraw = FALSE)
	rose_diagram_plunge(plunge_and_bearing(xyz)$plunge, main="Plunge")
	par(mai=t_mai)
}

#################################################################################################################
# Calc. color codes based on Benn values of nearest objects
#
spatial_benn = function(xyz, nearest = 40, maximum_distance = NA) {
	
  library(fields)
  library(colorspace)
  
	benn = matrix(NA, nrow = nrow(xyz), ncol = 2, 
	              dimnames = list(rownames(xyz), c("elongation","isotropy")))
	
	# For each artifact, get the nearest artifacts and compute Benn values
	for (k in 1:nrow(xyz)) {
		centerx = xyz$X1[k]
		centery = xyz$Y1[k]
		d = sqrt((centerx-xyz$X1)^2 + (centery-xyz$Y1)^2)
		sorted_pos = order(d)
		if (!is.na(maximum_distance)) sorted_pos = sorted_pos[d[sorted_pos] <= maximum_distance]
		xyz_subsample = xyz[sorted_pos[1:nearest], ]
		benn [k,] = benn(xyz_subsample, min_sample = nearest)[,c("EL","IS")]
	}

	# Calculate where the points would fall on Benn diagram so colors can be assigned
	b = benn_coords(cbind(elongation = benn[,"elongation"], isotropy = benn[,"isotropy"]))
	xp = b[,1]
	yp = b[,2]
	
	# Assign colors
	color_range = two.colors(n=21, start="green", end="red", middle="yellow")
	point_color = color_range[xp * 20 + 1]
	rgb_color = hex2RGB(point_color)
	hsv_color = t(rgb2hsv(rgb_color@coords[,1], rgb_color@coords[,2], rgb_color@coords[,3]))
	point_color = hsv(hsv_color[,1],1-yp,1)
	
	return(list(point_color, benn)) }

#################################################################################################################
# Create figures with color coding based on Benn valuse artifacts in immediate vicinity
# Significant probabilities are plotted as filled circles
#
spatial_benn_plots = function(xyz,  main="", limits = NULL, overlays = NULL,
                              background_image = NULL,
                              nearest = 40,...) {
	
	color_codes = spatial_benn(xyz, nearest = nearest)
	
	# Prepare a plan view graph of the results
	t_mai = par(c("mai"))
	par(mai = c(.5, .5, .5, .05))
  plot_2shot(xyz, main = main,
             limits = limits,
             color_code = color_codes[[1]],
             background_image = background_image)
	points(xyz$X1[color_codes[[2]][,2]>.4], xyz$Y1[color_codes[[2]][,2]>.4],
	       pch = 21, col = "black", bg = color_codes[[1]][color_codes[[2]][,2]>.4])
	
	if (!is.null(overlays)) add_overlays(overlays, col = "white")

	# Show a Benn diagram color coded as well and with the individual artifact computed values
	par(mai=c(.05,.05,.05,.05))
	benn_diagram(x = color_codes[[2]], main = "Color Key", cex = .8,
	             border = "black", labels = "outside", bg = rgb(0,0,0),
	             draw_grid = TRUE, newpage = TRUE,
	             colorcoded = TRUE,...)
	par(mai=t_mai)
	
}
	
#################################################################################################################
# Helping function to place points on Benn diagram
# Expects an object with columns named elongation and isotropy
#
benn_coords = function(benn) 	return(cbind(X = benn[,"elongation"] + (.5 * benn[,"isotropy"]), Y = benn[,"isotropy"] * .866))

#################################################################################################################
# Return list containing pairs of samples that are significantly different
#
benn_permutations = function(xyz1, xyz2, resampling = 100, min_sample = 0) {
	xyz_resample= rbind(xyz1, xyz2)
	benn_1_base = benn(xyz1, min_sample = min_sample)
	benn_2_base = benn(xyz2, min_sample = min_sample)
	d_base = sqrt((benn_1_base[,"IS"] - benn_2_base[,"IS"])^2 + (benn_1_base[,"EL"] - benn_2_base[,"EL"])^2 )
	d = vector(mode = "numeric", length = resampling)
	for (m in 1:resampling) {
		xyz_resample$group = "1"
		xyz_resample$group[sample(nrow(xyz_resample), nrow(xyz2), replace = FALSE)] = "2"
		benn_1 = benn(subset(xyz_resample, group=="1"), min_sample = min_sample)
		benn_2 = benn(subset(xyz_resample, group=="2"), min_sample = min_sample)
		d[m] = sqrt((benn_1[,"IS"] - benn_2[,"IS"])^2 + (benn_1[,"EL"] - benn_2[,"EL"])^2 ) 	}
	return((resampling - sum(d_base > d)) / resampling) }

#################################################################################################################
# Do all pairwise comparisons and return a list with results and with segments
# ready to be added to a Benn diagram using the segments() function.
#
benn_permutations_by_level = function(xyz, level, min_sample = 30, resampling = 100, p = .05) {
  
  benn_indices = benn(xyz = xyz, level = factor(level), min_sample = min_sample)

  perm_results = matrix(NA, nrow = nrow(benn_indices), ncol = nrow(benn_indices),
                        dimnames = list(rownames(benn_indices), rownames(benn_indices)))
  line_segments = NULL
  
  levelnames = rownames(benn_indices)

  for (k in 1:(length(levelnames) - 1)) {
    for (l in (k+1):length(levelnames)) {
      if (!is.na(benn_indices[levelnames[k],"EL"]) & !is.na(benn_indices[levelnames[l],"EL"])) {
        perm_results[levelnames[k], levelnames[l]] =
          benn_permutations(subset(xyz, level == levelnames[k]), 
                            subset(xyz, level == levelnames[l]), resampling = resampling)
        if (perm_results[levelnames[k],levelnames[l]] > p) {
          x0y0 = benn_coords(cbind(elongation = benn_indices[levelnames[k],"EL"],
                                   isotropy = benn_indices[levelnames[k],"IS"]))
          x1y1 = benn_coords(cbind(elongation=benn_indices[levelnames[l],"EL"],
                                   isotropy = benn_indices[levelnames[l],"IS"]))
          line_segments = rbind(line_segments, c(x0y0,x1y1)) } } } }
  
  return(list(segments = line_segments, results = perm_results))
}
  
#################################################################################################################
# Return list containing contour of resampled probability of Benn values
# For an explanation of swaps see Ringrose 1996 and McPherron 2018 (PLOS One)
#
benn_resampling = function(xyz, resampling = 10000, p = .95,
                           min_sample = 30, remove_swaps = FALSE)  {

	resampling_matrix = matrix(data = 0, nrow = 101, ncol = 101)
  
	if (remove_swaps) e = eigen_values(vector_normals(xyz))
  
	swaps = 0
  
	for (k in 1:resampling) {
	  s = eigen_values(vector_normals(xyz[sample(nrow(xyz), replace = TRUE),]))
		if (remove_swaps) keep = (sum(ringrose_swap(e, s)==c(1,2,3))==3) else keep = TRUE
		if (keep) {
		  isotropy = s$values[3] / s$values[1]
		  elongation = 1 - (s$values[2] / s$values[1])
		
  		resampling_matrix[floor(elongation * 100) + 1, floor(isotropy * 100) + 1] = 
  		  resampling_matrix[floor(elongation * 100) + 1,floor(isotropy * 100) + 1] + 1 } 
		else {
		  swaps = swaps + 1 }
		}
	
  resampling = resampling - swaps
  
	for (k in 1:resampling) {
		sample_portion = sum(resampling_matrix[resampling_matrix > k])
		if (sample_portion < (p * resampling)) {
			sample_portion_prev = sum(resampling_matrix[resampling_matrix > (k - 1)])
			contour_interval = (k-1) + (sample_portion_prev - (p * resampling)) /
			  (sample_portion_prev - sample_portion)
			return(contourLines(x = seq(0, 1, length.out = nrow(resampling_matrix)),
			                    y = seq(0, 1, length.out = ncol(resampling_matrix)),
			                    resampling_matrix, levels = c(contour_interval)))	} } }
			
#################################################################################################################
#  Return Benn statistics for an XYZ object
#
benn = function(xyz, level = "All Points", min_sample = 40) {

	benn = matrix(nrow = length(unique(level)), ncol = 6, 
	              dimnames = list(unique(level),c("N","E1","E2","E3","IS","EL")))

	for (l in unique(level)) {

		xyz_level = subset(xyz, level==l)
		
		if (nrow(xyz_level) < min_sample) {
			benn[l,] = c(nrow(xyz_level), rep(NA,5))
		
		} else {
	
			# Normalize and compute eigen values 
			e = eigen_values(vector_normals(xyz_level))

			# Compute shape indices for Benn Diagram
			isotropy = e$values[3] / e$values[1]
			elongation = 1 - (e$values[2] / e$values[1])
	
			benn[l,] = c(nrow(xyz_level), e$values[1], e$values[2], e$values[3], isotropy, elongation) } } 

	return(benn) }

#################################################################################################################
#  Normalize vectors 
#
vector_normals = function(xyz) {

	l = compute_lengths(xyz)
	xnorm = ifelse(l != 0, (xyz$X1 - xyz$X2) / l, 0)
	ynorm = ifelse(l != 0, (xyz$Y1 - xyz$Y2) / l, 0)
	znorm = ifelse(l != 0, (xyz$Z1 - xyz$Z2) / l, 0) 
	
	return(cbind(xnorm, ynorm, znorm))  }

#################################################################################################################
#  Compute eigen values from normalized vectors
#
eigen_values = function(xyz) {

	l = xyz[, 1]	# X
	m = xyz[, 2]	# Y
	n = xyz[, 3]	# Z

	# Build a matrix prior to computing eigen values
	M11 = sum(l ^ 2)
	M12 = sum(l * m)
	M13 = sum(l * n)
	M21 = sum(m * l)
	M22 = sum(m ^ 2)
	M23 = sum(m * n)
	M31 = sum(n * l)
	M32 = sum(n * m)
	M33 = sum(n ^ 2)
	M = matrix(c(M11,M12,M13,M21,M22,M23,M31,M32,M33), nrow = 3, ncol = 3)

	# Compute eigen values on matrix normalized for sample size
	n = nrow(xyz)
	return(eigen(M / n))    }

####################################################################################
# Function written to identify swapovers as defined by Ringrose and Benn 1997, 
# Benn and Ringrose 2001, and Ringrose 1996.
#
ringrose_swap = function(e, s) {
  
  # Setup the six possible comparisons  
  combinations = matrix(NA, nrow = 6 , ncol = 3)
  combinations[1,] = c(1, 2, 3)
  combinations[2,] = c(2, 1, 3)
  combinations[3,] = c(3, 2, 1)
  combinations[4,] = c(2, 3, 1)
  combinations[5,] = c(3, 1, 2)
  combinations[6,] = c(1, 3, 2)

  # Compute the similarity of the six possible combinations of population and sample vectors
  scores = vector(length = 6)
  for (i in 1:6) {
    a = combinations[i,1]
    b = combinations[i,2]
    c = combinations[i,3]
    scores[i] = abs(e$vectors[,1] %*% s$vectors[,a]) + abs(e$vectors[,2] %*% s$vectors[,b]) + abs(e$vectors[,3] %*% s$vectors[,c])
  }

  # Any value other than 1, 2, 3 means a swap is warrented
  return(combinations[which(scores == max(scores)),]) 

}


#################################################################################################################
#  Return vector(artifact) lengths
#
compute_lengths = function(xyz) return( sqrt((xyz$X1-xyz$X2)^2 + (xyz$Y1-xyz$Y2)^2 + (xyz$Z1-xyz$Z2)^2) )

#################################################################################################################
# Ternary plot taken from VCD and heavily modify to create Benn Diagrams
#	Removed used of Grid functions
#	Simplied some arguements and added additional ones
#	X can either be eigenvalues or elongation/isotrophy values
#	Can do a special version with a color coded background useful in conjunction with color coded plots of artifact orientations
#
benn_diagram = function(x,  plot_points = TRUE, main = "Benn Diagram", dim_names = c("Planar","Linear","Isotropic"), 
                        dimnames_position = c("corner", "edge", "none"), dimnames_color = "black",
                        coordinates = FALSE, id = NULL, id_color = "black", id_just = c("center", "center"),
                        drawhull = FALSE, colorcoded = FALSE,
                        legend_names = NULL, legend_colors = NULL, 
                        draw_grid = TRUE, grid_color = "gray", grid_labels = TRUE,
                        labels = "outside", labels_color = "darkgray",
                        new_page = TRUE, border = "black", bg = "white", pch = 19, cex = 1,
                        prop_size = FALSE, col = "black", ...) 
{
	
  library(grid)
  
  scale = 1
	if (!labels %in% c("inside", "outside", "none")) stop("Labels must be 'inside','outside', or 'none'")
	if (coordinates)  id <- paste("(", round(x[, 1], 2), ",", round(x[, 2], 2), ",", round(x[, 3], 2), ")", sep = "")
	dimnames_position <- match.arg(dimnames_position)
	if (is.null(dim_names) && dimnames_position != "none")  dim_names <- colnames(x)
	if (is.logical(prop_size) && prop_size)  prop_size <- 3
	if (any(x < 0, na.rm = TRUE))  stop("EL and IS values must be between 0 and 1.")
	
	if ("IS"  %in% colnames(x) & ("EL"  %in% colnames(x))) 	x = cbind(x[,"EL"],x[,"IS"])
	
	s <- rowSums(x)
	if (any(s <= 0, na.rm = TRUE)) x=x[which(s!=0),]

	top <- .866						# This makes the left and right side a unit length

	xlim <- c(-0.04, 1.04)			# Set up plot limits that give a litte extra room for labels
	ylim <- c(-.04, 1.04)

	eps <- 0.01

	temp_margins = par(c("mai"))    	# Since we don't need margins for axis labels etc., save the current margins and
	par(mai = c(.01,.01,.01,.01))		  # push the margins to the limit
	
	# Set up an empty plot frame with a 1:1 aspect ratio
	if (new_page) {
  	plot(x = 0, y = 0, xlim = xlim, ylim = ylim, main = "",
  	     type = "n", asp = 1, xlab = "", ylab = "",
  	     frame.plot = FALSE, axes = FALSE) 
	
  	# If there is a title, put it at the top of the triangle
  	if (!is.null(main)) {
  	  if (dimnames_position == "corner") { 
  	          text(main, x = .5, y = top + .15,  col = dimnames_color, cex = cex)
  	  } else {
  	          text(main, x = .5, y = top + .1,  col = dimnames_color, cex = cex) } }
  
  	# Create the triangle itself with an optional fill
  	polygon(c(0, 0.5, 1), c(0, top, 0), col = bg, border = border,...)
      
  	# If there are corner names, and they are to be positioned in the corners - do it
  	if (dimnames_position == "corner") text(x = c(0, 1, 0.5), y = c(-0.04, -0.04, top + 0.04), label = dim_names, col = dimnames_color, cex = cex)
  
  	# If there are corner names, and they are to be on the edges - do it
  	if (dimnames_position == "edge") {
  		shift <- eps * if (labels == "outside") 8  else 0
  		text(x = 0.25 - 2 * eps - shift, y = 0.5 * top + shift, label = dim_names[2], rot = 60, col = dimnames_color, cex = cex)
  		text(x = 0.75 + 3 * eps + shift, y = 0.5 * top + shift, label = dim_names[1], rot = -60, col = dimnames_color, cex = cex)
  		text(x = 0.5, y = -0.02 - shift, label = dim_names[3], col = dimnames_color, cex = cex)
  	}
  	
  	# If a grid is asked for, draw it
  	if (draw_grid) {
  		for (i in 1:4 * 0.2) {
  			if (!colorcoded) lines(c(1 - i, 1 - i + i/2), c(0, i) * top, lty = "dotted", col = grid_color)
  			if (!colorcoded) lines(c(i/2, 1 - i + i/2), c(i, i) * top, lty = "dotted", col = grid_color)
  			if (labels == "inside") {
  				text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 * top, label = i * scale, col = labels_color, srt = 120, cex = cex * .8)
  				text(x = 1 - i + i/4 + eps, y = i/2 * top - eps, label = (1 - i) * scale, col = labels_color, srt = -120, cex = cex * .8)
  			}
  			if (labels == "outside") {
  				text(x = (1 - i)/2 - 6 * eps, y = (1 - i) * top, label = (1 - i) * scale, col = labels_color, cex = cex * .8)
  				text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 -  i) * top + 5 * eps, label = i * scale, srt = 60, col = labels_color, cex = cex * .8)
  			}
  		}
  		if (grid_labels) {
  		  text(x=.88,y=.5,label="Elongation Ratio 1-(E2/E1)",srt=-60,col = labels_color, cex = cex * .8 )
  		  text(x=.12,y=.5,label="Isotropy Ratio E3/E1",srt=60,col = labels_color, cex = cex * .8 ) } }
      
  	size = unit(if (prop_size)  prop_size * (s/max(s))  else cex, "lines")
  
  	# This routine creates the color coded Benn Diagram - i.e. a colored coded background and then black dots for the actual values
  	if (colorcoded) {
  	
  		# Make a color key for the Benn Diagram by first making a set of points to cover the diagram
  		d=matrix(nrow=sum(0:101),ncol=2)
  		k=0
  		for (x1 in 0:100) {
  			for (y1 in 0:(100-x1)) {
  				k=k+1
  				d[k,1]=x1/100
  				d[k,2]=y1/100 } }
  
  		b = benn_coords(cbind(elongation = d[,1], isotropy = d[,2]))
  		xp = b[,1]
  		yp = b[,2]
  	
  		# Now color code these points
  		
  		# One way to make colors
  		# point_color=rgb(1-elong,1-isotroph,0)
  		
  		# Another way to make colors
  		color_range=two.colors(n = 21, start = "green", end = "red", middle = "yellow")
  		point_color=color_range[xp * 20 + 1]
  		rgb_color=hex2RGB(point_color)
  		hsv_color=t(rgb2hsv(rgb_color@coords[,1], rgb_color@coords[,2], rgb_color@coords[,3]))
  		point_color=hsv(hsv_color[,1], 1 - yp, 1)
  
  		# And plot them - here the scaling is tricky to achieve the right look - adjust as necessary
  		points(xp, yp, pch = pch, col = point_color, cex =  cex * (par('fin')[1] / 6.92) * 1.1)
  
  	}
  }
	# Plot the actual points
	# If ncol=3 then eigen values were passed to this routine
	# If ncol=2 then elongation and isotrophy values were passed to this routine
	if (ncol(x)==3) {
		b = benn_coords(cbind(elongation = (x[,2]/x[,1]), isotropy = x[,3]/x[,1] ))
	} else {
		b = benn_coords(cbind(elongation = x[,1], isotropy = x[,2])) }

	xp = b[,1]
	yp = b[,2]

	if (plot_points) {
  	if (!colorcoded) {
  	  points(xp, yp, pch = pch, col = col, cex = cex * .8)
  	} else {
  	  points(xp, yp, pch = pch, col = col, cex = cex * (par('fin')[1] / 6.92) * .8) }
	} 
	
	# If specified, put a convex polygon around each group of points or around the whole thing.
	if (drawhull == TRUE) {
		if (length(table(col)) > 1) {
			for (colname in unique(col)) {
				hull=chull(xp[col==colname], yp[col==colname])
				if (length(hull)>2) { polygon(xp[col==colname][hull], yp[col==colname][hull], col = colname, border = NA) } } }
		else {
			hull = chull(xp, yp)
			polygon(xp[hull], yp[hull])
		}
	}
 
	# If labeled points are asked for, plot them
	if (!is.null(id))  text(x = xp, y = yp, label = as.character(id), col = id_color, cex = cex * .8, pos = 3)
			
	# Put a legend in the top right if a legend has been provided
	# Note that the points are colored on unique col but should be using l[2] - but this doesn't work properly
	#if (!is.null(l)) {
	#		text(x=.94, y=(1:length(l))*4/100+.65,
	#		     label=as.character(l),
	#		     col = "black", cex = cex)
	#		points(rep(.92,length(l)),(1:length(l))*4/100+.65,
	#		       pch = pch, col = unique(col), cex = size) }
	
	if (!is.null(legend_names)) legend(x = .68, y = 1.05, 
	                        legend = legend_names,
	                        fill = legend_colors, bty = 'n',
	                        cex = cex)
	
	par(mai=temp_margins)
}

#################################################################################################################
#  Basic bearing and plunge analysis
#
benn_orientations = function(xyz, level = "All Points", min_sample = 50, resampling = 1000,
                             permutations = FALSE, p = .95, ...) {
	
	level = factor(level)
	
	if (resampling>0) resampling_contours = vector("list", length(levels(level)))  

	benn = round(benn(xyz, level, min_sample), 3)
	
	# Prepare for two sets of three figures per page
	m=matrix(c(1,3,5,7,2,4,6,8), nrow = 4, ncol = 2) 	
	layout(m)
	
	for (l in levels(level)) {

		xyz_level = subset(xyz, level==l)
		
		if (nrow(xyz_level) >= min_sample) {
			
			benn_diagram(cbind(benn[l,"EL"],benn[l,"IS"]), id = l,
			             main = "Benn Diagram", cex = .7, labels = "outside")

			if (resampling>0) {
				resampling_contours[[which(levels(level)==l)]] = benn_resampling(xyz_level, resampling = resampling, p = p) 
				for (a_contour in resampling_contours[[which(levels(level)==l)]]) lines(benn_coords(cbind(elongation = a_contour$x, isotropy = a_contour$y))) } } }

	# Prepare for two sets of three figures per page
	m=matrix(c(1), nrow=1, ncol=1) 	
	layout(m)
				
	if (length(levels(level)) > 1) benn_diagram(cbind(benn[,"EL"], benn[,"IS"]), id = rownames(benn),
	                                            main = "Benn Diagram", cex = .7, labels = "outside")
	
	if (resampling > 0) {
		benn_diagram(cbind(benn[,"EL"], benn[,"IS"]), id = rownames(benn),
		             main = "Benn Diagram", cex = .7, labels = "outside")
		for (a_level in resampling_contours) {
			for (a_contour in a_level) {
					lines(benn_coords(cbind(elongation = a_contour$x, isotropy = a_contour$y))) } } }

	if (permutations & length(levels(level))>1) {
		benn_diagram(cbind(benn[,"EL"], benn[,"IS"]), id = rownames(benn),  main = "Benn Diagram", cex = .7, labels = "outside")
		perm_results = matrix(NA, nrow = length(levels(level)), ncol = length(levels(level)), 
		                      dimnames = list(levels(level), levels(level)))
		for (k in 1:(length(levels(level))-1)) {
			for (l in (k+1):length(levels(level))) {
				perm_results[levels(level)[k],levels(level)[l]] = benn_permutations(subset(xyz, level==levels(level)[k]),subset(xyz, level==levels(level)[l]))
				if (perm_results[levels(level)[k],levels(level)[l]] < p) {
					x0y0 = benn_coords(cbind(elongation=benn[levels(level)[k],"EL"],
					                         isotropy = benn[levels(level)[k],"IS"]))
					x1y1 = benn_coords(cbind(elongation=benn[levels(level)[l],"EL"],
					                         isotropy = benn[levels(level)[l],"IS"]))
					segments(x0y0[1],x0y0[2],x1y1[1],x1y1[2] ) } } } 
		return(list(benn,perm_results)) }
	
	return(benn)
}
