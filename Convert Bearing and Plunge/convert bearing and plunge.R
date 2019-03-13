# Example of how to convert bearing and plunge data set to 
# a set of vector suitable for analysis with the orientations.R
# library of fuctions.

# Note that for now the example data file is not included
# but it is a simple CSV file with three columns (Code, Dip, Orientation)

# For filtering/cleaning data
library(dplyr)

# From McPherron 2018
source('orientations.R') 

 # Read the data
data = read.csv('my site.csv',
                stringsAsFactors = FALSE)

# Little filtering and renaming of columns
data = data %>%
  filter(!is.na(Dip) & !is.na(Dip)) %>%
  select(code = Code, plunge = Dip, bearing = Orientation)

# A function to convert degrees to radians
rad <- function(deg) {(deg * pi) / (180)}

# Put all bearing angles to 0-360
data$bearing = data$bearing %% 360

# Convert plunge+bearing to direction vectors
# as if they came from a total station
X1 = rep(0, nrow(data))
Y1 = rep(0, nrow(data))
Z1 = rep(0, nrow(data))
X2 = sin(rad(data$bearing))
Y2 = cos(rad(data$bearing))
Z2 = -tan(rad(data$plunge))

# Organize in a data frame with names
# corresponding to McPherron 2018 functions
data_df = data.frame(code = data$code, X1, Y1, Z1, X2, Y2, Z2)

# Compute Benn values
benn_results = benn(data_df)

# Draw a Benn diagram with a 95% confidence interval
benn_diagram(benn_results, id = '',
             main = "My Site", cex = .7, labels = "outside")

# Add a 95% confidence interval based on 10000x resampling
resampling_contours = benn_resampling(data_df,
                                      resampling = 10000) 
for (a_contour in resampling_contours) {
  lines(benn_coords(cbind(elongation = a_contour$x,
                          isotropy = a_contour$y))) }

