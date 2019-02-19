# Haversine function for great circle distance ----------------------------
# Modified from https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
haversine <- function(lat1, lat2, long1, long2) {
  deg2rad <- function(deg) deg*pi/180
  
  R <- 6371 # Mean earth radius (km)
  
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)
  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  
  delta.lat <- lat1 - lat2
  delta.long <- long1 - long2
  
  a <- sin(delta.lat/2)^2 + cos(lat1)*cos(lat2)*sin(delta.long/2)^2
  c <- 2*asin(purrr::map_dbl(a, function(x) min(1, sqrt(x))))
  c*R
}
