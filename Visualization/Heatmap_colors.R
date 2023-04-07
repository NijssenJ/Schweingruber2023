
###### Non-scaled colors
col_nonscaled <- function(data, cellheight, aspect.ratio, by = 0.2, cent = (low+high) / 2, low, high) {
  cellheight = cellheight
  aspect.ratio = aspect.ratio
  cellwidth = ((nrow(data)*cellheight)/ncol(data))/aspect.ratio
  center <- (low+high) / 2
  if(length(cent) > 0){
    center <- cent
  }
  breaks <- c(seq(low, center, by = by), seq(center + by,high, by= by))
  color1 = colorRampPalette(c("Navy", "lightyellow1"))(length(seq(low-by, center, by = by)))
  color2 = colorRampPalette(c("lightyellow", "Firebrick4"))(length(seq(center+by,high+by, by = by)))
  color <- c(color1, color2)
  assign("color",color,.GlobalEnv)
  assign("breaks",breaks,.GlobalEnv)
  assign("cellheight",cellheight,.GlobalEnv)
  assign("cellwidth",cellwidth,.GlobalEnv)
}

# data: your final data file of which to plot the heatmap (typically log2 transformed RPKMs)
# cellheight: cellheight in pixels of each row in the heatmap.
#             Typically you have about 300 pixels fit nicely on one page, so if you have 100 rows, choose cellheight 3
# aspect.ratio: final aspect ratio of the plot: 2 means 2 high, 1 wide (long format). 0.5 means 2 wide, 1 high (wide format).
# by: interval of the color scale inside the heatmap (0.2 is a good default, only if your range becomes very small you might want to change it)


###### Scaled colors
col_scaled <- function(data, cellheight, aspect.ratio, cutoff = 1) {
  cellheight = cellheight
  aspect.ratio = aspect.ratio
  cellwidth= ((nrow(data)*cellheight)/ncol(data))/aspect.ratio
  low = "Navy"
  mid = "lightyellow1"
  high= "firebrick4"
  length = 50
  scaled <- t(data)
  scaled <- scale(scaled)
  edge = max(abs(min(scaled)), max(scaled))
  cut = cutoff * edge
  breaks = unique(c(-edge,  seq(-cut, 0, length=length),0, seq(0,cut, length=length), edge))
  
  col1 = colorRampPalette(c(low,mid))(length(seq(-cut, 0, length=length)))
  col2 = colorRampPalette(c(mid,high))(length(seq(0,cut, length=length)))
  color <- c(low, col1, col2, high)
  assign("cellheight",cellheight,.GlobalEnv)
  assign("cellwidth",cellwidth,.GlobalEnv)
  assign("aspect.ratio",aspect.ratio,.GlobalEnv)
  assign("color",color,.GlobalEnv)
  assign("breaks",breaks,.GlobalEnv)
}

# cutoff: in a scaled plot, the color range will span symmetrically around 0. If you want to bring the most intense colors
#         in a bit (from the edges, for better visualisation), you can play with this. A setting of 0.75 means that as of 
#         75% of the max value on either edge the most intense color is shown. Use values 0.5 - 1.