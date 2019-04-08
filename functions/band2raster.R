### Read band to raster ----
#' @param path_h5 path to h5 file
#' @param band integer. spectral band to read
#' @param ncols integer. number of columns in the original image
#' @param nrows integer. number of rows in the original image
#' @param xstart integer. first pixel in x of output image, defaults to 1
#' @param ystart integer. first pixel in y of output image, defaults to 1
#' @param xstop integer. last pixel in x of output image, defaults ncols
#' @param ystop integer. last pixel in y of output image, defaults to nrows
#' @param crs CRS. coordinate reference system
#' @param res integer. pixel resolution in units of the CRS
#' @param noDataValue numeric. value in original image for no data
#' @return a one band raster object

band2raster <- function(path_h5, band, nrows, ncols, xstart=1, ystart=1, 
                        xstop=ncols, ystop=nrows, res, crs, xMin, yMax, NoDataValue){
  
  nam <- paste(h5ls(path_h5)$group, h5ls(path_h5)$name, sep="/")
  xMax <- (xMin + ncols*res) ### xMax = left edge + (no of cols*x pixel)
  yMin <- (yMax -  nrows*res)
  out <- h5read(path_h5,nam[grep("Reflectance_Data",nam)],
                index=list(band,xstart:xstop,ystart:ystop))
  
  #Convert from array to matrix
  out <- (out[1,,])
  out[out == NoDataValue] <- NA
  
  #turn the out object into a raster
  outr <- raster(t(out),crs=crs)
  
  # define the extents for the raster
  xMax <- xMin + (outr@ncols * res)
  yMin <- yMax - (outr@nrows * res)
  
  #create extents class
  rasExt  <- extent(xMin,xMax,yMin,yMax)
  
  #assign the extents to the raster
  extent(outr) <- rasExt
  
  #return the raster object
  return(outr)
}
