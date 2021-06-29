##########--------------------------------------#####>
##_FUN_##____Rededge image processing___________#####
##########--------------------------------------#####>
## AUTHOR      : Giacomo Crucil
## VERSION     : 20210504
## DESCRIPTION : Application of the radiometric calibration model and image correction: from raw Rededge images to calibrated reflectance
##               This script is to be used with the "Rededge-M" camera model with firmware from 2.10 version
## INPUT       : - "dir_CPR" = /text/ path to directory with CRP pictures in .tif format
##               - "dir_images" = /text/ path to directory with pictures to be converted to reflectance in .tif format
## OUTPUT      : /_r.tif/ processed images saved in subfolder "/reflectance"
## REQUIREMENTS : library(exifr)
##                library(raster)
##                library(rgdal)
##                library(spatstat)
##                configure_exiftool("YOUR-PATH-TO-exiftool.exe"/exiftool.exe")

RededgeM_processing <- function(dir_CRP, dir_images, metadata_trsf = T) {
  
#####@___Libraries___@#####
library(exifr)
library(raster)
library(rgdal)
library(spatstat)
######---------------#####>

#########------------------------------------------------#####>
##_subFUN_##____Rededge radiometric calibration model ___#####
#########------------------------------------------------#####>
## NAME        : rededge_corr
## AUTHOR      : Giacomo Crucil
## VERSION     : 20210504
## DESCRIPTION : RedEdge radiometric calibration converts the raw pixel values of an image into absolute spectral radiance values,
##               with units of W/m2/sr/nm. It compensates for sensor black-level, the sensitivity of the sensor,
##               sensor gain and exposure settings, and lens vignette effects. For info see:
##               https://support.micasense.com/hc/en-us/articles/115000351194-Rededge-Camera-Radiometric-Calibration-Model
##               This script is to be used with the "Rededge-M" camera model with firmware from 2.10 version
## INPUT       : - "input 1" = /text/ path to directory with pictures in .tif format
## OUTPUT      : /.tif/ corrected images saved in subfolder "/corr"

rededge_corr <- function(dir) { 
  image_list <- list.files(dir, pattern = ".tif", recursive = TRUE, full.names = TRUE)
  dir.create(paste0(dir, "/corr"))
  pb <- txtProgressBar(min = 0, max = length(image_list), style = 3)
    for (i in 1:length(image_list)) {
    metadata <- read_exif(image_list[i])
    image <- raster(image_list[i])
    ### collect parameters for image correction
    a <- as.numeric(metadata$RadiometricCalibration[[1]])
    Te <- metadata$ExposureTime
    g <- metadata$ISOSpeed/100
    pbl <- mean(as.numeric(strsplit(metadata$BlackLevel, " ")[[1]]))/65536
    cx <- metadata$VignettingCenter[[1]][1]
    cy <- metadata$VignettingCenter[[1]][2]
    pol <- as.numeric(metadata$VignettingPolynomial[[1]])
    ### calculate radiance image
    x_coord <- image
    x_coord[] <- xFromCell(image, 1:ncell(image))
    y_coord <- image
    y_coord[] <- yFromCell(image, 1:ncell(image))
    row_n <- image
    row_n[] <- rowFromCell(image, 1:ncell(image))
    image <- stack(image, x_coord, y_coord, row_n)
    image_L <- overlay(image, fun = function(image, x_coord, y_coord, row_n) {
      r <- sqrt((x_coord - cx)^2 + (y_coord - cy)^2)
      k <- 1 + pol[1]*r + pol[2]*r^2 + pol[3]*r^3 + pol[4]*r^4 + pol[5]*r^5 + pol[6]*r^6
      image_L <- 1/k * a[1]/g * (image/65536 - pbl)/(Te + a[2]*row_n - a[3]*Te*row_n)
    })
    ### assign (WGS84) crs to output images
    crs(image_L) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    ### write raster to disk in .tif
    writeRaster(image_L, filename = paste0(dir, "/corr/", substr(image_list[i], (nchar( image_list[1])-13), (nchar( image_list[1])-4)), "_corr.tif"),
                format = "GTiff")
    setTxtProgressBar(pb, i)
  }
}

#####@---------------@#############################################################################################>
##### Work-flow from image raw-DN to reflectance #####
#####@---------------@#############################################################################################>

### Rededge-M CRP serial: RP04-1808074-SC
CRP_Blue <- 0.4927
CRP_Green <- 0.4930
CRP_Red <- 0.4913
CRP_NIR <- 0.4874
CRP_Rededge <- 0.4901
CRP_reflectance <- c(CRP_Blue, CRP_Green, CRP_Red, CRP_NIR, CRP_Rededge)

#####___1) CRP data extraction #####
###################################################################################################################>

### Apply correction to CRP images
rededge_corr(dir_CRP)
### Select with mouse clicks (4 corners) the area of CRP for computing average Radiance L_avg
image_corr_list <- list.files(dir_CRP, pattern = "corr", recursive = TRUE, full.names = TRUE)
L_avg <- c()
for (i in 1:5) {
  plot(raster(image_corr_list[i]), main = paste0("click 4 corners inside each CRP: " , i, "/5" ))
  poly <- clickpoly(add = T, nv = 4, np = 1)
  ### convert clicked area to SpatialPolygon
  poly <- Polygon(matrix(unlist(poly$bdry), nrow = 4), hole = F)
  polys <- Polygons(list(poly), 1)
  sp_poly <- SpatialPolygons(list(polys))
  L_avg_i <- as.numeric(extract(raster(image_corr_list[i]), sp_poly, fun = mean, method='bilinear'))
  L_avg <- c(L_avg, L_avg_i)
}
print("CRP selection done. Now applying radiometric and vignetting correction...")
### Radiance_to_Reflectance "RtR" transfer function
RtR <- CRP_reflectance/L_avg

#####___2) Reflectance calculation #####
###################################################################################################################>

### Apply image correction function to all other images
image_list <- list.files(dir_images, pattern = ".tif", recursive = TRUE, full.names = TRUE)
rededge_corr(dir_images)
image_corr_list <- list.files(dir_images, pattern = "corr", recursive = TRUE, full.names = TRUE)
print("Radiometric and vignetting correction done. Now converting values to reflectance...")
### Apply RtR to all other images and save them in subfolder "/reflectance"
dir.create(paste0(dir_images, "/reflectance"))
dir_ref <- (paste0(dir_images, "/reflectance"))
pb <- txtProgressBar(min = 0, max = length(image_list), style = 3)
count <- 1
for (i in 1:5) { # i = n? of bands = n? of RtR factors
  for (j in seq(i,length(image_corr_list), 5)) { # j = sequential n? of the picture in the list, 5 = n? of bands
    image_R <- raster(image_corr_list[j])*RtR[i]
    writeRaster(image_R, filename = paste0(dir_images, "/reflectance/", substr(image_list[j], (nchar( image_list[1])-13), (nchar( image_list[1])-4)), "_r.tif"), format = "GTiff")
    setTxtProgressBar(pb, count)
    count <- count + 1
  }
}

### delete subfolder "/corr" with the images within (images in Radiance)
unlink(paste0(dir_CRP, "/corr"), recursive = T)
unlink(paste0(dir_images, "/corr"), recursive = T)

#####___3) Copy metadata #####
###################################################################################################################>

if(metadata_trsf == T) {
### Copy metadata from original to reflectance pictures
image_ref_list <- list.files(dir_ref, pattern = "_r.tif", full.names = TRUE)

print("Reflectance images saved on hard disk. Now transfering metadata...")
pb <- txtProgressBar(min = 0, max = length(image_list), style = 3)
for (i in 1:length(image_list)) {
  exiftool_call(args = c("-all=", "-TagsFromFile", image_list[i], "-all:all", "-xmp", image_ref_list[i], "-q", '-q'),
                quiet = T)
  setTxtProgressBar(pb, i)
}
file.remove(list.files(path = dir_ref, pattern = "_original", full.names = T)) # optional: remove original file copy
  print("Metadata transfered. Process complete")
} else {
  print("Reflectance images saved on hard disk. Process complete")
}

} # end of function




