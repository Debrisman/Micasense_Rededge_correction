#####@___Libraries___@#####
library(doSNOW)
library(doParallel)
#
library(exifr)
library(raster)
library(suncalc)
######---------------#####>

### parallel comp setup
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, type="PSOCK")  
registerDoSNOW(cl) 

#####@---------------@#############################################################################################>
##### Pix4D calibrated parameters #####
#####@---------------@#############################################################################################>

#####!!! USER INPUT data from P4D /params folder
params_dir <- "/media/jack/BLACK_LOADER/20190420_Ernage2_Rededge/params" # params folder
#####!!! USER INPUT reflectance images folder
dir <- "/media/jack/BLACK_LOADER/20190420_Ernage2_Rededge/reflectance" # dir of reflectance images
#####!!! USER INPUT model calibrated from previous script
Input_model <- readRDS("/media/jack/BLACK_LOADER/Angles_model/PolyModels_ZenSenZen.rds")

# automatic input data selection
Cam_pos_file <- list.files(params_dir, pattern = "calibrated_external_camera_parameters_wgs84.txt", full.names = TRUE)
image_list <- list.files(dir, pattern = ".tif")
F_file <-  list.files(params_dir, pattern = "pix4d_calibrated_internal_camera_parameters.cam", full.names = TRUE)
F_par <- readLines(F_file)[c(3, 18, 33, 48, 63)] # extract calibrated Focal length

### Camera Focal distance F (mm) (calibrated from Pix4D)
F <- signif(as.numeric(substr(F_par, 3, 24)), 10)
Fp <- F*1280/4.8 # Focal lenght in pixels
### Cu Cv centers of images in image reference sys
ucvc <- matrix(c(1280/2, 960/2), ncol = 2, byrow = TRUE)
# rownames(ucvc) <- c("Cam1", "Cam2", "Cam3", "Cam4", "Cam5")
colnames(ucvc) <- c("uc", "vc")
### Camera positions and attitute (calibrated from Pix4D)
Cam_pos <- read.table(Cam_pos_file, header = T)
print(paste0(cat("!!!WARNING!!! Check if image name in 'params' is consistent with images in your 'dir' folder and correct accordingly\n"), 
             "In params folder = '", head(Cam_pos)[1, 1], "' // in dir folder =  '", image_list[1], "'"))

#####!!! USER INPUT imageName modification
# Cam_pos$imageName <- paste0(substr(Cam_pos$imageName, 1, 20), substr(Cam_pos$imageName, 26, 29)) # for Sicy
# Cam_pos$imageName <- paste0(substr(Cam_pos$imageName, 1, 10), "_r", substr(Cam_pos$imageName, 11, 14)) # for Thorembais
# Cam_pos$imageName <- paste0(substr(Cam_pos$imageName, 1, 18), "_r", substr(Cam_pos$imageName, 19, 22)) # Geeste_1, PetitSeumoy
# Cam_pos$imageName[1:3840] <- paste0("0000SET_", Cam_pos$imageName[1:3840])  # Ernage 2
# Cam_pos$imageName[3841:5635] <- paste0("0001SET_", Cam_pos$imageName[3841:5635])  # Ernage 2
Cam_pos$imageName <- paste0(substr(Cam_pos$imageName, 1, 18), "_r", substr(Cam_pos$imageName, 19, 22)) # Ernage 2
# Cam_pos$imageName <- paste0(substr(Cam_pos$imageName, 1, 10), "_r", substr(Cam_pos$imageName, 11, 14)) # for Gembloux 20210401/0423, Villeroux 20190827F3
#####!!!

Cam_pos <- SpatialPointsDataFrame(coords = coordinates(Cam_pos[, c("longitude", "latitude")]), data = Cam_pos,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
Cam_pos$Fp <- rep(Fp, nrow(Cam_pos)/5)
### Calc of nadir position in image reference sys
## u (x columns)
Cam_pos$ucuh <- Cam_pos$Fp*tan(Cam_pos$Phi*pi/180) # x pixel distance from image center (Uc) to sensor's vertical (Uh) 
Cam_pos$uc <- rep(ucvc[, 1], nrow(Cam_pos))
Cam_pos$uCal <- Cam_pos$uc + Cam_pos$ucuh # x pixel position of calibrated sensor's vertical (Ucal)
## v (y rows)
Cam_pos$vcvh <- Cam_pos$Fp*tan(Cam_pos$Omega*pi/180) # y pixel distance from image center Vc to sensor's vertical Vh 
Cam_pos$vc <- rep(ucvc[, 2], nrow(Cam_pos))
Cam_pos$vCal <- Cam_pos$vc + Cam_pos$vcvh # y pixel position of calibrated sensor's vertical (Vcal)

#####@---------------@#############################################################################################>
##### Image lists #####
#####@---------------@#############################################################################################>

### all images for post-processing
image_list <- list.files(dir, recursive = TRUE, pattern = ".tif")
image_list <- image_list[match(Cam_pos$imageName, image_list)]
image_list_full <- paste0(dir, "/", image_list)
head(image_list_full)[1]

#####@---------------@#############################################################################################>
##### Functions for angles #####
#####@---------------@#############################################################################################>

### fun_img_az to calc sensor azimut map, 0 to 2pi clockwise, respect to calibrated image coord syst
fun_img_az360 <- function(u_coord, v_coord, uCal_dist, vCal_dist) {
  ifelse(u_coord >= uCal & v_coord >= vCal, atan(uCal_dist/vCal_dist),
         ifelse(u_coord >= uCal & v_coord < vCal, atan(abs(vCal_dist)/uCal_dist) + pi/2, 
                ifelse(u_coord < uCal & v_coord < vCal, atan(abs(uCal_dist)/abs(vCal_dist)) + pi, 
                       atan(vCal_dist/abs(uCal_dist)) + 3*pi/2)))}
### fun_sen_az360 to calc sensor azimuth in a 360 degree system
fun_sen_az360 <- function(img_az360, kappa) {
  ifelse(img_az360 >= kappa, img_az360 - kappa, 2*pi + img_az360 - kappa)}
### fun_sen_az180 to calc sensor azimuth in a 180 degree system
fun_sen_az180 <- function(sen_az360) {
  ifelse(sen_az360 >= pi, 2*pi - sen_az360, sen_az360)}
### fun for alpha angle calculation (angle between calibrated sensor zenith position and pixels)
alpha <- function (u_coord, v_coord, uc_dist, vc_dist, uCal_dist, vCal_dist) {
  f <- 5 # test select F for cam
  TN  <- Fp[f]/cos(Cam_pos$Omega[img_num]*pi/180)
  CP <- sqrt(abs(uc_dist)^2 + abs(vc_dist)^2)
  TP <- sqrt(Fp[f]^2 + CP^2)
  NP <- sqrt(abs(uCal_dist)^2 + abs(vCal_dist)^2)
  cos.alpha <- (TP^2 + TN^2 - NP^2)/(2*TP*TN)
  alpha <- acos(cos.alpha)}
### fun_sen_az360_rel to calc relative azimuth in a 360 degreees system
fun_sen_az360_rel <- function(img_az360, sun_az, kappa) {
  ifelse(kappa < sun_az, 
         ifelse(img_az360 >= sun_az - kappa, img_az360 - sun_az + kappa, 2*pi + img_az360 - sun_az + kappa),
         ifelse(img_az360 >= 2*pi - (kappa - sun_az), img_az360 - (2*pi - (kappa - sun_az)), 2*pi + img_az360 - (2*pi - (kappa - sun_az))))}
### fun_sen_az180_rel to calc relative azimuth in a 180 degreees system
fun_sen_az180_rel <- function(sen_az360_rel) {
  ifelse(sen_az360_rel >= pi, 2*pi - sen_az360_rel, sen_az360_rel)}

#####@---------------@#############################################################################################>
##### Start of correction cycle #####
#####@---------------@#############################################################################################>

dir.create(paste0(dir, "/corr"))
pb <- txtProgressBar(max = length(image_list_full), style = 3)
progress <- function(img_num) setTxtProgressBar(pb, img_num)
foreach (img_num = 1:length(image_list_full), # start of sampling cycle
         .options.snow = list(progress = progress),
         .packages = c("raster", "exifr", "suncalc")) %dopar% {

#####@---------------@#############################################################################################>
##### Solar angles #####
#####@---------------@#############################################################################################>

### calc of sun position for sample image
image <- raster(image_list_full[img_num])
metadata <- read_exif(image_list_full[img_num])
rig <- metadata$RigCameraIndex
img_coord <- as.numeric(Cam_pos@data[img_num, c("longitude", "latitude")])
DateTime <- as.POSIXct(metadata$DateTimeOriginal, tz = "UTC", format = "%Y:%m:%d %H:%M:%OS")
# sun_pos <- sunAngle(t = DateTime, latitude = img_coord[2], longitude = img_coord[1]) # azimuth in degrees eastward of north, from 0 to 360
sun_pos <- getSunlightPosition(date = DateTime, lat = img_coord[2], lon = img_coord[1])
sun_pos$azimuth <- ifelse(sun_pos$azimuth >= 0, sun_pos$azimuth + pi, pi - abs(sun_pos$azimuth))

#####@---------------@#############################################################################################>
##### Sensor angles #####
#####@---------------@#############################################################################################>

### common basemaps
u_coord <- image
u_coord[] <- xFromCell(image, 1:ncell(image)) # raster of u(x) coord
v_coord <- image
v_coord[] <- yFromCell(image, 1:ncell(image)) # raster of v(y) coord

uCal <- floor(Cam_pos$uCal[img_num]) # u(x) of calibrated nadir center
uCal_dist <- u_coord - uCal # raster of u(x) pixel distance from calibrated nadir center 
uc <- floor(Cam_pos$uc)[img_num] # u(x) of image center
uc_dist <- u_coord - uc # raster of u(x) pixel distance from image center 

vCal <- floor(Cam_pos$vCal[img_num]) # v(y) of calibrated nadir center
vCal_dist <- v_coord - vCal # raster of v(y) pixel distance from calibrated nadir center 
vc <- floor(Cam_pos$vc)[img_num] # v(y) of image center
vc_dist <- v_coord - vc # raster of v(y) pixel distance from image center 

#####___Sensor azimuth #####
###################################################################################################################>

### apply fun_img_az
topo_stack <- stack(u_coord, v_coord, uCal_dist, vCal_dist)
img_az360 <- overlay(topo_stack, fun = fun_img_az360) # output calibrated image azimuth raster
### calc sensor azimuth respect to world coord syst
kappa <- ifelse(Cam_pos$Kappa[img_num] < 0, abs(Cam_pos$Kappa[img_num]*pi/180),
                2*pi - Cam_pos$Kappa[img_num]*pi/180)
kappa <- raster(ncol = 1280, nrow = 960, v = kappa, ext = extent(img_az360))
sen_az360 <- stack(img_az360, kappa)
sen_az360 <- overlay(sen_az360, fun = fun_sen_az360)
sen_az180 <- overlay(sen_az360, fun = fun_sen_az180)

#####___Sensor zenith #####
###################################################################################################################>

### apply "alpha" fun 
sen_zen <- stack(u_coord, v_coord, uc_dist, vc_dist, uCal_dist, vCal_dist)
sen_zen <- overlay(sen_zen, fun = alpha)
sen_zen[is.na(sen_zen[])] <- 0

#####@---------------@#############################################################################################>
##### Relative sensor-sun angles #####
#####@---------------@#############################################################################################>

#####___Relative sensor-sun azimuth #####
###################################################################################################################>

rev_sun_pos <- ifelse(sun_pos$azimuth <= pi, sun_pos$azimuth + pi, sun_pos$azimuth - pi) # calc opposite of solar angle to have 0 when opposing sun, pi when facing it 
sun_az <- raster(ncol = 1280, nrow = 960, v = rev_sun_pos, ext = extent(image))
sen_az360_rel <- stack(img_az360, sun_az, kappa)
sen_az360_rel <- overlay(sen_az360_rel, fun = fun_sen_az360_rel)
sen_az180_rel <- overlay(sen_az360_rel, fun = fun_sen_az180_rel)

#####___Relative sensor-sun zenith #####
#################### ###############################################################################################>

sun_zen <- raster(ncol = 1280, nrow = 960, v = pi/2 - sun_pos$altitude, ext = extent(image))
sen_zen_rel <- sen_zen * -cos(sen_az180_rel) + sun_zen # '-cos' accounts for custom solar angle (0 = opposing sun)

#####@---------------@#############################################################################################>
##### Correction application #####
#####@---------------@#############################################################################################>

### data sampling ANIF vs angles
Angles_map <- stack(sen_az180_rel, sen_zen_rel, sen_zen, sun_zen)
names(Angles_map) <- c("Rel_az", "Rel_zen", "Sensor_zen", "Sun_zen")
i_fit <- rig + 1
ANIF_pred <- raster::predict(Angles_map, Input_model[[i_fit]])
image_corr <- image/ANIF_pred
writeRaster(image_corr, filename = paste0(dir, "/corr/", gsub(pattern = "_r.tif", replacement = "_corr.tif", basename(image_list_full[img_num]))),
            format = "GTiff")

#####@---------------@#############################################################################################>
##### Stop correction cycle #####
#####@---------------@#############################################################################################>

} # end of corr cycle

#####@---------------@#############################################################################################>
##### Metadata transfer #####
#####@---------------@#############################################################################################>

image_corr_list <- paste0(dir, "/corr/", image_list)
image_corr_list <- gsub(pattern = "_r.tif", replacement = "_corr.tif", image_corr_list)

pb <- txtProgressBar(max = length(image_corr_list), style = 3)
progress <- function(i) setTxtProgressBar(pb, i)
foreach (i = 1:length(image_corr_list), .options.snow = list(progress = progress),
         .packages = "exifr") %dopar% {
  exiftool_call(args = c("-all=", "-TagsFromFile", image_list_full[i], "-all:all", "-xmp", image_corr_list[i], "-q", "-q"),
                quiet = T)
}
file.remove(list.files(path = paste0(dir, "/corr"), pattern = "_original", full.names = T)) # optional: remove original file copy
stopCluster(cl)
print("Metadata transfered. Process complete")

