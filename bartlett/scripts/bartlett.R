##### Partitioning spectral diversity
# NEON Bartlett Forest site

### Open libraries ----
library(googledrive)
library(raster)
library(tidyverse)
library(hsdar)
library(gridExtra)
library(RStoolbox)
library(rasterVis)
library(viridis)
library(cowplot)
library(magick)



library(rgdal)
library(SphericalCubature)
library(signal)

### Read spectral diversity functions ----
source('functions/specdiv_functions.R')

### Download data ----

# Create folder to store data
if (!dir.exists('bartlett/data')) dir.create('bartlett/data')

# Hyperspectral data cube (Warning: 344 MB)
hyper_cube <- drive_download(file = as_id('155FcKLsrkrDSG24DVrzdkrQJXUFuVI-2'),
                             path = 'bartlett/data/hyper_cube.tif')

# Wavelenghts and full-width half maxima
wavelengths <- drive_download(file = as_id('18AXLw4KVdQwYtTGWW_aFLOfigUZU7V7V'),
               path = 'bartlett/data/wavelengths.csv')
fwhms <- drive_download(file = as_id('1bpwrfu7PYuddTO1AQBV7Qp8e0Hf0crkF'),
                              path = 'bartlett/data/fwhms.csv')

# NDVI
ndvi_tif <- drive_download(file = as_id('1Eumw5N7c6-uyFykNfr11Ho_J0b5DRed6'),
                           path = 'bartlett/data/ndvi.tif')

# Digital surface model (DSM)
dsm_tif <- drive_download(file = as_id('1EHanzWnvYWdOFyPptK8Jwxv3TccKd27u'),
                            path = 'bartlett/data/dsm.tif')

# High resolution RGB orthomosaic (Warning: 126 MB)
rgb_tif <- drive_download(file = as_id('1mKcq4B4ujmZNZ7CysU3EZ-kUPyI8b3SW'),
                          path = 'bartlett/data/rgb.tif')


#### Read data into R -----

# Wavelength data
wvl <- read_csv('bartlett/data/wavelengths.csv')
fwhm <- read_csv('bartlett/data/fwhms.csv')

# Hyperspectral data cube
cube <- brick('bartlett/data/hyper_cube.tif')
names(cube) <- wvl$band_name

# Show red, green and blue band wavelengths
(rgb_bands <- dplyr::filter(wvl, band %in% c(56, 28, 14)) )

# Plot RGB composite, with bands R = 56, G = 28, B = 14
png('bartlett/figures/cube_RGB.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube,
        r = 56, g = 28, b = 14, stretch = 'lin')
dev.off()

# Plot false colour composite, R = NIR (800 nm)
dplyr::filter(wvl, band %in% c(84, 56, 28))
png('bartlett/figures/cube_FCC.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube,
        r = 84, g = 56, b = 28, stretch = 'Lin')
dev.off()

# NDVI
ndvi <- raster('bartlett/data/ndvi.tif')
plot(ndvi)

# DSM
dsm <- raster('bartlett/data/dsm.tif')
plot(dsm)

# High resolution RGB
rgb <- brick('bartlett/data/rgb.tif')
png('bartlett/figures/hires_rgb.png', width = 2.8, height = 10, res = 1000, units = "in")
plotRGB(rgb)
dev.off()


### Pre-process spectra ----
# Remove bad bands, apply Savitzky-Golay filter

# Read cube as speclib
cube_speclib <- speclib(spectra = cube, wavelength = wvl$wvl, fwhm = fwhm$fwhm)

# Mask bad bands
band_mask <- data.frame(lb = c(383, 1340, 1790, 2400),
                        up = c(400, 1445, 1955, 2512))
hsdar::mask(cube_speclib) <- band_mask

# Apply Savitzky-Golay filter (Warning: takes a few minutes)
cube_sg <- noiseFiltering(cube_speclib, method = 'sgolay', n = 7)

# Transform to RasterBrick
cube_sg_brick <- brick(cube_sg)

# Identify bad bands that were removed
wvls <- wvl$wvl
window1 <- c(1340, 1445)
window2 <- c(1790, 1955)
window1_wvls <- which(wvls >= window1[1] & wvls <= window1[2])
window2_wvls <- which(wvls >= window2[1] & wvls <= window2[2])
below_400 <- which(wvls < 400)
above_2400 <- which(wvls > 2400)

# Get bad bands together
bad_bands <- c(below_400, window1_wvls, window2_wvls, above_2400)

wvl_bands <- wvl %>%
  dplyr::filter(!band %in% bad_bands)

# Add band names to cube_sg_brick
names(cube_sg_brick) <- wvl_bands$band_name


#### Apply NDVI mask -----
# Show NDVI < 0.8
ndvi_0.8 <- ndvi >= 0.8
plot(ndvi_0.8)

# Mask cube
cube_masked <- raster::mask(cube_sg_brick, ndvi_0.8, maskvalue = FALSE)

# Show red, green and blue band wavelengths
wvl_bands[c(10, 24, 52), ]

# Plot RGB composite, with bands R = 52, G = 24, B = 10
# Note: band indices have changed since some bands were removed
png('bartlett/figures/cube_RGB_masked.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube_masked, colNA = 'red',
        r = 52, g = 24, b = 10, stretch = 'lin')
dev.off()



#### Brightness normalization -----

# Apply brightness normalization
cube_norm <- calc(cube_masked, fun = bright_norm)
names(cube_norm) <- names(cube_masked)

# Plot RGB
png('bartlett/figures/cube_RGB_norm.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube_norm,
        r = 52, g = 24, b = 10, stretch = 'Lin')
dev.off()

# Plot false colour composite, R = NIR (800)
png('bartlett/figures/cube_FCC_norm.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube_norm,
        r = 80, g = 52, b = 24, stretch = 'Lin')
dev.off()

##### PCA-----
# Run PCA, scaling = 1 (preserves Euclidean distance) (Warning: takes several minutes)
pc_scaling1 <- pca(cube_norm, scaling = 1)

# Plot first three PCs as RGB composite
cube_pc_scaling1 <- pc_scaling1$cube_pc
png('bartlett/figures/cube_PCA_scaling1_RGB.png', width = 2.8, height = 10, res = 100, units = "in")
plotRGB(cube_pc_scaling1,
        r = 1, g = 2, b = 3, stretch = 'Lin')
dev.off()

# Plot PCs
cube_pc_plots_scaling1 <- build_pc_plot(cube_pc_scaling1)
png('bartlett/figures/cube_pc_plots_scaling1.png', width = 14, height = 10, res = 300, units = 'in')
show_pc_plot(cube_pc_plots_scaling1, nrow = 2)
dev.off()

# Get image cube with the 5 PCs
cube_pc_5_scaling1 <- dropLayer(cube_pc_scaling1, 6:17)

### Find number of pixels not masked or missing per plot (factor for plot = 40 x 40 pixels)
cube_npixels_scaling1 <- count_pixels(cube_pc_5_scaling1)
png('bartlett/figures/cube_npixels_scaling1.png', width = 4, height = 10, res = 100, units = 'in')
plot(cube_npixels_scaling1, 'prop') # proportion of pixels with non-missing values
dev.off()
plot(cube_npixels_scaling1, 'n') #number of pixels with non-missing values


### Calculate and partition spectral diversity using the 5 PCs
# Repeat 30 times (Warning: takes a few minutes)
cube_specdiv_scaling1 <- specdiv(cube_pc_5_scaling1, fact = 40, prop = 0.5, n = 30) # fact = expansion factor for plots (= 40 x 40 pixels)
sdiv_alpha_scaling1 <- cube_specdiv_scaling1$rasters$alpha_sdiv
lcsd_beta_scaling1 <- cube_specdiv_scaling1$rasters$beta_lcsd
fcsd_alpha_scaling1 <- cube_specdiv_scaling1$rasters$alpha_fcsd

## Plot the three PCs that contribute most to alpha and beta diversity
plotRGB(cube_pc_5_scaling1, r = 1, g = 2, b = 3, stretch = 'lin')
beta_PC_scaling1 <- ggRGB(cube_pc_5_scaling1, r = 1, g = 2, b = 3, stretch = 'lin') +
  theme_void() +
  ggtitle('PC 1,2,3') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10) )
beta_PC_scaling1

alpha_plot_scaling1 <- gplot(sdiv_alpha_scaling1) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(SD[alpha])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
alpha_plot_scaling1

beta_plot_scaling1 <- gplot(lcsd_beta_scaling1) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(LCSD[beta])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
beta_plot_scaling1

alpha_fcsd_plot_scaling1 <- gplot(fcsd_alpha_scaling1) +
  geom_raster(aes(fill = value)) + 
  facet_wrap(~variable, nrow = 1) +
  theme_void() + 
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_fill_viridis(direction = -1, option = 'magma', name = expression(FCSD[alpha]))
alpha_fcsd_plot_scaling1
ggsave('bartlett/figures/alpha_fcsd_plot_scaling1.png', 
       plot = alpha_fcsd_plot_scaling1, width = 5, height = 4, dpi = 600)

# Hi-res RGB plot (Warning: takes a few minutes)
hisres_rgb <- image_read('bartlett/figures/hires_rgb.png')
gplotRGB <- image_ggplot(hisres_rgb) +
  ggtitle('RGB') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10) )
gplotRGB 


# Make multi-panel plot with RGB, 3PCs, LCSD, SDalpha
# Warning: takes a few minutes
png('bartlett/figures/RGB_3PCs_LCSD_SD_scaling1.png', width = 4, height = 4, res = 1200, units = 'in')
plot_grid(gplotRGB, beta_PC_scaling1, beta_plot_scaling1, alpha_plot_scaling1, align = 'hv', nrow = 1)
dev.off()



### Try with shadows masked -----
# Create shadow layer, using the solar elevation/azimuth data from the different spectrometer flightlines

# Create slope and aspect from DSM
dsm_slope <- terrain(dsm, opt = 'slope')
dsm_aspect <- terrain(dsm, opt = 'aspect')

# Create shadow layer, using the solar elevation/azimuth data from the different spectrometer flightlines (obtained from cube metadata)
dsm_shade <- hillShade(dsm_slope, dsm_aspect,
                       angle = mean(c(34.41875, 33.70717, 33.091766, 32.477367, 31.954811) ),
                       direction = mean(c(143.19498, 145.84555, 148.40707, 151.21599, 153.92741)))
plot(dsm_shade)

# Create shade mask
shade_mask <- dsm_shade >= 0.4 
png('bartlett/figures/shade_mask.png', width = 4, height = 10, res = 100, units = "in")
plot(shade_mask, col=grey(c(0,1)), axes = F)
dev.off()

# Mask RGB
cube_no_shade <- raster::mask(cube, shade_mask, maskvalue = 0)
plotRGB(cube,
        r = 56, g = 28, b = 14, stretch = 'lin')
plotRGB(cube_no_shade,
        r = 56, g = 28, b = 14, stretch = 'lin')

# Save the plots
cube_with_shade_plot <- ggRGB(cube,
                         r = 56, g = 28, b = 14,
                         stretch = 'lin') +
  ggtitle('With shadows') +
  theme(plot.title = element_text(hjust = 0.5, size = 8) )
cube_with_shade_plot
cube_no_shade_plot <- ggRGB(cube_no_shade,
                         r = 56, g = 28, b = 14,
                         stretch = 'lin') +
  ggtitle('Shadows removed') +
  theme(plot.title = element_text(hjust = 0.5, size = 8) )
cube_no_shade_plot

# Make multi-panel plot shade and no shade
png('bartlett/figures/with_without_shade.png', width = 2, height = 4, res = 600, units = 'in')
plot_grid(cube_with_shade_plot, cube_no_shade_plot, align = 'hv', nrow = 1)
dev.off()


# Mask with shade
cube_norm_no_shade <- raster::mask(cube_norm, shade_mask, maskvalue = 0)

### Run analyses from no-shade cube
# Run PCA, scaling = 1 (preserves Euclidean distance) (Warning: takes several minutes)
pc_scaling1_no_shade <- pca(cube_norm_no_shade, scaling = 1)


# Plot first three PCs as RGB composite
cube_pc_scaling1_no_shade <- pc_scaling1_no_shade$cube_pc
plotRGB(cube_pc_scaling1_no_shade,
        r = 1, g = 2, b = 3, stretch = 'Lin')

# Plot PCs
cube_pc_scaling1_no_shade_plots <- build_pc_plot(cube_pc_scaling1_no_shade)
png('bartlett/figures/cube_pc_plots_no_shade_scaling1.png', width = 14, height = 10, res = 300, units = 'in')
show_pc_plot(cube_pc_scaling1_no_shade_plots, nrow = 2)
dev.off()


# Get image cube with the 5 PCs
cube_pc_5_scaling1_no_shade <- dropLayer(cube_pc_scaling1_no_shade, 6:17)

### Find number of pixels not masked or missing per plot (factor for plot = 40 x 40 pixels)
cube_npixels_scaling1_no_shade <- count_pixels(cube_pc_5_scaling1_no_shade)
plot(cube_npixels_scaling1_no_shade, 'prop') # proportion of pixels with non-missing values
setMinMax(cube_npixels_scaling1_no_shade)
plot(subset(cube_npixels_scaling1_no_shade, 'prop') > 0.35)

### Calculate and partition spectral diversity using the 5 PCs
cube_specdiv_scaling1_no_shade <- specdiv(cube_pc_5_scaling1_no_shade, fact = 40, prop = 0.35, n = 30) # fact = expansion factor for plots (= 40 x 40 pixels)
sdiv_alpha_scaling1_no_shade <- cube_specdiv_scaling1_no_shade$rasters$alpha_sdiv
lcsd_beta_scaling1_no_shade <- cube_specdiv_scaling1_no_shade$rasters$beta_lcsd
fcsd_alpha_scaling1_no_shade <- cube_specdiv_scaling1_no_shade$rasters$alpha_fcsd

# Plot the raster
plot(sdiv_alpha_scaling1_no_shade)
plot(sdiv_alpha_scaling1) # with shadows
plot(lcsd_beta_scaling1_no_shade)
plot(lcsd_beta_scaling1) # with shadows

# Look at relationship between sdiv alpha masked and unmasked
sdiv_scaling1_no_shade <- rasterToPoints(sdiv_alpha_scaling1_no_shade)
sdiv_scaling1 <- rasterToPoints(sdiv_alpha_scaling1)
plot(sdiv_scaling1_no_shade[, 'sdiv_alpha'], sdiv_scaling1[, 'sdiv_alpha'])
cor.test(sdiv_scaling1_no_shade[, 'sdiv_alpha'], sdiv_scaling1[, 'sdiv_alpha'])

sdiv_scaling1_both <- as_tibble(sdiv_scaling1_no_shade) %>% 
  left_join(as_tibble(sdiv_scaling1), by = c('x', 'y')) %>% 
  rename(sdiv_masked = sdiv_alpha.x, sdiv_shade = sdiv_alpha.y)

sdiv_scaling1_plot <- ggplot(sdiv_scaling1_both, aes(x = sdiv_shade, y = sdiv_masked)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab(expression(SD[alpha]~'(with shadows)')) +
  ylab(expression(SD[alpha]~'(shadows removed)')) +
  annotate("text", x = 0.002, y = 0.003, label = expression(italic(r)==0.93~'  '~italic(P)<'0.0001'))
sdiv_scaling1_plot 
ggsave('bartlett/figures/sdiv_shadows_vs_noshadows.png', plot = sdiv_scaling1_plot, width = 4, height = 4, dpi = 600)


# Look at relationship between LCSD beta masked and unmasked
lcsd_scaling1_no_shade <- rasterToPoints(lcsd_beta_scaling1_no_shade)
lcsd_scaling1 <- rasterToPoints(lcsd_beta_scaling1)
plot(lcsd_scaling1_no_shade[, 'lcsd_beta'], lcsd_scaling1[, 'lcsd_beta'])
cor.test(lcsd_scaling1_no_shade[, 'lcsd_beta'], lcsd_scaling1[, 'lcsd_beta'])

lcsd_scaling1_both <- as_tibble(lcsd_scaling1_no_shade) %>% 
  left_join(as_tibble(lcsd_scaling1), by = c('x', 'y')) %>% 
  rename(lcsd_masked = lcsd_beta.x, lcsd_shade = lcsd_beta.y)

lcsd_scaling1_plot <- ggplot(lcsd_scaling1_both, aes(x = lcsd_shade, y = lcsd_masked)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab(expression(LCSD[beta]~'(with shadows)')) +
  ylab(expression(LCSD[beta]~'(shadows removed)')) +
  annotate("text", x = 0.01, y = 0.025, label = expression(italic(r)==0.97~'  '~italic(P)<'0.0001'))
lcsd_scaling1_plot 
ggsave('bartlett/figures/lcsd_shadows_vs_noshadows.png', plot = lcsd_scaling1_plot, width = 4, height = 4, dpi = 600)


# Multi-panel figure showing the two relationships
sdiv_scaling1_plot2 <- sdiv_scaling1_plot +
  ggtitle('a') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0, size = 20, face = 'bold') )
lcsd_scaling1_plot2 <- lcsd_scaling1_plot +
  ggtitle('b') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0, size = 20, face = 'bold') )
png('bartlett/figures/sdiv_lcsd_with_without_shade.png', width = 7, height = 3.5, res = 300, units = 'in')
plot_grid(sdiv_scaling1_plot2, lcsd_scaling1_plot2, align = 'hv', nrow = 1)
dev.off()


## Plot the three PCs that contribute most to diversity
plotRGB(cube_pc_5_scaling1_no_shade, r = 1, g = 2, b = 3, stretch = 'Lin')
beta_PC_scaling1_no_shade <- ggRGB(cube_pc_5_scaling1_no_shade, r = 1, g = 2, b = 3, stretch = 'Lin') +
  theme_void() +
  ggtitle('PC 1,2,3') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10) )
beta_PC_scaling1_no_shade

alpha_plot_scaling1_no_shade <- gplot(sdiv_alpha_scaling1_no_shade) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(SD[alpha])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
alpha_plot_scaling1_no_shade

beta_plot_scaling1_no_shade <- gplot(lcsd_beta_scaling1_no_shade) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(LCSD[beta])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
beta_plot_scaling1_no_shade

# Make multi-panel plot with RGB, 3PCs, LCSD, SDalpha
png('bartlett/figures/RGB_3PCs_LCSD_SD_shadows_masked.png', width = 4, height = 4, res = 1200, units = 'in')
plot_grid(gplotRGB, beta_PC_scaling1_no_shade, beta_plot_scaling1_no_shade, alpha_plot_scaling1_no_shade, align = 'hv', nrow = 1)
dev.off()


##### Try with different plot sizes -----

# Define plot size ranging in width from 5 to 100 pixels (= 5 m to 1 ha wide)
plot_sizes <- seq(2, 140, by = 2)

# Number of resampling for each plot size
times <- 1

# For loop
sim_ss_table <- tibble()
lcsd_rasters <- list()
sdiv_rasters <- list()

for (i in 1:length(plot_sizes)) {
  # Run partitioning of spectral diversity
  specdiv_tmp <- specdiv(cube_pc_5_scaling1, fact = plot_sizes[i], prop = 0.5, n = times)
  
  # Store SS results
   ss_tmp <- specdiv_tmp$ss %>% 
    mutate(plot_size =  plot_sizes[i]) %>% 
    dplyr::filter(source != 'gamma')
  sim_ss_table <- bind_rows(sim_ss_table, ss_tmp)
  
  # Store LCSD-beta
  lcsd_rasters[[i]] <- specdiv_tmp$rasters$beta_lcsd
  
  # Store SD-alpha
  sdiv_rasters[[i]] <- specdiv_tmp$rasters$alpha_sdiv
}

#### Plot SS as a function of plot size
ss_size_plot <- ggplot(sim_ss_table, aes(x = plot_size, y = prop_gamma)) +
  geom_area(aes(fill = source)) +
  scale_fill_discrete(labels = c(expression(alpha), expression(beta)),
                    name = 'SS') +
  xlab('Plot size (number of pixels wide/high)') +
  ylab(expression('Proportion of'~SS[gamma])) +
  scale_x_continuous(breaks = seq(0, 140, by = 10)) +
  theme_minimal()
ss_size_plot
ggsave('bartlett/figures/prop_SS_gamma_plot_size.png', plot = ss_size_plot, width = 5, height = 3, units = 'in', dpi = 300)

### Plot LCSDs for selected plot sizes ------

# Add names corresponding to plot sizes
names(lcsd_rasters) <- plot_sizes

# Plot size = 2
lcsd_size2 <- gplot(lcsd_rasters$`2`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 2') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size2

# Plot size = 6
lcsd_size6 <- gplot(lcsd_rasters$`6`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 6') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size6


# Plot size = 10
lcsd_size10 <- gplot(lcsd_rasters$`10`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 10') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size10


# Plot size = 20
lcsd_size20 <- gplot(lcsd_rasters$`20`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 20') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size20

# Plot size = 30
lcsd_size30 <- gplot(lcsd_rasters$`30`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 30') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size30


# Plot size = 40
lcsd_size40 <- gplot(lcsd_rasters$`40`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 40') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size40


# Plot size = 60
lcsd_size60 <- gplot(lcsd_rasters$`60`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 60') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size60


# Plot size = 80
lcsd_size80 <- gplot(lcsd_rasters$`80`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 80') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size80


# Plot size = 100
lcsd_size100 <- gplot(lcsd_rasters$`100`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 100') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
lcsd_size100

### Multi-panel figure with LCSDs for different plot sizes
png('bartlett/figures/LCSDs_plot_sizes.png', width = 10, height = 4, res = 1200, units = 'in')
plot_grid(beta_PC_scaling1, lcsd_size2, lcsd_size6, lcsd_size10, lcsd_size20, lcsd_size30, lcsd_size40, lcsd_size60, lcsd_size100, align = 'hv', nrow = 1)
dev.off()



### Plot SD alpha for selected plot sizes ------

# Add names corresponding to plot sizes
names(sdiv_rasters) <- plot_sizes

# Plot size = 2
sdiv_size2 <- gplot(sdiv_rasters$`2`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 2') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size2

# Plot size = 6
sdiv_size6 <- gplot(sdiv_rasters$`6`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 6') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size6


# Plot size = 10
sdiv_size10 <- gplot(sdiv_rasters$`10`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 10') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size10


# Plot size = 20
sdiv_size20 <- gplot(sdiv_rasters$`20`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 20') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size20

# Plot size = 30
sdiv_size30 <- gplot(sdiv_rasters$`30`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 30') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size30


# Plot size = 40
sdiv_size40 <- gplot(sdiv_rasters$`40`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 40') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size40


# Plot size = 60
sdiv_size60 <- gplot(sdiv_rasters$`60`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 60') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size60


# Plot size = 80
sdiv_size80 <- gplot(sdiv_rasters$`80`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 80') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size80


# Plot size = 100
sdiv_size100 <- gplot(sdiv_rasters$`100`) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle('Size = 100') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
sdiv_size100

### Multi-panel figure with SD alpha for different plot sizes
png('bartlett/figures/SD_alpha_plot_sizes.png', width = 10, height = 4, res = 1200, units = 'in')
plot_grid(beta_PC_scaling1, sdiv_size2, sdiv_size6, sdiv_size10, sdiv_size20, sdiv_size30, sdiv_size40, sdiv_size60, sdiv_size100, align = 'hv', nrow = 1)
dev.off()


##### Try with PCA scaling = 2 -----

##### PCA-----
# Run PCA, scaling = 2 (preserves Mahalanobis distance) (Warning: takes several minutes)
pc_scaling2 <- pca(cube_norm, scaling = 2)

# Plot first three PCs as RGB composite
cube_pc_scaling2 <- pc_scaling2$cube_pc

# Plot PCs
cube_pc_plots_scaling2 <- build_pc_plot(cube_pc_scaling2)
show_pc_plot(cube_pc_plots_scaling2, nrow = 2)

# Get image cube with the 5 PCs
cube_pc_5_scaling2 <- dropLayer(cube_pc_scaling2, 6:17)

### Find number of pixels not masked or missing per plot (factor for plot = 40 x 40 pixels)
cube_npixels_scaling2 <- count_pixels(cube_pc_5_scaling2)
plot(cube_npixels_scaling2, 'prop') # proportion of pixels with non-missing values

### Calculate and partition spectral diversity using the 5 PCs
# Repeat 30 times (Warning: takes a few minutes)
cube_specdiv_scaling2 <- specdiv(cube_pc_5_scaling2, fact = 40, prop = 0.5, n = 30) # fact = expansion factor for plots (= 40 x 40 pixels)
sdiv_alpha_scaling2 <- cube_specdiv_scaling2$rasters$alpha_sdiv
lcsd_beta_scaling2 <- cube_specdiv_scaling2$rasters$beta_lcsd
fcsd_alpha_scaling2 <- cube_specdiv_scaling2$rasters$alpha_fcsd

## Plot the three PCs that contribute most to beta diversity
plotRGB(cube_pc_5_scaling2, r = 2, g = 3, b = 4, stretch = 'lin')
beta_PC_scaling2 <- ggRGB(cube_pc_5_scaling2, r = 2, g = 3, b = 4, stretch = 'lin') +
  theme_void() +
  ggtitle('PC 2,3,4') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10) )
beta_PC_scaling2

alpha_PC_scaling2 <- ggRGB(cube_pc_5_scaling2, r = 1, g = 3, b = 5, stretch = 'lin') +
  theme_void() +
  ggtitle('PC 1,3,5') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10) )
alpha_PC_scaling2

alpha_plot_scaling2 <- gplot(sdiv_alpha_scaling2) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(SD[alpha])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
alpha_plot_scaling2

beta_plot_scaling2 <- gplot(lcsd_beta_scaling2) +
  geom_raster(aes(fill = value)) + 
  theme_void() + 
  ggtitle(expression(LCSD[beta])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_viridis(direction = -1, option = 'magma')
beta_plot_scaling2

alpha_fcsd_plot_scaling2 <- gplot(fcsd_alpha_scaling2) +
  geom_raster(aes(fill = value)) + 
  facet_wrap(~variable, nrow = 1) +
  theme_void() + 
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_fill_viridis(direction = -1, option = 'magma', name = expression(FCSD[alpha]))
alpha_fcsd_plot_scaling2
ggsave('bartlett/figures/alpha_fcsd_plot_scaling2.png', 
       plot = alpha_fcsd_plot_scaling2, width = 5, height = 4, dpi = 600)


# Make multi-panel plot with RGB, 3PCs, LCSD, SDalpha
# Warning: takes a few minutes
png('bartlett/figures/RGB_3PCs_LCSD_SD_scaling2.png', width = 5, height = 4, res = 1200, units = 'in')
plot_grid(gplotRGB, beta_PC_scaling2, beta_plot_scaling2, alpha_PC_scaling2, alpha_plot_scaling2, align = 'hv', nrow = 1)
dev.off()

# Version without the RGB
png('bartlett/figures/3PCs_LCSD_SD_scaling2.png', width = 4, height = 4, res = 1200, units = 'in')
plot_grid(beta_PC_scaling2, beta_plot_scaling2, alpha_PC_scaling2, alpha_plot_scaling2, align = 'hv', nrow = 1)
dev.off()
