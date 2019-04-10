# Figure to show alpha, beta, gamma

#Load libraries
library(raster)
library(tidyverse)
library(RStoolbox)
library(rasterVis)
library(cowplot)

# Wavelength data
wvl <- read_csv('bartlett/data/wavelengths.csv')
fwhm <- read_csv('bartlett/data/fwhms.csv')

# Hyperspectral data cube
cube <- brick('bartlett/data/hyper_cube.tif')
names(cube) <- wvl$band_name

# Get extent of cube
cube_extent <- extent(cube)

# Crop cube
cube_extent@xmax <- 314730
cube_extent@ymax <- 4880120
cube_cropped <- crop(cube, cube_extent)

# Plot it
plotRGB(cube_cropped,
        r = 56, g = 28, b = 14, stretch = 'lin')

# Plot as FCC
plotRGB(cube_cropped,
        r = 80, g = 52, b = 24, stretch = 'Lin')

# Plot with bands with R = NIR, G = red, B = lignin feature (2300 nm)
plotRGB(cube_cropped,
        r = 80, g = 52, b = 384, stretch = 'Lin')

# Get true RGB hi-resolution
rgb <- brick('bartlett/data/rgb.tif')

# Crop RGB
rgb_cropped <- crop(rgb, cube_extent)
plotRGB(rgb_cropped)


# Crop cube to retain only the three bands
cube_cropped2 <- subset(cube_cropped, c(80, 52, 384))
plotRGB(cube_cropped2, stretch = 'lin')

# Calculate sum squares
cube_cropped_ss <- rasterToPoints(cube_cropped2) %>% 
  as_tibble() %>% 
  mutate_at(c("Band_80", "Band_52", "Band_384"), function(x) (x - mean(x))^2)
cube_cropped_ss_xy <- dplyr::select(cube_cropped_ss, x, y)
cube_cropped_ss_values <- dplyr::select(cube_cropped_ss, -x, -y)
cube_cropped_ss_df <- SpatialPixelsDataFrame(cube_cropped_ss_xy, cube_cropped_ss_values, proj4string = CRS(proj4string(cube_cropped2) ))
cube_cropped3 <- brick(cube_cropped_ss_df)
plotRGB(cube_cropped3, stretch = 'hist')

# new cube with plots
fact = 20
new_res <- res(cube_cropped2) * fact
cube_plots <- raster(crs = proj4string(cube_cropped2))
extent(cube_plots) <- extent(cube_cropped2)
res(cube_plots) <- new_res
cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
cube_plots_df <- rasterToPoints(cube_plots) %>% 
  as_tibble()

# Plot it with grid
cropped_cube_plot <- ggRGB(cube_cropped2, r = 1, g = 2, b = 3, stretch = 'lin') +
  theme_void() +
  geom_tile(data = cube_plots_df, aes(x = x, y = y), fill = NA, colour = 'white', size = 1)
cropped_cube_plot

# Plot RGB with grid
rgb_tile <- ggRGB(rgb_cropped, r = 1, g = 2, b = 3) +
  theme_void() +
  geom_tile(data = cube_plots_df, aes(x = x, y = y), fill = NA, colour = 'white', size = 1)
rgb_tile

# plot only plot centroids
cube_pixels_df <- rasterToPoints(disaggregate(cube_plots, fact) ) %>% 
  as_tibble()
cropped_cube_df <- rasterToPoints(cube_cropped2) %>% 
  as_tibble() %>% 
  left_join(cube_pixels_df) %>% 
  rename(group = layer)

# Plot beta component
cropped_cube_beta <- cropped_cube_df %>% 
  group_by(group) %>% 
  summarise_all(mean)
cropped_cube_beta_xy <- dplyr::select(cropped_cube_beta, x, y)
cropped_cube_beta_values <- dplyr::select(cropped_cube_beta, -x, -y, -group)
beta_df <- SpatialPixelsDataFrame(cropped_cube_beta_xy, cropped_cube_beta_values, proj4string = CRS(proj4string(cube_cropped2) ))
beta_brick <- brick(beta_df)
plotRGB(beta_brick, stretch = 'lin')

# Plot alpha component
subs_mean <- function(x) {
  xy <- x %>% 
    dplyr::select(x, y)
  values <- x %>% 
    dplyr::select(-x, -y, -group) %>%
    mutate_all(function(y) y - mean(y)) %>% 
    bind_cols(xy)
  values$x <- xy$x
  values$y <- xy$y
  return(values)
}

cropped_cube_alpha <- cropped_cube_df %>% 
  group_by(group) %>% 
  do(subs_mean(.)) %>% 
  ungroup() %>% 
  dplyr::select(-group)
cropped_cube_alpha_xy <- dplyr::select(cropped_cube_alpha, x, y)
cropped_cube_alpha_values <- dplyr::select(cropped_cube_alpha, -x, -y)
alpha_df <- SpatialPixelsDataFrame(cropped_cube_alpha_xy, cropped_cube_alpha_values, proj4string = CRS(proj4string(cube_cropped2) ))
alpha_brick <- brick(alpha_df)
plotRGB(alpha_brick, stretch = 'lin')

# test cube_cropped2 specdiv
cube_cropped_specdiv <- specdiv(cube_cropped2)


# Make ggplots
rgb_tile2 <- ggRGB(rgb_cropped, r = 1, g = 2, b = 3) +
  theme_void() +
  geom_tile(data = cube_plots_df, aes(x = x, y = y), fill = NA, colour = 'white', size = 0.5) +
  ggtitle('Split region into communities') +
  theme(plot.title = element_text(hjust = 0.5, size = 13))
rgb_tile2

gamma_plot <- ggRGB(cube_cropped2, r = 1, g = 2, b = 3, stretch = 'lin') +
  theme_void() +
  #geom_tile(data = cube_plots_df, aes(x = x, y = y), fill = NA, colour = 'white', size = 0.5) +
  ggtitle(expression(Spectral~gamma*'-diversity')) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))
gamma_plot

beta_plot <- ggRGB(beta_brick, r = 1, g = 2, b = 3, stretch = 'lin') +
  theme_void() +
  ggtitle(expression(beta~component)) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))
beta_plot


alpha_plot <- ggRGB(alpha_brick, r = 1, g = 2, b = 3, stretch = 'lin') +
  theme_void() +
  geom_tile(data = cube_plots_df, aes(x = x, y = y), fill = NA, colour = 'white', size = 0.5) +
  ggtitle(expression(alpha~component)) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))
alpha_plot

# Make multi-panel figure
png('bartlett/figures/conceptual_figure_partitioning.png', width = 11, height = 2.5, res = 1200, units = 'in')
plot_grid(rgb_tile2, gamma_plot, beta_plot, alpha_plot, align = 'hv', nrow = 1)
dev.off()
