###### Partitioning spectral diversity
# Simulated example with leaf-level data

### Set working directory ----
setwd("/Users/etienne/Documents/ProjetsRecherche/CABO/specdiv")

### Source specdiv functions ----
source('functions/specdiv_functions.R')

### Load libraries ----
library(tidyverse)
library(raster)
library(sp)
library(vegan)
library(ggpubr)
library(viridis)
library(rasterVis)


### Load spectra data ----
refl_trees <- read_csv('simulated_example/data/refl_trees.csv')

# Check data
refl_trees

### Brightness normalization of the spectra ----
refl_trees <- refl_trees %>% 
  group_by(plantID) %>% 
  do(bright_norm_df(.)) %>% 
  ungroup()

# check data again
refl_trees


#### Principal Component Analysis (PCA) ----

# Convert data table to wide format for PCA
refl_wide <- refl_trees %>%
  dplyr::select(plantID, wvl_nm, refl_norm) %>% 
  spread(wvl_nm, refl_norm)

# Convert to matrix and move plantID as row names
refl_mat <- as.matrix(refl_wide[, names(refl_wide) != 'plantID'] )
rownames(refl_mat) <- refl_wide$plantID

# Perform PCA, type-1 and type-2 scaling
refl_pca <- pca_mat(refl_mat, scaling = 1, p = .99999) # type-1 scaling, retain PCs contributing to >99.999% total variation
refl_pca2 <- pca_mat(refl_mat, scaling = 2, p = .99999) # type-2 scaling, retain PCs contributing to >99.999% total variation


# Extract eigenvalues %
pc1 <- round(refl_pca$prop[1] * 100, 0) # percentage
pc2 <- round(refl_pca$prop[2] * 100, 0) # percentage
pc3 <- round(refl_pca$prop[3] * 100, 0) # percentage

# Extract the PCA site scores
site_scores <- bind_cols(plantID = refl_wide$plantID,
                          data.frame(refl_pca$obj) ) %>% 

# Create RGB colours using the first 3 PCs
  mutate(green = ((PC1 - min(PC1)) / (max(PC1) - min(PC1)) ),
       red = ((PC2 - min(PC2)) / (max(PC2) - min(PC2)) ),
       blue = ((PC3 - min(PC3)) / (max(PC3) - min(PC3)) ),
       RGB_col = rgb(red, green, blue) )

# Add the species names
tree_sp <- refl_trees %>%
  dplyr::select(plantID, scientificName) %>% 
  distinct()
site_scores <- site_scores %>% 
  left_join(tree_sp, by = 'plantID')

# Check data
site_scores

# Get mean colour for polygons
col_means <- site_scores %>% 
  group_by(scientificName) %>% 
  summarise(red = mean(red),
            green = mean(green),
            blue = mean(blue)) %>% 
  mutate(RGB_col = rgb(red, green, blue) )

# Draw PCA biplot, type 1 scaling
pca_plot <- ggplot(site_scores, aes(x = PC1, y = PC2)) +
  stat_chull(aes(group = scientificName, fill = scientificName), alpha = 0.5, geom = 'polygon') +
  geom_point(colour = site_scores$RGB_col, size = 2.5) +
  scale_fill_manual(breaks = col_means$scientificName, values = col_means$RGB_col, name = '') +
  xlab(bquote('PC 1 (' * .(pc1) *'%)')) +
  ylab(bquote('PC 2 (' * .(pc2) *'%)')) +
  coord_equal() +
 # xlim(c(-2.7, 3.7)) +
  theme_bw() +
  theme(legend.justification=c(.95,0.05), legend.position=c(.95,0.05),
        panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

# Convert x-axis label
pca_plot2 <- pca_plot + xlab('PC 1 (71.0%)')


# Save plot as PDF and PNG
ggsave('simulated_example/figures/pca_plot.png', plot = pca_plot, width = 4, height = 3, dpi = 600)
ggsave('simulated_example/figures/pca_plot.pdf', plot = pca_plot, width = 4, height = 3)


### Plot the band scores ----
band_scores <- bind_cols(wvl_nm = as.integer(as.character(rownames(refl_pca2$descript))),
                         data.frame(refl_pca2$descript))

# Convert to long format
band_long <- band_scores %>% 
  gather(key = PC, value = loading, -wvl_nm) %>% 
  dplyr::filter(PC %in% c('PC1', 'PC2', 'PC3')) %>% 
  mutate(PC = factor(PC, levels = c('PC1', 'PC2', 'PC3'),
                     labels = c('PC 1 (71.0%)',
                                'PC 2 (22.9%)',
                                'PC 3 (3.2%)')))

# Plot the PC loadings, with wavelength as colours
loadings_plot_colours <- ggplot(band_long, aes(x = wvl_nm, y = loading)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_bar(aes(fill = wvl_nm), stat = 'identity', width = 10) +
  facet_wrap(~ PC, ncol = 7) +
  ylab('PC score') +
  xlab('Wavelength (nm)') +
  scale_x_continuous(breaks = seq(400, 2400, 400)) +
  scale_y_continuous(breaks = seq(-0.006, .006, by = .002)) +
  scale_fill_viridis(option = 'magma', name = 'Wvl (nm)', direction = -1) +
  #scale_fill_gradientn(colours = rainbow(20),
  #                     name = 'Wvl (nm)') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
loadings_plot_colours
  
# Save plot as PDF and PNG
ggsave('simulated_example/figures/loadings.png', plot = loadings_plot_colours, width = 8, height = 3, dpi = 600)
ggsave('simulated_example/figures/loadings.pdf', plot = loadings_plot_colours, width = 8, height = 3) 

### PCA plot of loadings ---- 
PCA_loading_plot <- ggplot(band_scores, aes(x = 0, y = 0)) +
  geom_hline(yintercept = 0, linetype = 'dotted', colour = 'grey20') +
  geom_vline(xintercept = 0, linetype = 'dotted', colour = 'grey20') +
  geom_segment(aes(xend = PC1, yend = PC2, colour = wvl_nm ) ) +
  xlab(bquote('PC 1 (' * .(pc1) *'%)')) +
  ylab("") +
  coord_equal() +
  scale_colour_viridis(option = 'magma', name = 'Wvl (nm)', direction = -1) +
 # scale_colour_gradientn(colours = rainbow(20),
  #                       name = 'Wvl (nm)') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
PCA_loading_plot

# add titles
pca_plot2 <- pca_plot + ggtitle('a') +
  theme(plot.title = element_text(size = 20, face = 'bold'))
PCA_loading_plot2 <- PCA_loading_plot + ggtitle('b') +
  theme(plot.title = element_text(size = 20, face = 'bold'))
loadings <- loadings_plot_colours + ggtitle('c') +
  theme(plot.title = element_text(size = 20, face = 'bold'))

pca_panel_plot <- ggarrange(ggarrange(pca_plot2, PCA_loading_plot2, widths = c(2.98,3)),  loadings, nrow = 2)
ggsave('simulated_example/figures/pca_panel.png', plot = pca_panel_plot, width = 8.5, height = 7, dpi = 600)
ggsave('simulated_example/figures/pca_panel.pdf', plot = pca_panel_plot, width = 8.5, height = 7)


### Plot the spectra ----

# Extract the colour for each plant from the PCs
cols <- site_scores %>% 
  dplyr::select(plantID, RGB_col)

# Merge to spectra and plot and add colours
refl_long <- refl_trees %>% 
  rename(Absolute = refl, Normalized = refl_norm) %>% 
  gather(key = refl_type, value = reflectance, -wvl_nm, -scientificName, -plantID)
  
# Draw spectra plot
spectra_plot <- ggplot(refl_long, aes(x = wvl_nm, y = reflectance)) +
  geom_line(aes(group = plantID, colour = plantID)) +
  facet_grid(refl_type ~ scientificName, scales = 'free_y') +
  xlab('Wavelength (nm)') +
  ylab('Reflectance') +
  scale_colour_manual(breaks = cols$plantID, values = cols$RGB_col) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
spectra_plot

# Save plot as PDF and PNG
ggsave('simulated_example/figures/spectra.png', plot = spectra_plot, width = 7, height = 4, dpi = 600)
ggsave('simulated_example/figures/spectra.pdf', plot = spectra_plot, width = 5, height = 5)



##### Simulate landscapes ----

# Create XY positions
len <- 25 # grid length (5 x 5)
pixels <- tibble(X = rep(1:len, len),
                     Y = rep(1:len, each = len),
                     comm = c(rep(rep(1:5, each = 5), 5),
                              rep(rep(6:10, each = 5), 5),
                              rep(rep(11:15, each = 5), 5),
                              rep(rep(16:20, each = 5), 5),
                              rep(rep(21:25, each = 5), 5)) ) 

# Calculate mean XY positions of communities (plots)
plots <- pixels %>% 
  group_by(comm) %>% 
  summarise(X = mean(X),
            Y = mean(Y))

# Get spectra for each plant
spectra <- refl_trees %>%
  ungroup() %>% 
  dplyr::select(plantID, scientificName, wvl_nm, refl_norm) %>% 
  spread(wvl_nm, refl_norm) %>% 
  arrange(scientificName, plantID)


#### Scenario with high beta but low alpha spectral diversity ----
set.seed(9999)
max_beta <- function(x) {
  sp <- sample(unique(spectra$scientificName), size = 1, prob = c(0.6,0.35,0.05))
  pixel_sel <- sample_n(subset(spectra, scientificName == sp), size = nrow(x), replace = T)
  res <- bind_cols(x, pixel_sel)
  return(res)
}

# Get high beta diversity scenario
beta_high_scenario <- pixels %>% 
  group_by(comm) %>% 
  do(max_beta(.)) %>% 
  ungroup()


#### Principal Component Analysis (PCA) of all pixels ----

# Convert data table to wide format for PCA
beta_spectra <- beta_high_scenario %>%
  dplyr::select(`410`:`2400`)


# Perform PCA, type-1 scaling
pca_pixels <- pca_mat(beta_spectra, scaling = 1, p = .99999) 

# Extract the PCA site scores
pixels_scores <- bind_cols(dplyr::select(beta_high_scenario, X, Y, comm, plantID, scientificName),
                          data.frame(pca_pixels$obj) )

# Add colours
col_RGB <- site_scores %>% 
  dplyr::select(plantID, green, red, blue, RGB_col)
beta_high <- pixels_scores %>%
  left_join(col_RGB, by = 'plantID')

# Plot it
max_beta_plot <- ggplot(beta_high, aes(x = X, y = Y)) +
  geom_tile(fill = beta_high$RGB_col, alpha = 0.75) +
  geom_tile(data = plots, fill = NA, colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
  geom_text(data = plots, aes(label = comm), colour = 'black', size = 5) +
  coord_equal() +
  ggtitle(expression('High'~beta*', low'~alpha)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 20) )
max_beta_plot


#### Scenario with low beta but high alpha spectral diversity ----
set.seed(1001)
XY <- sample_n(pixels, size = nrow(pixels), replace = F)
alpha_high <- beta_high %>% 
  dplyr::select(-X, -Y, -comm) %>% 
  bind_cols(XY)

# Plot it
max_alpha_plot <- ggplot(alpha_high, aes(x = X, y = Y)) +
  geom_tile(fill = alpha_high$RGB_col, alpha = 0.75) +
  geom_tile(data = plots, fill = NA, colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
  geom_text(data = plots, aes(label = comm), colour = 'black', size = 5) +
  coord_equal() +
  ggtitle(expression('Low'~beta*', high'~alpha)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 20) )
max_alpha_plot

# Plot the two scenarios next to each other
scenarios_plot <- ggarrange(max_beta_plot, max_alpha_plot)
ggsave('simulated_example/figures/scenarios.png', plot = scenarios_plot, width = 7, height = 4, dpi = 600)
ggsave('simulated_example/figures/scenarios.pdf', plot = scenarios_plot, width = 7, height = 4)


#### Transform scenarios to Raster Bricks ----

# Beta high scenario
beta_high_xy <- beta_high %>%
  dplyr::select(X, Y)
beta_high_PCs <- beta_high %>%
  dplyr::select(PC1:PC14)
beta_high_points <- SpatialPixelsDataFrame(beta_high_xy, beta_high_PCs)
beta_high_brick <- brick(beta_high_points)


# Alpha high scenario
alpha_high_xy <- alpha_high %>%
  dplyr::select(X, Y)
alpha_high_PCs <- alpha_high %>%
  dplyr::select(PC1:PC14)
alpha_high_points <- SpatialPixelsDataFrame(alpha_high_xy, beta_high_PCs)
alpha_high_brick <- brick(alpha_high_points)


### Calculate spectral diversity for both scenarios ----
beta_high_sdiv <- specdiv(beta_high_brick, fact = 5)
alpha_high_sdiv <- specdiv(alpha_high_brick, fact = 5)

# Merge the two alpha and beta SS components tables
ss_beta_high <- beta_high_sdiv$ss %>% 
  mutate(scenario = "'High'~beta*', low'~alpha") %>% 
  dplyr::filter(source != 'gamma')
ss_alpha_high <- alpha_high_sdiv$ss %>% 
  mutate(scenario = "'Low'~beta*', high'~alpha") %>% 
  dplyr::filter(source != 'gamma')
ss_all <- ss_beta_high %>% 
  bind_rows(ss_alpha_high)
  

# Donut plot to show proportion of alpha and beta ----
ss_donut <- ggdonutchart(ss_all, 'prop_gamma', label = 'source', lab.pos = 'out', fill = 'source') +
  facet_wrap(~scenario, labeller = label_parsed) +
   theme_void() +
  scale_fill_brewer(labels = c(expression(alpha), expression(beta)),
                    palette = 'PuOr',
                    name = 'SS') +
  theme(legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        legend.position = 'right',
        legend.text = element_text(size = 16))
ss_donut
ss_donut_no_titles <- ss_donut + theme(strip.text = element_blank(),
                                       legend.position = 'top') 

# Plot Donut SS plot
ggsave('simulated_example/figures/ss_donut.png', plot = ss_donut, width = 6, height = 2.8, dpi = 600)
ggsave('simulated_example/figures/ss_donut.pdf', plot = ss_donut, width = 6, height = 2.8)


# Plot LCSS beta ----
beta_high_lcss <- as_tibble(rasterToPoints(beta_high_sdiv$rasters$beta_lcss)) %>% 
  mutate(scenario = "'High'~beta*', low'~alpha") %>% 
  arrange(y, x)
alpha_high_lcss <- as_tibble(rasterToPoints(alpha_high_sdiv$rasters$beta_lcss)) %>% 
  mutate(scenario = "'Low'~beta*', high'~alpha") %>% 
  arrange(y, x)

# Merge the two tables
lcss_all <- beta_high_lcss %>% 
  bind_rows(alpha_high_lcss) %>% 
  rename(X = x, Y = y)

### Plot it
lcss_plot <- ggplot(lcss_all, aes(x = X, y = Y)) +
  geom_tile(aes(fill = lcss_beta), colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
  # geom_text(data = plots, aes(label = comm), colour = 'black', size = 5) +
  facet_wrap(~scenario, labeller = label_parsed) +
  coord_equal() +
  scale_fill_viridis(direction = -1, option = 'magma',
                     name = expression(LCSS[beta])) +
  theme_void() +
  theme(strip.text.x = element_text(size = 20))
lcss_plot


# Plot LCSD beta ----
beta_high_lcsd <- as_tibble(rasterToPoints(beta_high_sdiv$rasters$beta_lcsd)) %>% 
  mutate(scenario = "'High'~beta*', low'~alpha") %>% 
  arrange(y, x)
alpha_high_lcsd <- as_tibble(rasterToPoints(alpha_high_sdiv$rasters$beta_lcsd)) %>% 
  mutate(scenario = "'Low'~beta*', high'~alpha") %>% 
  arrange(y, x)

# Merge the two tables
lcsd_all <- beta_high_lcsd %>% 
  bind_rows(alpha_high_lcsd) %>% 
  rename(X = x, Y = y)

### Plot it
lcsd_plot <- ggplot(lcsd_all, aes(x = X, y = Y)) +
  geom_tile(aes(fill = lcsd_beta), colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
 #geom_text(data = plots, aes(label = comm), colour = 'white', size = 5) +
  facet_wrap(~scenario, labeller = label_parsed) +
  coord_equal() +
  scale_fill_viridis(direction = -1, option = 'magma',
                     name = expression(LCSD[beta])) +
  theme_void() +
  theme(strip.text.x = element_text(size = 20))
lcsd_plot

# Plot LCSD
ggsave('simulated_example/figures/lcsd.png', plot = lcsd_plot, width = 9, height = 4, dpi = 600)
ggsave('simulated_example/figures/lcsd.pdf', plot = lcsd_plot, width = 9, height = 4)

#Plot LCSS and LCSD together
lcsd_plot2 <- lcsd_plot + 
  theme(strip.text.x = element_blank() )
lcss_lcsd_plot <- ggarrange(lcss_plot, lcsd_plot2,
                            nrow = 2, align = 'hv')
ggsave('simulated_example/figures/lcss_lcsd.png', plot = lcss_lcsd_plot, width = 7, height = 6, dpi = 600)
ggsave('simulated_example/figures/lcss_lcsd.pdf', plot = lcss_lcsd_plot, width = 7, height = 6)


### Plot alpha spectral diversity ----

beta_high_alpha <- as_tibble(rasterToPoints(beta_high_sdiv$rasters$alpha_sdiv)) %>% 
  mutate(scenario = "'High'~beta*', low'~alpha") %>% 
  arrange(y, x)
alpha_high_alpha <- as_tibble(rasterToPoints(alpha_high_sdiv$rasters$alpha_sdiv)) %>% 
  mutate(scenario = "'Low'~beta*', high'~alpha") %>% 
  arrange(y, x)

# Merge the two tables
alpha_all <- beta_high_alpha %>% 
  bind_rows(alpha_high_alpha) %>% 
  rename(X = x, Y = y)

### Plot it
alpha_plot <- ggplot(alpha_all, aes(x = X, y = Y)) +
  geom_tile(aes(fill = sdiv_alpha), colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
  #geom_text(data = plots, aes(label = comm), colour = 'white', size = 5) +
  facet_wrap(~scenario, labeller = label_parsed) +
  coord_equal() +
  scale_fill_viridis(direction = -1, option = 'magma',
                     name = expression(SD[alpha])) +
  theme_void() +
  theme(strip.text.x = element_text(size = 20))
alpha_plot

# Plot LCSD
ggsave('simulated_example/figures/sdiv.png', plot = alpha_plot, width = 9, height = 4, dpi = 600)
ggsave('simulated_example/figures/sdiv.pdf', plot = alpha_plot, width = 9, height = 4)


### Plot feature contribution to spectral alpha diversity ----

beta_high_fcsd <- as_tibble(rasterToPoints(beta_high_sdiv$rasters$alpha_fcsd)) %>% 
  dplyr::select(x:PC2) %>% 
  gather(key = 'PC', value = 'value', -x, -y) %>% 
  mutate(scenario = "'High'~beta*', low'~alpha") %>% 
  arrange(y, x)
alpha_high_fcsd <- as_tibble(rasterToPoints(alpha_high_sdiv$rasters$alpha_fcsd)) %>% 
  dplyr::select(x:PC2) %>% 
  gather(key = 'PC', value = 'value', -x, -y) %>% 
  mutate(scenario = "'Low'~beta*', high'~alpha") %>% 
  arrange(y, x)

# Merge the two tables
fcsd_all <- beta_high_fcsd %>% 
  bind_rows(alpha_high_fcsd) %>% 
  rename(X = x, Y = y)

### Plot it
fcsd_plot <- ggplot(fcsd_all, aes(x = X, y = Y)) +
  geom_tile(aes(fill = value), colour = 'black', width = sqrt(len), height = sqrt(len), size = 1) +
  #geom_text(data = plots, aes(label = comm), colour = 'white', size = 5) +
  facet_grid(PC ~ scenario, labeller = label_parsed) +
  coord_equal() +
  scale_fill_viridis(direction = -1, option = 'magma',
                     name = expression(FCSD[alpha])) +
  theme_void() +
  theme(strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20,
                                    angle = 270))
fcsd_plot

# Plot LCSD
ggsave('simulated_example/figures/fcsd.png', plot = fcsd_plot, width = 9, height = 6, dpi = 600)
ggsave('simulated_example/figures/fcsd.pdf', plot = fcsd_plot, width = 9, height = 6)
