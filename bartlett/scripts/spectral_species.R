####### Apply our method to spectral species-----

cube_data <- getValues(cube_pc_scaling1)
i <- which(!is.na(cube_data))
cube_data <- na.omit(cube_data)

## kmeans classification with 10 'spectral species'
E <- kmeans(cube_data, 10, iter.max = 100, nstart = 10)
kmeans_raster <- raster(cube_pc_scaling1)
kmeans_raster[i] <- E$cluster
plot(kmeans_raster)

# Plot 10 spectral species
spectral_species_plot <- gplot(kmeans_raster) +
  geom_raster(aes(fill = as.factor(value))) + 
  theme_void() + 
  ggtitle('10 spectral species') +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none') + 
  scale_fill_brewer(palette = 'RdYlBu')
spectral_species_plot
ggsave('10_spectral_species.png')

# Multi-panel figure showing PCs and 10 spectral species
png('bartlett/figures/PCs_ten_spectral_species.png', width = 2.5, height = 4, res = 1200, units = 'in')
plot_grid(beta_PC_scaling1, spectral_species_plot, align = 'hv', nrow = 1)
dev.off()

