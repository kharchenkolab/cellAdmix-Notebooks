# Requires notebook 4

library(grid)
library(ComplexHeatmap)

hm1 <- CachePath('nsclc_scaled_dat11', 'bridge_hmap.rds') %>% readRDS()

cairo_pdf(OutputPath('nsclc_bridge_frem.pdf'), family="DejaVu Sans", width=4.5, height=3.5)
hm1
dev.off()
