

#compute the density
df_endo<-as.data.frame(seu_ssc_sep@meta.data)
dens <- density(seu_ssc_sep$pseudotime[which(seu_ssc_sep$cell_type_l1_all=='endo_fib')], bw = "nrd0",from = min(seu_ssc_sep$pseudotime), to   = max(seu_ssc_sep$pseudotime))
density_df <- data.frame(
  pseudotime = dens$x,
  density = dens$y
)
density_df$density_scaled <- density_df$density / max(density_df$density)


#plot out the figure 5 correlations between distance to fibrotic endothelial cells and pseudo times
seu_ssc_sep$seurat_clusters2<-seu_ssc_all$seurat_clusters[match(colnames(seu_ssc_sep),colnames(seu_ssc_all))]
seu_ssc_sep@meta.data$cell_type_l1_all<-ifelse(seu_ssc_sep$seurat_clusters2==2,'endo_fib','Endothelium')
seu_ls_sep$seurat_clusters2<-seu_ssc_all$seurat_clusters[match(colnames(seu_ls_sep),colnames(seu_ssc_all))]
seu_ls_sep@meta.data$cell_type_l1_all<-ifelse(seu_ls_sep$seurat_clusters2==2,'endo_fib','Endothelium')

library(dplyr)
df<-seu_ssc_sep@meta.data
df<-df[which(!is.na(df$min_dis_to_endo_fib)&df$min_dis_to_endo_fib<4000),]
df$mature_numeric <- ifelse(df$cell_type_l1_all == "endo_fib", 1, 0)
df$min_dis_to_endo_fib<-df$min_dis_to_endo_fib/4
library(mgcv)
gam_fit <- gam(mature_numeric ~ s(pseudotime, k = 50), 
               data = df, 
               family = binomial)
pt_grid <- seq(min(df$pseudotime), max(df$pseudotime), length.out = 2000)
pred_df <- data.frame(
  pseudotime = pt_grid,
  frac_mature = predict(gam_fit, newdata = data.frame(pseudotime = pt_grid), type = "response")
)
strip_height <- 0.05 * max(df$min_dis_to_endo_fib)
tile_width <- (max(df$pseudotime, na.rm = TRUE) -
                 min(df$pseudotime, na.rm = TRUE)) / 2000


ggplot(df, aes(y = min_dis_to_endo_fib, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  theme_classic(base_size = 14) +
  geom_tile(
    data = pred_df,
    aes(
      x = pseudotime,
      y= -strip_height,
      fill = frac_mature
    ),
    height = strip_height,
    width = (max(df$pseudotime) - min(df$pseudotime)) / 2000   # matches the grid spacing
  ) +
  scale_y_reverse(
    name = "Left axis (High → Low)",
    sec.axis = sec_axis(~ max(df$min_dis_to_endo_fib) + min(seu_ssc@meta.data$min_dis_to_endo_fib) - ., name = "Right axis (Low → High)")
  ) +
  stat_cor(method = "pearson", label.x = 0.5)+
  scale_fill_gradient(limits = c(0, 1),low = "white", high = "red", name = "% mature (smooth)")->p
ggsave('p.png',p,width = 8.5,height = 6,dpi = 300)








