library(mclust)
merged_filter_exclude_epi<-merged_filter_exclude_epi2
test_dist(c('Macrophage','Neutrophil','T_cell','B_Plasma'))->merged_filter_exclude_epi$min_dis_to_immune

# distances vector (in microns)
d <- merged_filter_exclude_epi$min_dis_to_immune[which(merged_filter_exclude_epi$cell_type_l1_all%in%c('Endothelium','Fibroblasts'))]
d<-d[which(!is.na(d))]
d<-d[which(d<1000)]

ggplot(data.frame(min_dis_to_immune=d),aes(x = min_dis_to_immune/4)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  theme_classic()+
  theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))->p

set.seed(100)
fit <- Mclust(d, G = 2:4)
means <- fit$parameters$mean        # mean distance of each component
ord   <- order(means)               # sort components by mean distance

prox_comp <- ord[1]                 # smallest mean = proximal cluster
far_comp  <- ord[length(ord)]       # largest mean = distant/background

grid <- seq(min(d), max(d), length.out = 1000)

# For 1D, newdata must be a matrix with 1 column
pred <- predict(fit, newdata = matrix(grid, ncol = 1))
post <- pred$z                      # posterior probs: rows = grid points, cols = components

p_prox <- post[, prox_comp]
p_far  <- post[, far_comp]

idx    <- which.min(abs(p_prox - p_far))  # where |Pprox - Pfar| is smallest
d_star <- grid[idx]

d_star
ggplot(data.frame(min_dis_to_immune=d),aes(x = min_dis_to_immune/4)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  geom_vline(xintercept =d_star/4)+
  theme_classic()+
  theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))->p
ggsave('p.png',p,width = 5.5,height = 5.5,dpi = 300)

