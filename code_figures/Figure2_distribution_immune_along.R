merged_filter_exclude_epi2 <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")

#the distribution of immune cells along the distance to different types of parenchymal cells
test_dist(c('Fibroblasts','Endothelium'))->merged_filter_exclude_epi$min_dis_to_fib_endo
merged_filter_exclude_epi_filter<-merged_filter_exclude_epi2
merged_filter_exclude_epi_filter@meta.data %>%
  filter(cell_type_l1_all %in% c('B_Plasma','Macrophage','T_cell','Neutrophil'))->df1
df_long <- df1[,grep('min_dis_to',colnames(df1))] %>%
  mutate(rowname = rownames(df1[,grep('min_dis_to',colnames(df1))])) %>%    
  pivot_longer(cols = -rowname, 
               names_to = "original_col",
               values_to = "value")
df_long<-df_long[which(!is.na(df_long$value)),]
df_long[which(df_long$original_col%in%c('min_dis_to_endo','min_dis_to_fib','min_dis_to_gland','min_dis_to_kera')),]->target_df
target_df %>% group_by(original_col) %>% summarise(pct_less_25 = mean(value < 220) * 100)
ggplot(target_df[which(target_df$value<4000),], aes(x = value/4, 
                                                    color = original_col, 
                                                    fill = original_col)) +
  geom_density(alpha = 0.3) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  theme_classic() +
  geom_vline(xintercept = 33)+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )->p
ggsave('p.png',p,width = 6.5,height = 5.5,dpi = 300)
ggplot(target_df[which(target_df$value<4000),], aes(x = value/4, colour = original_col)) +
  stat_ecdf(geom = "step", size = 1) +
  labs(
    x = "Distance threshold",
    y = "Cumulative density  P(distance ≤ threshold)",
    colour = "original_col"
  ) +
  theme_classic()->p
ggsave('p.png',p,width = 6.5,height = 5.5,dpi = 300)

#make the pie plot in figure 2
merged_filter_exclude_epi2->merged_filter_exclude_epi
merged_filter_exclude_epi@meta.data[,c('min_dis_to_gland','min_dis_to_kera','min_dis_to_fib','min_dis_to_endo')]->df
df<-df[which(merged_filter_exclude_epi$cell_type_l1_all%in%c('T_cell','B_Plasma','Neutrophil','Macrophage')),]
df<-df[-grep('G_1',rownames(df)),]
unlist(apply(df,1,function(x) colnames(df)[which.min(x)]))->a
#pie plot based on nearest distance
ggplot(data.frame(nearest_cell_com=a), aes(x = "", fill = nearest_cell_com)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Pie chart of labels") +
  theme_void()

#confusing figure showing individual distance
data.frame(min_label=a,min_value=apply(df, 1, min)/4)->min_df
ggplot(min_df, aes(x = min_label, y = log2(min_value))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # remove extreme outliers if desired
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +  # show raw points
  theme_classic() +
  labs(x = "Label", y = "Value")

#another plot do not consider the min
library(tidyverse)
df_long <- df %>%
  mutate(rowname = rownames(df)) %>%    
  pivot_longer(cols = -rowname, 
               names_to = "original_col",
               values_to = "value")
ggplot(df_long, aes(x = original_col, y = log2(value),fill = original_col)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA)+  # remove extreme outliers if desired
  theme_classic() +
  theme(axis.text.y = element_text(size = 13))+
  labs(x = "Label", y = "Value")->p
ggsave('p.png',p,width = 5.5,height = 4,dpi = 300)


unlist(apply(df,1,function(x) {
  if (min(x)<130 ){
    return(colnames(df)[which.min(x)])
  } else {
    return('other')
  }
}))->a

#check the distribution of min distance to one of fib,endo,kera,gland
ggplot(data.frame(value=apply(df, 1, min)/4), aes(x = value)) +
  geom_density(alpha = 0.3) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  theme_classic() +
  geom_vline(xintercept = 32.5)+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )->p





breaks<-seq(0, 4000, by = 10)
df_cdf <- target_df[which(target_df$value<4000),] %>%
  group_by(original_col) %>%
  summarise(
    thresh = list(breaks),
    cdf = list(sapply(breaks, function(t) mean(value <= t))),
    .groups = "drop"
  ) %>%
  unnest(c(thresh, cdf))
df_combined <- df_cdf %>%
  pivot_wider(names_from = original_col, values_from = cdf) %>%
  mutate(
    derived = (min_dis_to_endo + min_dis_to_fib) -
      (min_dis_to_gland + min_dis_to_kera)
  )

merged_filter_exclude_epi_filter@meta.data %>%
  filter(cell_type_l1_all %in% c('Keratinocyte'))->df2
merged_filter_exclude_epi_filter@meta.data %>%
  filter(cell_type_l1_all %in% c('Glandular_epithelium'))->df3


dfs <- list(
  sample1 = df1[,c(22,23)],
  sample2 = df2[,c(22,23)],
  sample3 = df3[,c(22,23)]
)

# Bind them together and add a column "dataset"
big_df <- bind_rows(dfs, .id = "dataset")

ggplot(big_df, aes(x = min_dis_to_fib_endo/4, 
                   color = dataset, 
                   fill = dataset)) +
  geom_density(alpha = 0.3) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  geom_vline(xintercept = 25) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 1),
    axis.text.y = element_text(size = 1)
  )
ggplot(df,aes(x = min_dis_to_fib_endo/4)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density of Fibroblasts vs Distance",
       x = "Distance",
       y = "Density") +
  geom_vline(xintercept =25) +
  theme(axis.text.x = element_text(size=100),axis.text.y = element_text(size=100))+
  theme_classic()->p
ggsave('p1.png',p,width = 6,height = 6,dpi = 300)
merged_filter_exclude_epi<-merged_filter_exclude_epi[which(merged_filter_exclude_epi$orig.ident!='G_1'),]
merged_filter_exclude_epi@meta.data$min_dis_to_fib[merged_filter_exclude_epi@meta.data$cell_type_l1_all%in% c('B_Plasma','Macrophage','T_cell','Neutrophil')]->dis
density(dis[which(!is.na(dis))])->den
den$x[which.max(den$y)]
den$y[den$x>54&den$x<55]
