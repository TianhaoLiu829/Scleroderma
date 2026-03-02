readRDS('/ix1/wchen/liutianhao/result_l1.rds')->result_l1
readRDS('/ix1/wchen/liutianhao/result_l1_8um.rds')->result_l1_8um
result_l1->object
prop_table<-list()
for (sample in names(object)) {
  as.character(object[[sample]]@meta.data$spot_class)->class
  class<-class[which(class!='reject')]
  class<-factor(class,levels = c('singlet','doublet_certain','doublet_uncertain'))
  prop.table(table(class))->prop_table[[sample]]
}
as.data.frame(do.call(rbind,prop_table))->prop_table

library(tidyverse)

df_long <- prop_table %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )
df_long$variable<-factor(df_long$variable,levels = c('singlet','doublet_certain','doublet_uncertain'))
summary_df <- df_long %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n)
  )

ggplot() +
  # Bars
  geom_col(
    data = summary_df,
    aes(x = variable, y = mean),
    fill = "skyblue", width = 0.7
  ) +
  
  # Error bar (SE)
  geom_errorbar(
    data = summary_df,
    aes(x = variable, y = mean, ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  
  # Individual sample points
  geom_jitter(
    data = df_long,
    aes(x = variable, y = value),
    width = 0.1, size = 2, alpha = 0.7
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15,angle = 45, hjust = 1),axis.text.y = element_text(size = 15))+
  labs(x = "", y = "Value")->p
ggsave('p.png',p,width = 7,height = 7,dpi = 300)




as.data.frame(skin_ST_combine_harmony@meta.data[,c('orig.ident','spot_class')])->df
df_long<-df %>% group_by(orig.ident,spot_class) %>% summarise(n = n(), .groups = "drop_last") %>% mutate(prop = n / sum(n))
df_long<-df_long[which(df_long$spot_class%in%c('singlet','doublet_certain','doublet_uncertain')),]
df_long$spot_class<-factor(df_long$spot_class,levels = c('singlet','doublet_certain','doublet_uncertain'))
summary_df <- df_long %>%
  group_by(spot_class) %>%
  summarise(
    mean = mean(prop, na.rm = TRUE),
    sd   = sd(prop, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n)
  )
ggplot() +
  # Bars
  geom_col(
    data = summary_df,
    aes(x = spot_class, y = mean),
    fill = "skyblue", width = 0.7
  ) +
  
  # Error bar (SE)
  geom_errorbar(
    data = summary_df,
    aes(x = spot_class, y = mean, ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  
  # Individual sample points
  geom_jitter(
    data = df_long,
    aes(x = spot_class, y = prop),
    width = 0.1, size = 2, alpha = 0.7
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15,angle = 45, hjust = 1),axis.text.y = element_text(size = 15))+
  labs(x = "", y = "Value")->p
ggsave('p.png',p,width = 7,height = 7,dpi = 300)


