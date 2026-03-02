min<-min(min(df_num_inf_ls_streng),min(df_num_inf_ssc_streng))
max<-max(max(df_num_inf_ls_streng),max(df_num_inf_ssc_streng))
max_num<-max(max(df_num_inf_ls),max(df_num_inf_ssc))
dot_breaks<-c(0, 1, 5, 10, 20,30,40)
df_streng<-df_num_inf_ls_streng
df_inf<-df_num_inf_ls

#df_streng<-df_num_inf_ls_streng
#df_inf<-df_num_inf_ls


df <- as.data.frame(as.table(as.matrix(df_inf))) %>%
  rename(sender = Var1, receiver = Var2, n_signals = Freq) %>%
  left_join(
    as.data.frame(as.table(as.matrix(df_streng))) %>%
      rename(sender = Var1, receiver = Var2, comm_strength = Freq),
    by = c("sender", "receiver")
  ) %>%
  mutate(
    sender = as.character(sender),
    receiver = as.character(receiver),
    n_signals = as.numeric(n_signals),
    comm_strength = as.numeric(comm_strength)
  )

# Optional: remove pairs with no signals (keeps plot less cluttered)
df_plot <- df %>% filter(n_signals >= 0)


df_plot$sender<-factor(df_plot$sender,levels = c('Macrophage','T_cell','B_Plasma'))
df_plot$receiver<-factor(df_plot$receiver,levels = c('inf_endo','inf_fib'))

p <- ggplot(df_plot, aes(x =  sender, y =receiver)) +
  geom_point(aes(size = n_signals, color = comm_strength), alpha = 0.9) +
  scale_size(
    range = c(1.5, 8),                  # visual size in mm
    limits = c(0,max_num),
    breaks = dot_breaks,
    name = "# significant\nsignals"
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "orange",
    high = "red",
    midpoint = (min + max) / 2,
    limits = c(min, max),
    oob = scales::squish,   # values < min or > max are clipped
    name = "Summed\ncommunication\nstrength"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Receiver", y = "Sender")
ggsave('p_non_inf.png',p,width = 6,height = 2.5,dpi = 300)
ggsave('p_non_inf1.png',p,width = 6,height = 5,dpi = 300)
ggsave('p_inf_ssc.png',p,width = 6,height = 2.5,dpi = 300)
ggsave('p_inf_ssc1.png',p,width = 6,height = 5,dpi = 300)
