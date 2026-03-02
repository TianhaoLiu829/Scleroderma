# install.packages("dbscan") # if needed
library(dbscan)
library(dplyr)
library(tidyr)
library(rlang)
library(purrr)
# df1: data frame with columns x, y  (spots of interest)
# df2: data frame with columns x, y, and a cell-type column (default "cell_type")
# targets: character vector of cell types to sum; can be NULL or character(0)
# id_col: optional column name in df1 to carry through as an identifier
celltype_composition_by_radius <- function(df1, df2, radius = 100,
                                           celltype_col = "cell_type",
                                           targets = character(0),
                                           id_col = NULL) {
  stopifnot(all(c("x","y") %in% names(df1)),
            all(c("x","y", celltype_col) %in% names(df2)))
  
  # Queries and reference matrices
  q <- as.matrix(df1[, c("x","y")])
  r <- as.matrix(df2[, c("x","y")])
  
  # Fixed-radius neighbors of df1 points in df2
  nn <- frNN(x = r, eps = radius, query = q, sort = TRUE)
  
  # Long table of neighbor pairs
  pairs_long <- lapply(seq_along(nn$id), function(i) {
    if (length(nn$id[[i]]) == 0) return(NULL)
    tibble(
      df1_idx  = i,
      df2_idx  = nn$id[[i]],
      dist     = nn$dist[[i]]
    )
  }) |> list_rbind()
  
  # If no neighbors at all, return a zero-filled skeleton
  all_celltypes <- df2[[celltype_col]] |> as.character() |> unique() |> sort()
  if (is.null(pairs_long)) {
    out <- tibble(df1_idx = seq_len(nrow(df1)),
                  n_neighbors = 0)
    # add per-CT props (all zero)
    if (length(all_celltypes)) {
      zeros <- as_tibble(matrix(0, nrow(df1), ncol = length(all_celltypes),
                                dimnames = list(NULL, paste0("prop_", all_celltypes))))
      out <- bind_cols(out, zeros)
    }
    out$prop_sum_targets <- 0
    if (!is.null(id_col) && id_col %in% names(df1)) {
      out[[id_col]] <- df1[[id_col]]
      out <- out |> relocate(any_of(id_col), .before = 1)
    }
    return(out)
  }
  
  # Attach cell type to df2 indices
  pairs_long <- pairs_long |>
    mutate(!!celltype_col := df2[[celltype_col]][df2_idx] |> as.character())
  
  # Count by df1 spot × cell type; compute proportions
  counts <- pairs_long |>
    count(df1_idx, !!sym(celltype_col), name = "n") |>
    group_by(df1_idx) |>
    mutate(n_neighbors = sum(n), prop = n / n_neighbors) |>
    ungroup()
  
  # Ensure zero rows for missing cell types per df1 (complete)
  counts_complete <- counts |>
    complete(df1_idx, !!sym(celltype_col) := all_celltypes,
             fill = list(n = 0, prop = 0, n_neighbors = 0))
  
  # One row per df1 with wide prop columns
  props_wide <- counts_complete |>
    mutate(ct_col = paste0("prop_", .data[[celltype_col]])) |>
    select(df1_idx, n_neighbors, ct_col, prop) |>
    distinct() |>
    pivot_wider(names_from = ct_col, values_from = prop, values_fill = 0) |>
    group_by(df1_idx) |>
    summarise(across(everything(), ~ dplyr::first(.x)), .groups = "drop")
  
  # Add sum of target CT proportions (0 if targets empty or absent)
  if (length(targets) > 0) {
    target_cols <- paste0("prop_", targets)
    present <- intersect(target_cols, names(props_wide))
    props_wide <- props_wide |>
      mutate(prop_sum_targets = if (length(present)) rowSums(across(all_of(present)), na.rm = TRUE) else 0)
  } else {
    props_wide <- props_wide |>
      mutate(prop_sum_targets = 0)
  }
  
  # Attach optional ID and return
  if (!is.null(id_col) && id_col %in% names(df1)) {
    props_wide[[id_col]] <- df1[[id_col]][props_wide$df1_idx]
    props_wide <- props_wide |> relocate(any_of(id_col), .before = 1)
  }
  
  props_wide
}

celltype_composition_by_radius(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium'),c('x','y')],merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all!='Endothelium'),c('x','y','cell_type_l1_all')],celltype_col = "cell_type_l1_all",targets = c('B_Plasma','T_cell','Macrophage','Neutrophil'),radius = 220)->out
out$df1_idx<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')][out$df1_idx]
out$sample<-substr(out$df1_idx,1,3)
out$condition<-merged_filter_exclude_epi$condition[match(out$sample,merged_filter_exclude_epi$orig.ident)]
#pre/post treatment compare
out_treat<-out[which(out$condition=='SSC'&out$sample!='A_2'),]
out_treat$condition<-'pre'
out_treat$condition[grepl('2',out_treat$sample)]<-'post'
pre<-out_treat$prop_sum_targets[which(out_treat$condition=='pre')]
post<-out_treat$prop_sum_targets[which(out_treat$condition=='post')]
pre<-out_treat$prop_sum_targets[which(out_treat$condition=='pre'&grepl('K',out_treat$sample))]
post<-out_treat$prop_sum_targets[which(out_treat$condition=='post'&grepl('K',out_treat$sample))]
t.test(pre[which(pre>0)],post[which(post>0)])

out_treat[which(out_treat$prop_sum_targets>0&grepl('H',out_treat$sample)),]->plot
plot$condition<-factor(plot$condition,levels = c('pre','post'))
plot<-plot[-25,]
ggplot(plot, aes(x = condition, y = prop_sum_targets, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.5) +
  geom_boxplot(
    width = 0.15,          # thinner than violin
    color = "black",
    fill = "white"
  ) +
  # 3️⃣ Jitter layer — show individual points
  geom_jitter(
    width = 0.1,
    size = 0.5,
    color = "black",
    alpha = 0.8
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("pre" = "#E9967A", "post" = "#66C2A5")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13)
  )->p

#LS/SSc compare
out<-out[grep('1',out$sample),]
out<-out[which(out$sample!='G_1'),]
out$prop_sum_targets[which(out$condition=='SSC')]->SSC
out$prop_sum_targets[which(out$condition=='LS')]->LS
wilcox.test(SSC[which(SSC>0)],LS[which(LS>0)])
library(ggsignif)
ggplot(out[which(out$prop_sum_targets>0),], aes(x = condition, y = prop_sum_targets, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.5) +
  geom_boxplot(
    width = 0.15,          # thinner than violin
    color = "black",
    fill = "white"
  ) +
  # 3️⃣ Jitter layer — show individual points
  geom_jitter(
    width = 0.1,
    size = 0.5,
    color = "black",
    alpha = 0.8
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("LS" = "#E9967A", "SSC" = "#66C2A5")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13)
  )->p
ggsave('p.png',p,width = 5,height = 6,dpi = 300)
