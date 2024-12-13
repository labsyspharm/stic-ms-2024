# If we want to see correlation with gene signature, then we have to z-transformed our data
# first we need to read the gene signature (from published gene set)

gene_sets <- list(
  #Immune_cells = read_excel("Gene_signatures_ratio_calc.xlsx", "Immune_cells"),
  IFN_IRDS = read_excel("Gene_signatures_ratio_calc.xlsx", "IFN_IRDS")
) %>%
  map(
    \(x) pivot_longer(x, everything(), names_to = "gene_set", values_to = "gene")
  ) %>%
  bind_rows(.id = "gene_set_group") %>%
  drop_na() %>%
  mutate(
    gene_unique = paste(gene, gene_set, sep = "_")
  )

gene_sets %>%
  dplyr::count(gene) %>%
  arrange(desc(n))

setdiff(gene_sets$gene, count_matrix$TargetName)
# 3 genes in gene set not in count matrix. Filtered out during QC

#only chosing the columns matched with count matrix file
selected_samples <- ROI_file %>%
  filter(
    sample_ID_geomx %in% colnames(count_matrix)
  )

selected_gene_sets <- gene_sets %>%
  filter(
    gene %in% count_matrix$TargetName, 
    gene_set %in%c("IRDS") # we can add multiple sets
  )

mat <- count_matrix_mat_df %>%
  inner_join(
    selected_gene_sets,
    by = c("TargetName" = "gene")
  ) %>%
  select(gene_unique, all_of(selected_samples$sample_ID_geomx)) %>%
  column_to_rownames("gene_unique") %>%
  as.matrix() %>% {
    .[selected_gene_sets$gene_unique,]
  }

#first log transformed and then normalised

scaled_mat <- mat %>%
  log10() %>%
  t() %>%
  scale() %>%
  t()


# prior to correlation with a gene set, we calculate mean_expression of that gene set/sample_ID
#dataset is z-transformed and scaled
#as_tibble Creates a new column called "gene_unique" that contains the row names of the original matrix

mean_expression_df <- scaled_mat %>%
  as_tibble(rownames = "gene_unique") %>%
  pivot_longer(-gene_unique, names_to = "sample_id", values_to = "count") %>%
  inner_join(
    selected_gene_sets,
    by = "gene_unique"
  ) %>%
  group_by(gene_set, sample_id) %>%
  summarize(mean_count = mean(count), .groups = "drop")

#Then use this mean_expression df to move forward with the correlation


gene_strata_plots <- mean_expression_df %>%
  power_inner_join(
    combined_geomxROI %>%
      distinct(
        sample_ID_geomx, Cat, Cell_Types
      ),
    by = c("sample_id" = "sample_ID_geomx"),
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  group_by(gene_set, Cell_Types) %>%
  summarize(
    p = list(
      ggplot(
        cur_data(),
        aes(Cat, mean_count)
      ) +
        geom_boxplot(outliers = FALSE, fill = NA) +
        geom_quasirandom()
    ),
    .groups = "drop"
  )

#Group the data by Cat and Cell_Types
#Create a new binary category column bin based on whether mean_count is above or below the median

gene_strata <- mean_expression_df %>%
  power_inner_join(
    combined_geomxROI %>%
      distinct(
        sample_ID_geomx, Cat, Cell_Types
      ),
    by = c("sample_id" = "sample_ID_geomx"),
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  group_by(gene_set, Cat, Cell_Types) %>%
  filter(n() >= 4) %>%
  mutate(
    bin = cut(mean_count, breaks = c(-Inf, median(mean_count), Inf), labels = c("low", "high"))
    ) %>%
  ungroup()

RNA_protein_relation <- gene_strata %>%
  power_inner_join(
    combined_geomxROI %>% 
      distinct(sample_ID_geomx, GroupCount, mean_CD8ap),  #here other markers can be added after mean_CD8ap (CD8+ T cells)
    by = c("sample_id" = "sample_ID_geomx"),
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  pivot_longer(
    c(mean_CD8ap),  #here we have to add markers which we added before
    names_to = "cycif_marker", values_to = "proportion"
  ) %>%
  group_by(cycif_marker) %>%
  mutate(
    cell_count = GroupCount * proportion,
    zero_replacement = min(proportion[proportion > 0]) * .5,
    proportion_zero_replaced = if_else(proportion == 0, zero_replacement, proportion)
  ) %>%
  ungroup()

cat_order <- c("STIC.I", "STIC.C", "Inv Cancer")


#Now plot the relation
p <- RNA_protein_relation %>%
  filter(Cat %in% cat_order) %>%
  filter(Cell_Types %in% c("epithelial")) %>%
  mutate(
    across(Cat, \(x) factor(x, levels = cat_order))
  ) %>%
  ggplot(
    aes(Cat, proportion_zero_replaced, color = bin)
  ) +
  geom_hline(
    aes(yintercept = zero_replacement),
    data = \(x) distinct(x, Cat, gene_set, zero_replacement),
    linetype = "dashed"
  ) +
  geom_quasirandom(
    dodge.width = .75
  ) +
  geom_boxplot(outliers = FALSE, fill = NA, position = position_dodge2(width = .5, padding = 0)) +
  facet_grid(cycif_marker~gene_set) +
  #scale_y_log10(labels = scales::label_percent()) +
  scale_y_log10()+
  theme_light() +
  labs(y = "cells expressing marker") +
  ggtitle("title of the plot")
p


RNA_protein_fisher_res <- RNA_protein_relation %>%
  group_by(
    gene_set, Cell_Types, Cat
  ) %>%
  summarize(
    cont_table = cur_data() %>%
      transmute(
        bin,
        proportion_greater_zero = if_else(proportion > 0, "greater", "zero") %>%
          factor(levels = c("greater", "zero"))
      ) %>%
      table() %>%
      list(),
    .groups = "drop"
  ) %>%
  mutate(
    res = map(
      cont_table,
      \(x) fisher.test(as.matrix(x))
    ),
    res_df = map(res, broom::tidy)
  )


RNA_protein_fisher_res_df <- RNA_protein_fisher_res %>%
  select(-res) %>%
  unnest(res_df)

# Check contingency table
RNA_protein_fisher_res_df$cont_table[[6]]


