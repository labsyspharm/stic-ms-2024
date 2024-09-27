library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(powerjoin)
library(broom)
library(broom.mixed)
library(emmeans)


masterfile <- read_csv("masterfile_join_slidename_ROI.csv")

# To count the number of cells, group by these columns
group_cols <- c("slideName", "BRCA", "Stage", "Incidental", "Cat", "ROIcat")
# Select columns with positive / negative calls
marker_cols <- colnames(masterfile) %>%
  str_subset(".+[pn]$") %>%
  setdiff(c("Orientation", "Region", "ROI_n"))

marker_cols

marker_combinations <- tribble(
  ~marker_1, ~marker_2,
  "pSTAT1p", "pTBK1p",
  "pTBK1p", "HLA_Ap",
  "pTBK1p", "HLA_Ep",
  "pSTAT3p", "pTBK1p"
)

make_contingency_tables <- function(
  marker_1, marker_2
) {
  marker_1_sym <- sym(marker_1)
  marker_2_sym <- sym(marker_2)
  masterfile %>%
    count(
      across(all_of(c(group_cols, marker_1, marker_2)))
    ) %>%
    # Rename columns x, y to marker_status_1, marker_status_2
    rename(
      marker_status_1 = !!marker_1_sym,
      marker_status_2 = !!marker_2_sym
    )
}

cell_populations <- marker_combinations %>%
  mutate(
    data = map2(
      marker_1, marker_2,
      make_contingency_tables
    )
  )

lesion_type_order <- c(
  "Fallopian Tube",
  "Fimbriae",
  "p53 signature",
  "STIC",
  "Carcinoma"
)


cell_populations_pivot <- cell_populations %>%
  mutate(
    data = map(
      data,
      \(x) pivot_wider(
        x,
        names_from = marker_status_2,
        names_prefix = "marker_2_",
        values_from = n,
        values_fill = 0
      )
    )
  ) %>%
  unnest(data)  %>%
  mutate(
    marker_status_1 = factor(as.character(marker_status_1)) %>%
      fct_relevel("0"),
  ) %>%
  group_nest(marker_1, marker_2, Incidental) %>%
  mutate(
    data = map(
      data,
      \(x)
        x %>%
        filter(
          ROIcat %in% lesion_type_order
        ) %>%
        mutate(
          ROIcat = factor(
            ROIcat,
            levels = intersect(lesion_type_order, unique(ROIcat))
        )
      )
    )
  )

cell_population_fisher_res <- cell_populations_pivot %>%
  mutate(
    data = map(
      data,
      \(x) group_by(
        x,
        Cat, ROIcat,
        marker_status_1
      ) %>%
        summarize(
          across(c(marker_2_0, marker_2_1), sum),
          .groups = "drop"
        ) %>%
        group_nest(Cat, ROIcat)
    )
  ) %>%
  unnest(data) %>%
  mutate(
    cont_table = map(
      data,
      \(x) cbind(x$marker_2_0, x$marker_2_1) %>%
        as.matrix()
    ),
    fisher_res = map(
      cont_table,
      fisher.test
    ),
    fisher_df = map(
      fisher_res,
      tidy
    )
  )

cell_population_fisher_res_long <- cell_population_fisher_res %>%
  select(marker_1, marker_2, Incidental, Cat, ROIcat, fisher_df) %>%
  unnest(fisher_df)

# to confirm the fisher exact test is identical with glm

cell_population_glm_res <- cell_population_fisher_res %>%
  mutate(
    glm_res = map(
      data,
      \(x) glm(
        cbind(marker_2_1, marker_2_0) ~ marker_status_1,
        data = x,
        family = binomial
      )
    ),
    glm_df = map(
      glm_res,
      tidy
    )
  )

cell_population_glm_res_long <- cell_population_glm_res %>%
  select(marker_1, marker_2, Incidental, Cat, ROIcat, glm_df) %>%
  unnest(glm_df) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    estimate_or = exp(estimate)
  )


#results from glm and fisher was identical, now we can run this glm with random effect by patient ID
library(lme4)
library(lmerTest)
cell_populations_models <- cell_populations_pivot %>%
  mutate(
    model = map(
      data,
      \(x) glmer(
        cbind(marker_2_1, marker_2_0) ~ marker_status_1 * ROIcat + (1 | slideName),
        data = x,
        family = binomial
      )
    ),
    contrasts_across_stages = map(
      model,
      \(x) emmeans(
        x,
        specs = ~ ROIcat * marker_status_1
      ) %>%
        contrast(
          interaction = TRUE,
          method = "revpairwise", adjust = "bh", type = "response"
        ) %>%
        as_tibble() %>%
        extract(
          ROIcat_revpairwise,
          c("ROIcat_1", "ROIcat_2"),
          "(.*) / (.*)",
          remove = FALSE
        ) %>%
        mutate(
          across(
            c(ROIcat_1, ROIcat_2),
            ~ factor(.x, levels = intersect(lesion_type_order, unique(.x)))
          )
        )
    ),
    contrasts_by_stages = map(
      model,
      \(x) emmeans(
        x,
        specs = ~ marker_status_1 | ROIcat
      ) %>%
        contrast(method = "revpairwise", adjust = "bh", type = "response") %>% {
          inner_join(
            as_tibble(.),
            as_tibble(confint(.)) %>%
              select(contrast, ROIcat, starts_with("asymp")),
            by = c("contrast", "ROIcat")
          )
        }
    ),
    across(
      starts_with("contrast"),
      \(x)
      map(
        x,
        \(y) mutate(
          y,
          p_cut = cut(p.value, c(-Inf, .001, .01, .05, Inf), c("p<0.001", "p<0.01", "p<0.05", "ns"), ordered_result = TRUE),
          p_stars = cut(p.value, c(-Inf, .001, .01, .05, Inf), c("***", "**", "*", "ns"), ordered_result = TRUE)
        )
      )
    )
  )

cell_populations_across_stages <- cell_populations_models %>%
  select(
    marker_1, marker_2, Incidental, contrasts_across_stages
  ) %>%
  unnest(contrasts_across_stages)

cell_populations_by_stages <- cell_populations_models %>%
  select(
    marker_1, marker_2, Incidental, contrasts_by_stages
  ) %>%
  unnest(contrasts_by_stages)

cell_populations_models_plots <- cell_populations_models %>%
  mutate(
    lolipop_plot = pmap(
      list(contrasts_across_stages, contrasts_by_stages, marker_1, marker_2),
      \(contrast_across_stages, contrast_by_stages, marker_1, marker_2) {
        ggplot(
          contrast_by_stages,
          aes(
            x = ROIcat,
            y = odds.ratio,
            color = p_cut
          )
        ) +
          geom_hline(
            yintercept = 1,
            linetype = "dotted",
            color = "gray"
          ) +
          geom_segment(
            aes(
              xend = ROIcat,
              yend = 1
            ),
            linewidth = .5,
            color = "gray40"
          ) +
          geom_errorbar(
            aes(
              ymin = asymp.LCL,
              ymax = asymp.UCL
            ),
            color = "black",
            width = 0.2
          ) +
          geom_point(
            shape = 16,
            size = 4
          ) + {
            signif_across_stages <- contrast_across_stages %>%
              filter(p.value < 0.05)
            if (nrow(signif_across_stages) > 0) {
              ggpubr::geom_bracket(
                aes(
                  xmin = ROIcat_1,
                  xmax = ROIcat_2,
                  label = p_stars
                ),
                data = signif_across_stages,
                y.position = max(contrast_by_stages$asymp.UCL) * 1.1,
                step.increase = 0.1,
                inherit.aes = FALSE
              )
            }
          } +
          scale_y_log10() +
          scale_color_viridis_d(
            drop = FALSE,
            direction = -1
          ) +
          labs(
            title = paste(marker_1, "vs", marker_2),
            y = "Fewer double positive ⟵   Odds ratio   ⟶ More double positive",
            x = NULL
          )
      }
    )
  )
pmap(
  cell_populations_models_plots, 
  \(marker_1,marker_2,Incidental,lolipop_plot, ...){
    ggsave(paste0("lolipop_plots_", Incidental, "_", marker_1,"vs",marker_2,".pdf"),lolipop_plot, width=4, height=3)
  }
)


# To use lolipop script One has to plot the OR with CI:
#change the ROICat to .I or .C manually, saved just the ORs to make a dataframe and then plot it here:
OR_df <- read_csv("glm_OR_pTBK1_HLA_Ep_STAT3_pSTAT1_input_lolipop_v2.csv")
OR_df$Cat_collapse<- factor(OR_df$Cat_collapse, levels = c("FT.I", "Fim.I", "p53.I","STIC.I", "STIC.C", "Cancer"))
ggplot(OR_df, aes(x=`Cat_collapse`, y=OR, color= `marker`))+
  geom_segment( aes(x=`Cat_collapse`, xend=`Cat_collapse`, y=1, yend=OR), 
                linewidth = .2,
                color = "white"
                #color=ifelse(OR_df$`Cat_collapse` %in% c("STIC.C","Cancer"), "orange", "cyan"),
                #size=ifelse(OR_df$`Cat_collapse` %in% c("STIC.C","Cancer"), 2, 2)
  )+ 
  
  geom_point(
    shape = 16,
    size = 4,
    position = position_dodge(1)
  )+
  #geom_point(
    #color=ifelse(OR_df$`Cat_collapse` %in% c("STIC.C","Cancer"), "orange", "cyan"), 
    #size=ifelse(OR_df$`Cat_collapse` %in% c("STIC.C","Cancer"), 3, 3)
  #)+
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    color = "gray"
  ) +
  geom_errorbar(
    aes(
      ymin = LCL,
      ymax = UCL,
      color= `marker`
      ),
    #color = "black",
    width = 0.2,
    position = position_dodge(1)
  ) +
  theme_classic()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Lesion Types") +
  ylab("OR")+
  scale_y_log10() +
  scale_color_viridis_d(
    drop = FALSE,
    direction = -1
  )
 
  

ggsave(
  "pTBK1_HLA-E_pSTAT3_pSTAT1_correlation_OR.pdf",
  width = 10, height =5
)



