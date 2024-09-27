library(tidyverse)
library(readxl)
library(seriation)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(hrbrthemes)
library(viridis)
library(powerjoin)
library(broom)
library(broom.mixed)
library(emmeans)

marker_data <- read_csv("panckposmean.csv")

marker_data_long <- marker_data %>%
  select(
    slideName, ROIid, GroupCount,
    BRCA, Stage, Incidental, Cat,
    matches("^mean_[A-Z0-9_]+[pn]$")
  ) %>%
  pivot_longer(
    cols = matches("^mean_[A-Z0-9_]+[pn]$"),
    names_to = "population",
    values_to = "proportion"
  ) %>%
  filter(!Cat %in% c("Inc.", "Non")) %>%
  mutate(
    population = str_replace(population, fixed("mean_"), ""),
    n = round(GroupCount * proportion)
  )

marker_data_long_subset <- marker_data_long %>%
  mutate(
    Cat_collapse=fct_relabel(Cat, \(x) str_replace(x, fixed(".(Inc)"), ".I")) %>%
      fct_relabel(\(x) str_replace(x, fixed(".(Non)"), ".C"))
  )

marker_data_long_subset2 <- marker_data_long_subset %>%
  filter(!Cat_collapse %in% c("", "Fim.C", "FT.C", "Inc.", "Non", "p53.C", "NA"))
marker_data_long_subset2$Cat_collapse<- factor(marker_data_long_subset2$Cat_collapse, levels = c("FT.I", "Fim.I", "p53.I","STIC.I", "STIC.C", "Cancer"))


markers_of_interest <- c(
  "CD4pPD1p", "CD4pLAG3p", "CD4pFoxP3p", "CD4pHLADRp", "CD4pPD1pLAG3p"
)
marker_data_long_subset3 <- marker_data_long_subset2 %>%
  filter(population %in% markers_of_interest)
marker_data_long_cd4p <- marker_data_long_subset2 %>%
  filter(population %in% c("CD4p")) %>%
  group_by(slideName, ROIid, Cat_collapse) %>%
  summarize(proportion = sum(proportion), .groups = "drop")%>%
  group_by(Cat_collapse) %>%
  summarize(mean_proportion = mean(proportion), .groups = "drop")
mean_marker_data_long_subset3 <- marker_data_long_subset3 %>%
  group_by(Cat_collapse, population) %>%
  summarise(mean_proportion = mean(proportion)) %>%
  left_join(
    marker_data_long_cd4p %>%
      select(Cat_collapse, mean_proportion_CD4p = mean_proportion)
  ) %>%
  mutate(
    proportion_of_CD4p = mean_proportion / mean_proportion_CD4p
  )

mean_marker_data_long_subset3$population<- factor(mean_marker_data_long_subset3$population, levels = c("CD4pHLADRp", "CD4pPD1p", "CD4pLAG3p","CD4pFoxP3p", "CD4pPD1pLAG3p"))


#summarise the ratio of each subset of cell types 
ggplot(data = mean_marker_data_long_subset3, aes(x = `Cat_collapse`, y = `proportion_of_CD4p`, Cat = factor(`Cat_collapse`), fill = factor(`population`))) + 
  geom_col(width=0.5) + 
  ylab("Ratio: subset of CD4+ T cells/CD4+ T cell (%)") + 
  scale_y_continuous(labels=scales::percent, expand = c(0, 0)) +
  xlab("Type of Lesions") +
  labs(fill = "Population") +
  theme_classic()+ 
  ggokabeito::scale_fill_okabe_ito() +
  theme(axis.text = element_text(color="black"), axis.title =element_text(color="black"))+
  theme(legend.position = c(0.2, 0.8))


ggsave(
  "ratio_CD4+ T cell in total CD4T_epithelia.pdf",
  width = 5, height = 4
)

#LMM stat:
marker_data_long_cd4p_by_roi <- marker_data_long_subset2 %>%
  filter(population %in% c("CD4p")) %>%
  group_by(slideName, ROIid, Cat_collapse) %>%
  summarize(proportion_CD4p = sum(proportion), .groups = "drop")

marker_data_long_subset3_with_CD4p <- marker_data_long_subset3 %>%
  power_inner_join(
    marker_data_long_cd4p_by_roi %>%
      select(slideName, ROIid, proportion_CD4p),
    by = c("slideName", "ROIid"),
    check = check_specs(
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  ) %>%
  mutate(
    n_CD4p = round(proportion_CD4p * GroupCount)
  )

glm_of_interest_random <- marker_data_long_subset3_with_CD4p %>%
  filter(population %in% markers_of_interest) %>%
  group_nest(population) %>%
  mutate(
    mod = map(
      data,
      \(x) lme4::glmer(
        cbind(n, n_CD4p - n) ~ Cat + (1 | slideName),
        data = x,
        family = binomial
      )
    ),
    mod_glance = map(mod, glance)
  )

glm_of_interest_random_res <- glm_of_interest_random %>%
  mutate(
    mod_emmeans = map(
      mod,
      \(x) emmeans(
        x,
        ~ Cat,
        type = "response"
      )
    ),
    mod_contrast = map(
      mod_emmeans,
      \(x) contrast(x, method = "revpairwise", adjust = "fdr")
    ),
    mod_contrast_ci = map(
      mod_contrast, confint
    )
  )

glm_of_interest_random_tidy <- glm_of_interest_random_res %>%
  mutate(
    mod_contrast_tidy = map2(
      mod_contrast, mod_contrast_ci,
      \(x, y) tidy(x) %>%
        power_inner_join(
          y %>%
            select(contrast, odds.ratio.lower.ci = asymp.LCL, odds.ratio.upper.ci = asymp.UCL),
          by = "contrast",
          check = check_specs(
            unmatched_keys_left = "warn",
            unmatched_keys_right = "warn"
          )
        )
    )
  ) %>%
  select(population, mod_contrast_tidy) %>%
  unnest(mod_contrast_tidy)

#dir.create("cell_type_proportions", showWarnings = FALSE)
write_csv(
  glm_of_interest_random_tidy,
  "CD4_subtypesratio_panckpos_binomial_glm_random_p_values.csv"
)
#here is the panckneg population

marker_data <- read_csv("pancknegmean.csv")
# replace the panck pos popualtions with negative and run the script again. 