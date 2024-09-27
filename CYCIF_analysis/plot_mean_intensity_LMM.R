library(tidyverse)
library(seriation)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(powerjoin)
library(tidyr)
library(plotly)
library(ggfortify)
library(GGally)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(vioplot)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(brms)
library(emmeans)
library(tidybayes)
library(broom)
library(broom.mixed)



#LMMmodeliing for stat test

marker_data <- read_csv("panckposmean.csv")

marker_data_long <- marker_data %>%
  select(
    slideName, ROIid, GroupCount,
    BRCA, Stage, Incidental, Cat,
    matches("^mean_[A-Z0-9_]")
  ) %>%
  pivot_longer(
    cols = matches("^mean_[A-Z0-9_]"),
    names_to = "marker",
    values_to = "intensity"
  ) %>%
  filter(!Cat %in% c("Inc.", "Non")) %>%
  mutate(
    marker = str_replace(marker, fixed("mean_"), "")
    )

markers_of_interest <- c(
  "p53")


lm_of_interest_random <- marker_data_long %>%
  filter(marker %in% markers_of_interest) %>%
  group_nest(marker) %>%
  mutate(
    mod = map(
      data,
      \(x) lme4::lmer(
        log10(intensity) ~ Cat + (1 | slideName),
        data = x
      )
    ),
    mod_glance = map(mod, glance)
  )

lm_of_interest_random_res <- lm_of_interest_random %>%
  mutate(
    mod_emmeans = map(
      mod,
      \(x) emmeans(
        x,
        ~ Cat,
        type = "response"
      ) %>%
        contrast(method = "revpairwise", adjust = "fdr") %>%
        tidy()
    )
  ) %>%
  select(marker, mod_emmeans) %>%
  unnest(mod_emmeans)

#dir.create("cell_type_proportions", showWarnings = FALSE)
write_csv(
  lm_of_interest_random_res,
  "intensity_p53_lm_random_p_values.csv"
)








