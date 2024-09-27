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


#read the files we need to plot from: 

panCKpos <- read_csv("panckposmean.csv")
panCKneg <- read_csv("pancknegmean.csv")

#Then change the name of the categories

panCKposcollapse <- panCKpos %>%
  mutate(
    Cat_collapse=fct_relabel(Cat, \(x) str_replace(x, fixed(".(Inc)"), ".I")) %>%
      fct_relabel(\(x) str_replace(x, fixed(".(Non)"), ".C"))
  )

panCKnegcollapse <- panCKneg %>%
  mutate(
    Cat_collapse=fct_relabel(Cat, \(x) str_replace(x, fixed(".(Inc)"), ".I")) %>%
      fct_relabel(\(x) str_replace(x, fixed(".(Non)"), ".C"))
  )



#This is to check the variable under each colummn: 
#panCKposcollapse$Cat_collapse %>% table()


#All panckpos markers first: 


#subsetting the data basically to remove NA/blank from the graph:
#filter to remove the datapoint we dont want to plot
panCKpos_collapse_subset <- panCKposcollapse %>%
  filter(!Cat_collapse %in% c("", "Fim.C", "FT.C", "Inc.", "Non", "p53.C", "NA"))
panCKpos_collapse_subset$Cat_collapse<- factor(panCKpos_collapse_subset$Cat_collapse, levels = c("FT.I", "Fim.I", "p53.I","STIC.I", "STIC.C", "Cancer"))



#This is just for few epithelial marker without the log scale, 
#HLA-E, Ki67, HLA-A, gamma H2aX, p-STAT3+

plot <- ggplot(panCKpos_collapse_subset, aes(x = `Cat_collapse`, y = `mean_Ki67p`, Cat = factor(`Cat_collapse`), color=`Cat_collapse`))+ 
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_quasirandom(width = 0.0001) +
  #coord_cartesian(ylim=c(0,0.55))+
  scale_y_continuous(labels=scales::percent, expand = c(0, 0))+
  ylab(" cells in epithelia (%)") + 
  xlab("Type of Lesions") +
  theme_classic()+ 
  ggokabeito::scale_fill_okabe_ito() +
  theme(axis.text = element_text(color="black"), axis.title =element_text(color="black"))+
  theme(legend.position = "none")+
  ggtitle ("Ki67+ epithelial cells (%)")
plot 
ggsave(
  "cells with Ki67+ in epithelia.pdf",
  width = 4.5, height = 4
)


#immune cells populations are small, so we needed to do logscale

#writing functions to convert the yscale for log transformed

replace_zeroes <- function(df, col) {
  val <- zero_location(df, {{col}})
  mutate(df, {{col}} := if_else({{col}} == 0, val, {{col}}))
}

zero_location <- function(df, col) {
  val <- pull(df, {{col}}) %>%
    .[. > 0] %>%
    min()
  val / 2
}

#plot with boxplot+ggbeeswarm: PanCKpos: small proportion of cells log transformed
plot <- replace_zeroes(panCKpos_collapse_subset, mean_CD68p) %>%
  ggplot(aes(x = `Cat_collapse`, y = `mean_CD68p`, Cat = factor(`Cat_collapse`), fill= `Cat_collapse`, color = 'Cat_collapse'))+  
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_quasirandom(width = 0.0001) +
  geom_hline(yintercept = zero_location(panCKpos_collapse_subset, mean_CD68p), linetype = "dashed") +
  scale_y_continuous(labels=scales::percent, trans="log10") +
  ylab("CD68+ in epithelia (%)") + 
  xlab("Type of Lesions") +
  theme_classic()+ 
  ggokabeito::scale_fill_okabe_ito() +
  theme(axis.text = element_text(color="black"), axis.title =element_text(color="black"))+
  theme(legend.position = "none")+
  ggtitle ("CD68+ cells in epithelia")
plot 
ggsave(
  "CD68+ in epithelial_log scale.pdf",
  width = 4.5, height = 4
)

#LMMmodeliing for stat test

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
markers_of_interest <- c(
  "HLA_Ep", "HLA_Ap", 
  "Ki67p", "gH2Axp",
  "pSTAT3p", "pTBK1p"
)
glm_of_interest_random <- marker_data_long %>%
  filter(population %in% markers_of_interest) %>%
  group_nest(population) %>%
  mutate(
    mod = map(
      data,
      \(x) lme4::glmer(
        cbind(n, GroupCount - n) ~ Cat + (1 | slideName),
        data = x,
        family = binomial
      )
    ),
    mod_glance = map(mod, glance),
    mod_emmeans = map(
      mod,
      \(x) emmeans(
        x,
        ~ Cat,
        type = "response"
      ) %>%
        contrast(method = "revpairwise", adjust = "fdr")
    )
  )

glm_of_interest_random_res <- glm_of_interest_random %>%
  transmute(population, mod_emmeans = map(mod_emmeans, tidy)) %>%
  unnest(mod_emmeans)

#dir.create("cell_type_proportions", showWarnings = FALSE)
write_csv(
  glm_of_interest_random_res,
  "cell_type_epithelial_proportions_binomial_glm_random_pvalue.csv"
)


#panCKneg marker

#replace the panckpos file name with panCK neg file name and run the script. 

