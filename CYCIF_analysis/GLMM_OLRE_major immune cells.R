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
library(lme4)
library(glmmTMB)
library(DHARMa)
library(rlang)



#read the files we need to plot from: 

panCKpos <- read_csv("panckposmean.csv")
panCKneg <- read_csv("pancknegmean.csv")

#Then change the name of the categories the way we want to plot:

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

#Binomial GLMMmodeliing with patient & observational random effect (OLRE) for stat test
#PanCKpos file: 
model_funs <- list(
  # nb_g = \(x) glmer.nb(
  #   n ~ Cat + offset(log(GroupCount)) + (1 | slideName),
  #   data = x,
  #   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e6))
  # ),
  #nb_1 = \(x) glmmTMB(
  #n ~ Cat + offset(log(GroupCount)) + (1 | slideName),
  #data = x,
  #family = nbinom1(link = "log")
  #),
  #nb_2 = \(x) glmmTMB(
  #n ~ Cat + offset(log(GroupCount)) + (1 | slideName),
  #data = x,
  #family = nbinom2(link = "log")
  #),
  #nb_2_zi = \(x) glmmTMB(
  #n ~ Cat + offset(log(GroupCount)) + (1 | slideName),
  #data = x,
  #ziformula = ~ Cat + (1 | slideName),
  #family = nbinom2(link = "log")
  #),
  #bb = \(x) glmmTMB(
  #cbind(n, GroupCount - n) ~ Cat + (1 | slideName),
  #data = x,
  #family = betabinomial(link = "logit")
  #),
  #b = \(x) glmmTMB(
  #cbind(n, GroupCount - n) ~ Cat + (1 | slideName),
  #data = x,
  #family = binomial(link = "logit")
  #),
  b_olre = \(x) glmmTMB(
    cbind(n, GroupCount - n) ~ Cat + (1 | slideName) + (1 | measurement_id),
    data = x,
    family = binomial(link = "logit")
  )#,
  #b_olre_wo_patient = \(x) glmmTMB(
  #cbind(n, GroupCount - n) ~ Cat + (1 | measurement_id),
  #data = x,
  #family = binomial(link = "logit")
  #)
)

marker_data <- read_csv("panckposmean.csv")

stage_levels <- c(
  "FT.(Inc)", "FT.(Non)", "Fim.(Inc)", "Fim.(Non)", "p53.(Inc)",
  "p53.(Non)", "STIC.(Inc)", "STIC.(Non)", "Cancer"
)

abbreviate_level <- function(x) {
  str_replace(x, fixed(".(Inc)"), ".I") %>%
    str_replace(fixed(".(Non)"), ".C")
}

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
    n = round(GroupCount * proportion),
    Cat = factor(
      Cat,
      levels = stage_levels
    ),
    measurement_id = paste0("measurement_", seq_len(n()))
  )

markers_of_interest <- c(
  "CD45pCD8ap", "CD45pCD4p", "CD20p",
  "CD8apCD103n", "CD8apCD103p",  "CD68pCD11cp", "CD11cpCD103pCD68n", "CD163p", "CD68p", "HLADRp", "CD11cp"
)



#Once we have markers of interest, then run this drop

  glm_of_interest_random <- marker_data_long %>%
    filter(population %in% markers_of_interest) %>%
    group_nest(population) %>%
    crossing(model_name = names(model_funs)) %>%
    mutate(
      mod_list = pmap(
        list(data, population, model_name),
        \(x, y, z) {
          message("Processing ", y, " with ", z)
          quietly(model_funs[[z]])(x)
        }
      ),
      mod = map(mod_list, "result"),
      mod_warnings = map(mod_list, "warnings"),
      mod_glance = map(mod, possibly(glance)),
      mod_emmeans = map(
        mod,
        \(x) possibly(emmeans)(
          x,
          ~ Cat,
          type = "response"
        )
      ),
      mod_contrasts = map(
        mod_emmeans,
        \(x) {
          set_names(c("trt.vs.ctrl1", "consec")) %>%
            map(\(y) possibly(contrast)(x, method = y, adjust = "none"))
        }
      ),
      mod_residuals = map(mod, possibly(simulateResiduals))
    )
  
  #library(qs)
  #qsave(
    #glm_of_interest_random,
    #"glm_of_interest_random_epithelia.qs"
  #)
  
  glm_of_interest_random %>%
    filter(!map_lgl(mod_warnings, is_empty)) %>%
    count(model_name)
  # None of the b_olre models have warnings!
  
  glm_of_interest_random_glances <- glm_of_interest_random %>%
    select(population, model_name, mod_glance) %>%
    unnest(mod_glance)
  
  glm_of_interest_random_res <- glm_of_interest_random %>%
    select(
      population, model_name, mod_contrasts
    ) %>%
    unnest_longer(mod_contrasts, indices_to = "contrast_type") %>%
    mutate(across(mod_contrasts, \(x) map(x, tidy))) %>%
    unnest(mod_contrasts) %>%
    group_by(population, contrast, model_name) %>%
    # Some comparisons made twice, so take the first one
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(population, model_name) %>%
    mutate(
      padj = p.adjust(p.value, method = "fdr"),
      padj_cut = cut(padj, breaks = c(-Inf, .01, 0.05, Inf), labels = c("**", "*", "ns"))
    ) %>%
    ungroup() %>%
    separate_wider_delim(
      contrast,
      " / ",
      names = c("Cat_1", "Cat_2"),
      cols_remove = FALSE
    ) %>%
    mutate(
      across(
        c(Cat_1, Cat_2),
        \(x) factor(x, levels = stage_levels)
      )
    ) %>%
    select(-term, -starts_with("null"), -df)
  
  write_csv(
    glm_of_interest_random_res,
    "cell_type_proportions_epithelial_all_glm_res.csv"
  )
  
  
  #stat summary for median, range etc
  panCKpos_collapse_subset %>%
    select(Cat_collapse, mean_HLA_Ep, mean_HLA_Ap, mean_Ki67p, mean_gH2Axp, 
           mean_CD45pCD8ap, mean_CD20p, mean_CD45pCD4p, mean_CD8apCD103n,
           mean_CD8apCD103p, mean_CD8apCD103nKi67p, mean_CD8apCD103nLAG3p, mean_CD8apCD103nPD1p,
           mean_CD8apCD103pKi67p, mean_CD8apCD103pLAG3p, mean_CD8apCD103pPD1p,
           mean_CD8apGZMBp, mean_CD4pFoxP3p, mean_CD4pHLADRp, mean_CD4pLAG3p,
           mean_CD4pPD1p, mean_CD68pHLADRp, mean_CD68pCD11cp, mean_CD11cpCD103pCD68n, mean_CD11cp, 
           mean_CD11cpHLADRpCD68n, mean_CD163p, mean_CD68p, mean_CD8apLAG3p, mean_CD8apPD1p, mean_CD8apKi67p, mean_pSTAT3p, mean_pTBK1p) %>%
    split(.$Cat_collapse) %>%
    map(summary)

  #panCKneg: replace panckpos file with negative and run the whole script as above.