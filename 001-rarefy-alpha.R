# need to run independently. The order of OTU is corrected for joining the metadata
library(tidyverse)
library(vegan) # tools for descriptive community ecology. Bray-Curtis and Jaccard analyses
library(readxl)
library(ggplot2)

# read rarefied ASV data
rarefy_V4 <- read.csv(here::here("analysis","data","derived_data",
                                 "V4_shared_format.1.subsample.60216.csv"), sep = ",")
rarefy_V1 <- read.csv(here::here("analysis","data","derived_data",
                                 "V1_shared_format.1.subsample.13658.csv"), sep = "")
rarefy_V6 <- read.csv(here::here("analysis","data","derived_data",
                                 "V6_shared_format.1.subsample.28402.csv"), sep = ",")
rarefy_ITS1 <- read.csv(here::here("analysis","data","derived_data",
                                   "ITS1_shared_format.1.subsample.387.csv"), sep = ",")
context <- read_excel(here::here("analysis","data", "raw_data", "GJB_amp_meta.xlsx"))

# to get the same order for all
rarefy_V4 <- rarefy_V4 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V1 <- rarefy_V1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V6 <- rarefy_V6 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_ITS1 <- rarefy_ITS1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)

# tidy meta data and get the same order as the rarefy data
context_tidy <-
  context %>%
  mutate(`adpater no` = tolower(`adpater no`),
         type = tolower(str_extract(`adpater no`, "[[:alpha:]]+(?=\\d)")),
         type = case_when(type == "con" ~ "control", TRUE ~ as.character(type)),
         area = case_when(area == "inside" ~ "interior", area == "outside" ~ "exterior"),
         group = if_else(!is.na(area)&!area%in% c("belly", "outside belly"),
                         paste0(type, "-", area), type),
         group = factor(group, levels = c("control","pot-interior","pot-exterior","soil"))) %>%
  arrange(`adpater no`)

# function to have OTU reads only
reads_only <- function(df) {
  df %>%
    select(-label, -Group, -numOtus)
}

rarefy_V4_reads_only  <- reads_only(rarefy_V4)
rarefy_V1_reads_only  <- reads_only(rarefy_V1)
rarefy_V6_reads_only  <- reads_only(rarefy_V6)
rarefy_ITS1_reads_only  <- reads_only(rarefy_ITS1)

# function to get context for each region
meta_region <- function(df, region){
  context_tidy %>%
    filter(Region == region) %>%
    arrange(`adpater no`)
}

V4_meta<- meta_region(context_tidy, "V4")
V1_meta<- meta_region(context_tidy, "V1")
V6_meta<- meta_region(context_tidy, "V6")
ITS1_meta<- meta_region(context_tidy, "ITS1")

# statistical testing for the richness, following https://rpubs.com/an-bui/vegan-cheat-sheet
richness_test <- function (df1, df2){
  rich <- specnumber(df1) # the number of species within each sample
  rich_aov <- aov(rich ~ group, data =df2) # analysis of variance between groups
  return(c(rich, summary(rich), summary(rich_aov)))
}

rich_V4_test  <- richness_test(rarefy_V4_reads_only, V4_meta) # rich_V4[27][[1]]$`Pr(>F)`[1]
rich_V1_test  <- richness_test(rarefy_V1_reads_only, V1_meta)
rich_V6_test  <- richness_test(rarefy_V6_reads_only, V6_meta)
rich_ITS1_test  <- richness_test(rarefy_ITS1_reads_only, ITS1_meta)

pair_richness_test <- function (df1, df2){
  rich <- specnumber(df1) # the number of species within each sample
  pair_rich_aov <- pairwise.t.test(rich ~ group, data =df2, p.adj = "none") # analysis of variance between groups
  return(c(rich, summary(rich), summary(rich_aov)))
}

sppdiv_aov_V4 <- aov(diversity(rarefy_V4_reads_only) ~ group, data =V4_meta)
sppdiv_aov_V1 <- aov(diversity(rarefy_V1_reads_only) ~ group, data =V1_meta)
sppdiv_aov_V6 <- aov(diversity(rarefy_V6_reads_only) ~ group, data =V6_meta)
sppdiv_aov_ITS1 <- aov(diversity(rarefy_ITS1_reads_only) ~ group, data =ITS1_meta)

# functions calculate richness, Chao1, evenness index, and Shannon index using vegan package
# follow https://scienceparkstudygroup.github.io/microbiome-lesson/04-alpha-diversity/index.html
alpha_cal_df <- function(df1, df2){
  richness <- estimateR(df1)
  evenness <- diversity(df1)/log(specnumber(df1))
  shannon <- diversity(df1, index = "shannon")
  alphadiv <- cbind(df2, t(richness), shannon, evenness) #combine all in one
  #rm(richness, evenness, shannon) #remove the unnecessary data/vector
  return(alphadiv)
}

V4_alphadiv <- alpha_cal_df(rarefy_V4_reads_only, V4_meta)
V1_alphadiv <- alpha_cal_df(rarefy_V1_reads_only, V1_meta)
V6_alphadiv <- alpha_cal_df(rarefy_V6_reads_only, V6_meta)
ITS1_alphadiv <- alpha_cal_df(rarefy_ITS1_reads_only, ITS1_meta)

# pairwise t-test on richness and shannon
library(rstatix)
rich_pair_V4 <- pairwise_t_test(V4_alphadiv, S.obs ~ group)
rich_pair_V1 <- pairwise_t_test(V1_alphadiv, S.obs ~ group)
rich_pair_V6 <- pairwise_t_test(V6_alphadiv, S.obs ~ group)
rich_pair_ITS1 <- pairwise_t_test(ITS1_alphadiv, S.obs ~ group)

options(scipen = 10, digits = 2)

# function to tidy data for alpha analysis, including richness, chao, evenness, and shannon
alphadiv_tidy <- function(df){
  df %>%
    select(`adpater no`, type, group, area, S.obs, S.chao1, evenness, shannon) %>%
    pivot_longer(- c(`adpater no`, type, group, area),
                 names_to = "variable", values_to = "value") %>%
    mutate(variable = case_when(variable == "S.obs" ~"Richness", variable == "S.chao1" ~"Chao1",
                                variable == "evenness" ~"Evenness", variable == "shannon" ~"Shannon")) %>%
    mutate(variable = factor(variable, levels = c("Richness","Evenness","Shannon","Chao1")))
}

V4_alphadiv_tidy <- alphadiv_tidy(V4_alphadiv)
V1_alphadiv_tidy <- alphadiv_tidy(V1_alphadiv)
V6_alphadiv_tidy <- alphadiv_tidy(V6_alphadiv)
ITS1_alphadiv_tidy <- alphadiv_tidy(ITS1_alphadiv)

# function to make boxplot
alphadiv_plot <-function(df, vec){
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.9) +
    geom_point(color = "grey30", size = 1) +
    #scale_x_discrete(labels = c("control\n(n = 3)", "pot-interior\n(n = 6)",
                               # "pot-exterior\n(n = 6)","soil\n(n = 5)")) +
    labs(x = NULL,
         y = NULL,
         title = vec) +
    scale_fill_viridis_d(name = "Sample group", #scale_fill_discrete
                         labels = c("control (n = 3)", "pot-interior (n = 6)", "pot-exterior (n = 6)", "soil (n = 5)")) +
    facet_wrap(~variable, scales = "free") +
    theme_minimal() +
    #theme(legend.position = "none") +
    theme(plot.background = element_rect(fill = "white", colour = "white"),
          strip.text = element_text(size= 12))
}

V4_alphadiv_plot <- alphadiv_plot(V4_alphadiv_tidy, "V4")
V1_alphadiv_plot <- alphadiv_plot(V1_alphadiv_tidy, "V1")
V6_alphadiv_plot <- alphadiv_plot(V6_alphadiv_tidy, "V6")
ITS1_alphadiv_plot <- alphadiv_plot(ITS1_alphadiv_tidy, "ITS1")

# save plots by regions
ggsave(here::here("analysis", "figures", "01-Alpha_ITS1_plots.png"),
       width = 8, height = 5, units = "in")
