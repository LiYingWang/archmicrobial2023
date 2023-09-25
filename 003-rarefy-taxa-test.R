# can run independently
# follow https://riffomonas.org/code_club/2021-07-02-wilcox.test
library(tidyverse)
library(broom)
library(Hmisc)
library(readxl)

# function to read data
shared_rarefy_long <- function(filename, vec){
  read.csv(here::here("analysis","data","derived_data", filename), sep = vec) %>%
    select(Group, starts_with("ASV")) %>%
    pivot_longer(-Group, names_to="otu", values_to="count")
}

shared_V4_rarefy <-shared_rarefy_long("V4_shared_format.1.subsample.60216.csv", ",")
shared_V1_rarefy <-shared_rarefy_long("V1_shared_format.1.subsample.13658.csv", "")
shared_V6_rarefy <-shared_rarefy_long("V6_shared_format.1.subsample.28402.csv", ",")
shared_ITS1_rarefy <-shared_rarefy_long("ITS1_shared_format.1.subsample.387.csv", ",")

# function to tidy taxonomy data
tidy_taxonomy <- function(filename){
  read.table(here::here("analysis","data","raw_data",filename), header = T, sep="\t") %>%
  select(OTU.ID, taxonomy) %>%
  rename(otu = OTU.ID) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, "; c__; o__; f__; g__; s__", ""),
         taxonomy = str_replace(taxonomy, "; o__; f__; g__; s__", ""),
         taxonomy = str_replace(taxonomy, "; f__; g__; s__", ""),
         taxonomy = str_replace(taxonomy, "; g__; s__", ""),
         taxonomy = str_replace(taxonomy, "(; s__)([^(; s__)]*)$", ""),
         taxonomy = str_replace(taxonomy, "(; f__)([^(; f__)]*)$", ""),
         taxonomy = str_replace(taxonomy, "(; g__)([^(; g__)]*)$", ""))
  }

taxonomy_V4 <- tidy_taxonomy("V4_asv_table.txt")
taxonomy_V1 <- tidy_taxonomy("V1_asv_table.txt")
taxonomy_V6 <- tidy_taxonomy("V6_asv_table.txt")
taxonomy_ITS1 <- tidy_taxonomy("ITS1_asv_table.txt")

# function to read metadata for each region
meta_region_extract <- function(filename, vec){
  read_excel(here::here("analysis","data", "raw_data", filename)) %>%
    filter(Region == vec) %>%
    rename(Group = `adpater no`) %>%
    mutate(type= str_extract(Group, "[[:alpha:]]+(?=\\d)"), type = tolower(type)) %>%
    mutate(group = if_else(!is.na(area)&!area%in% c("belly", "outside belly"),
                           paste0(type, "-", area), type),
           pot = case_when(group == "pot-inside" ~ "interior", group == "pot-outside" ~"exterior"),
           potex_soil = case_when(group == "pot-outside" ~ "pot-exterior", type== "soil" ~ "soil"),
           potin_soil = case_when(group == "pot-inside" ~ "pot-interior", type== "soil" ~ "soil"),
           pot_con = case_when(type== "pot" ~ "pot", type== "con" ~ "control"),
           potex_con = case_when(group == "pot-outside" ~ "pot-exterior", type== "con" ~ "control"),
           potin_con = case_when(group == "pot-inside" ~ "pot-interior", type== "con" ~ "control")) %>%
    relocate(Group, type, group, sample, area, site)
}

metadata_V4 <- meta_region_extract("GJB_amp_meta.xlsx", "V4")
metadata_V1 <- meta_region_extract("GJB_amp_meta.xlsx", "V1")
metadata_V6 <- meta_region_extract("GJB_amp_meta.xlsx", "V6")
metadata_ITS1 <- meta_region_extract("GJB_amp_meta.xlsx", "ITS1")

# function to calculate relative abundance for taxa level, not ASV level
composite_region <- function(df1, df2, df3){
  inner_join(df1, df2, by= "otu") %>%
    group_by(Group, taxonomy) %>%
    dplyr::summarize(count = sum(count), .groups="drop") %>%
    group_by(Group) %>%
    mutate(rel_abund = count / sum(count)) %>%
    ungroup() %>%
    select(-count) %>%
    inner_join(., df3, by = "Group") %>%
    filter(!taxonomy == "k__Bacteria")
}

composite_V4 <- composite_region(shared_V4_rarefy, taxonomy_V4, metadata_V4)
composite_V1 <- composite_region(shared_V1_rarefy, taxonomy_V1, metadata_V1)
composite_V6 <- composite_region(shared_V6_rarefy, taxonomy_V6, metadata_V6)
composite_ITS1 <- composite_region(shared_ITS1_rarefy, taxonomy_ITS1, metadata_ITS1)

# get significant genera
sig_genera_V4 <-
  composite_V4 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~  pot_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.1)
  #mutate(p.adjust = p.adjust(p.value, method="BH")) %>% #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
  #filter(p.adjust < 0.01) %>%
  #select(taxonomy, p.adjust)

sig_genera_V1 <-
  composite_V1 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~  pot_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.1)

sig_genera_V6 <-
  composite_V6 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~  pot_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.1)

sig_genera_ITS1 <-
  composite_ITS1 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~  pot_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.1)

# make a plot
library(ggtext)
composite_V4_plot <-
  composite_V4 %>%
  inner_join(sig_genera_V4, by="taxonomy") %>%
  filter(!is.na(pot)) %>%
  mutate(taxa = ifelse(str_detect(taxonomy, "g__"), str_extract(taxonomy, "(g__)([A-z]+)", group =2),
                       str_extract(taxonomy, "(f__)([A-z]+)", group =2)), #"g__[A-z]+
         taxa= str_replace(taxa, "(.*)", "*\\1*"),
         taxa = str_replace(taxa, "\\*Candidatus_(.*)\\*", # "\\*(.*)_Candidatus\\*"
                                "Candidatus<br>*\\1*")) %>%
  mutate(rel_abund = 100 * (rel_abund + 1/15000), # make it percentage and add a very small number
         pot = factor(pot, levels = c("exterior", "interior"))) %>%
  ggplot(aes(x=rel_abund, y=reorder(taxa, rel_abund), fill=pot)) +
  geom_boxplot(aes(fill = pot),show.legend = FALSE) +
  #geom_vline(xintercept = 100/8490, size=0.5, color="gray") + # limit of detection; https://www.nature.com/articles/s41467-023-38694-0
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.5), shape=21, size =2) +
  #stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               #geom="pointrange",
               #position = position_dodge(width=0.8),
               #color="black", show.legend = FALSE) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("interior", "exterior"),
                     values = c("#238A8DFF", "ivory3"),
                     labels = c("Interior", "Exterior")) +
  scale_fill_manual(NULL,
                    breaks = c("interior", "exterior"),
                    values = c("#238A8DFF", "ivory3"),
                    labels = c("Interior", "Exterior")) +
  labs(x= "Relative abundance (%)",y= NULL) +
  #theme_classic() +
  theme(axis.text.y = element_markdown())

# pairwise wilcox test
pairwise.wilcox.test(composite_V4$rel_abund, composite_V4$group, p.adjust.method = "BH")

multi_sig_genera <-
  composite_V4 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x=data, ~pairwise.wilcox.test(.x$rel_abund, .x$group, exact = FALSE) %>%  tidy)) %>%
  unnest(test) %>%
  filter(p.value < 0.05) %>%
  filter(group2=="pot-inside", group1=="pot-outside")
