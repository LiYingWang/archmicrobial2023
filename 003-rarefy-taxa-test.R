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
           pot_soil = case_when(type== "pot" ~ "pot", type== "soil" ~ "soil"),
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
sig_genera_V4 <- # replace the region
  composite_V4 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~ pot, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)
  #mutate(p.adjust = p.adjust(p.value, method="fdr")) %>% #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
  #filter(p.adjust < 0.1)
  #select(taxonomy, p.adjust)

sig_genera_V1 <-
  composite_V1 %>%
  nest(data = -taxonomy) %>% #c(-taxonomy, -otu)
  mutate(test = map(.x= data, ~wilcox.test(rel_abund~  pot_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

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
  filter(p.value < 0.05)

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
pairwise.wilcox.test(composite_V4$rel_abund, composite_V4$group, p.adjust.method = "fdr")

multi_sig_genera <-
  composite_V4 %>%
  nest(data = -taxonomy) %>% # c(-taxonomy, -otu)
  mutate(test = map(.x=data, ~pairwise.wilcox.test(.x$rel_abund, .x$group, exact = FALSE) %>%  tidy)) %>%
  unnest(test) %>%
  filter(p.value < 0.05) %>%
  filter(group2=="pot-inside", group1=="pot-outside")

# function to separate the taxonomy data, also in 000-tidy-data
tax_clean_sep <- function(df){
  df %>%
    separate(taxonomy,
             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
    mutate(Phylum = str_trim(Phylum), Class = str_trim(Class), Order = str_trim(Order),
           Family = str_trim(Family), Genus = str_trim(Genus), Species = str_trim(Species)) %>%
    mutate(Domain = ifelse(Domain == "Unassigned", "unclassified", Domain),
           Phylum = ifelse(is.na(Phylum), "unclassified", Phylum),
           Class = ifelse(is.na(Class), "unclassified", Class),
           Order = ifelse(is.na(Order), "unclassified", Order),
           Family = ifelse(is.na(Family), "unclassified", Family),
           Genus = ifelse(is.na(Genus), "unclassified", Genus),
           Species = ifelse(is.na(Species), "unclassified", Species)) %>%
    mutate(Phylum = ifelse(Phylum == "p__WPS_2", "p__Eremiobacterota", Phylum),
           Phylum = ifelse(Phylum == "p__Deinococcus_Thermus","p__Deinococcus", Phylum))
}

composite_V4_sep <- tax_clean_sep(composite_V4)
composite_V6_sep <- tax_clean_sep(composite_V6)
composite_V1_sep <- tax_clean_sep(composite_V1)
composite_ITS1_sep <- tax_clean_sep(composite_ITS1)

# function for adding up relative abundance based on different level
composite_sep_genera <- function(df, col){ # col= Family, Genus etc.
  df %>%
    group_by(Group, {{col}}) %>%
    mutate(total_genera_sum = sum(rel_abund)) %>%
    ungroup() %>%
    distinct(Group, {{col}}, group, total_genera_sum, pot_con, pot,
             potex_con, potin_con, potin_soil, potex_soil, pot_soil)
}

composite_V4_sep_family <- composite_sep_genera(composite_V4_sep, Family)
composite_V4_sep_genus <- composite_sep_genera(composite_V4_sep, Genus)
composite_V4_sep_species <- composite_sep_genera(composite_V4_sep, Species)
composite_V6_sep_family <- composite_sep_genera(composite_V6_sep, Family)
composite_V6_sep_genus <- composite_sep_genera(composite_V6_sep, Genus)
composite_V6_sep_species <- composite_sep_genera(composite_V6_sep, Species)
composite_V1_sep_family <- composite_sep_genera(composite_V1_sep, Family)
composite_V1_sep_genus <- composite_sep_genera(composite_V1_sep, Genus)
composite_V1_sep_species <- composite_sep_genera(composite_V1_sep, Species)
composite_ITS1_sep_family <- composite_sep_genera(composite_ITS1_sep, Family)
composite_ITS1_sep_genus <- composite_sep_genera(composite_ITS1_sep, Genus)
composite_ITS1_sep_species <- composite_sep_genera(composite_ITS1_sep, Species)

# get significant genera at family level between pottery and control
sig_family_ITS1 <-
  composite_ITS1_sep_family %>% # replace the region
  nest(data = -Family) %>%
  mutate(test = map(.x= data, ~wilcox.test(total_genera_sum~ potex_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

# get significant genera at family level between pottery
sig_family_ITS1_pot <-
  composite_ITS1_sep_family %>%
  nest(data = -Family) %>%
  mutate(test = map(.x= data, ~wilcox.test(total_genera_sum~ pot, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

# get significant genera at genus level between pottery and control
sig_genus_ITS1 <-
  composite_ITS1_sep_genus %>%
  nest(data = -Genus) %>%
  mutate(test = map(.x= data, ~wilcox.test(total_genera_sum~ potin_con, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

# get significant genera at genus level between pottery
sig_genus_ITS1_pot <-
  composite_ITS1_sep_genus %>%
  nest(data = -Genus) %>%
  mutate(test = map(.x= data, ~wilcox.test(total_genera_sum~ pot, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

# get significant genera at species level
sig_species_ITS1 <-
  composite_ITS1_sep_species %>%
  nest(data = -Species) %>%
  mutate(test = map(.x= data, ~wilcox.test(total_genera_sum~ pot, data=.x, exact = FALSE) %>% tidy)) %>% #run broom first
  unnest(test) %>%
  filter(p.value < 0.05)

# plot
composite_V4_pot_control <-
  composite_V4_sep_family %>%
  inner_join(sig_family_V4, by="Family") %>%
  mutate(Family = str_remove(Family, "f__"),
         Family = str_replace(Family, "(.*)", "*\\1*")) %>%
  filter(!is.na(pot_con)) %>%
  mutate(total_genera_sum = 100 * (total_genera_sum + 1/15000), # make it percentage and add a very small number
         pot_con = factor(pot_con, levels = c("control", "pot"))) %>%
  ggplot(aes(x= total_genera_sum, y = Family, fill= pot_con)) +
  geom_boxplot(aes(fill = pot_con), show.legend = TRUE) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("control", "pot"),
                     values = c("ivory3", "#238A8DFF"),
                     labels = c("control", "pottery")) +
  scale_fill_manual(NULL,
                    breaks = c("control", "pot"),
                    values = c("ivory3", "#238A8DFF"),
                    labels = c("control", "pottery")) +
  labs(x= "Relative abundance (%)",y= NULL) +
  theme(axis.text.y = element_markdown())

sig_family_all <- rbind(sig_family_V4, sig_family_V4_pot)
sig_family_name <- distinct(sig_family_all, Family)

composite_V4_combine_pot_con <-
  composite_V4_sep_family %>%
  inner_join(sig_family_name, by = "Family") %>%
  mutate(Family = str_remove(Family, "f__"),
         Family = str_replace(Family, "(.*)", "*\\1*")) %>%
  mutate(total_genera_sum = 100 * (total_genera_sum + 1/15000), # make it percentage and add a very small number
         pot = factor(group, levels = c("con", "pot-inside", "pot-outside", "soil"))) %>%
  ggplot(aes(x= total_genera_sum, y = Family, fill= group)) +
  geom_boxplot(aes(fill = group), show.legend = TRUE) +
  scale_x_log10() +
  scale_fill_viridis_d(breaks = c("con", "pot-inside", "pot-outside", "soil"),
                       labels = c("control","pot-interior", "pot-exterior", "soil"),
                       direction = -1) +
  scale_y_discrete(limits=rev) +
  labs(x= "Relative abundance (%)", y= NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(size= 12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

library(cowplot)
V4_rel_sig_taxa_combine <-
  plot_grid(V4_bac_phyla_rel_ab , composite_V4_combine_pot_con, # run 004-taxa-relative-ab to get the object
            labels = "AUTO", ncol = 1)

ggsave(here::here("analysis","figures","V4_rel_sig_taxa_combine.png"),
       width = 9.5, height = 7, dpi = 360, units = "in")
