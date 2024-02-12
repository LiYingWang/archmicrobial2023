# run packages, context(metadata), and taxonomy data in 000
# follow: # https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html#other_visualizations
rarefy_V4 <- read.csv(here::here("analysis","data","derived_data","V4_shared_format.1.subsample.60216.csv"), sep = ",")
rarefy_V1 <- read.csv(here::here("analysis","data","derived_data","V1_shared_format.1.subsample.13658.csv"), sep = "")
rarefy_V6 <- read.csv(here::here("analysis","data","derived_data","V6_shared_format.1.subsample.28402.csv"), sep = ",")
rarefy_ITS1 <- read.csv(here::here("analysis","data","derived_data","ITS1_shared_format.1.subsample.387.csv"), sep = ",")

# to get the same order for all
rarefy_V4 <- rarefy_V4 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V1 <- rarefy_V1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V6 <- rarefy_V6 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_ITS1 <- rarefy_ITS1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)

# get the specific region context
region_context <- function(df, vec){
  df %>%
    mutate(`adpater no` = tolower(`adpater no`)) %>%
    filter(Region == vec) %>%
    arrange(`adpater no`) %>%
    select(c(1:5, tail(names(.), 5)))
}

V4_context <- region_context(context, "V4")
V1_context <- region_context(context, "V1")
V6_context <- region_context(context, "V6")
ITS1_context <- region_context(context, "ITS1")

# rarefied taxonomy data with sample name and ASV name
tax_rarefy_df <- function(df1, df2, df3){
  df1 %>% #rarefy_V4
    select(-label, -numOtus) %>%
    left_join(df2, by = c("Group"= "adpater no")) %>% #V4_context
    arrange(type, part, order) %>%
    #mutate(label = factor(label, levels = label)) %>%
    select(-Group, -area, -Region) %>%
    pivot_longer(!c(sample, label, type, group, site, part, order), names_to = "name", values_to = "count") %>%
    left_join(df3, by = c("name"= "OTU.ID")) %>% #V4_tax_clean_sep
    mutate(Phylum = ifelse(Phylum == "p__WPS_2", "p__Eremiobacterota", Phylum),
           Phylum = ifelse(Phylum == "p__Deinococcus_Thermus","p__Deinococcus", Phylum))
}

V4_tax_rarefy_df <- tax_rarefy_df(rarefy_V4, V4_context, V4_tax_clean_sep)
V1_tax_rarefy_df <- tax_rarefy_df(rarefy_V1, V1_context, V1_tax_clean_sep)
V6_tax_rarefy_df <- tax_rarefy_df(rarefy_V6, V6_context, V6_tax_clean_sep)
ITS1_tax_rarefy_df <- tax_rarefy_df(rarefy_ITS1, ITS1_context, ITS1_tax_clean_sep)

# get percentage of taxa across group
options(scipen = 999) #convert scientific notation for numbers. If reverse: options(scipen = 0)
V4_taxa_percent_group <-
  V4_tax_rarefy_df %>%
  filter(!is.na(Phylum)) %>%
  group_by(group, Phylum) %>%
  dplyr::summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count) * 100) %>%
  ungroup()

# Archaeal only: relative abundance plot (the top ones with others)
library(RColorBrewer)
V6_archaeal_count_ab  <-
  V6_tax_rarefy_df %>%
  filter(Domain == "k__Archaea") %>%
  ggplot(aes(x = label, y = count, fill = Phylum)) +
  geom_bar(stat="identity", position = "stack", width = 0.8) + #"stack"
  #scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = brewer.pal(7,"Paired")) +
  labs(x = NULL, y = "Abundance") +
  facet_grid(~type, scale= "free", space = "free")

# get the top 11 Phylum based on the counts of 4 groups
V4_Phylum_group <- V4_tax_rarefy_df %>%
  group_by(group) %>%
  slice_max(count, n = 19) %>%
  ungroup() %>%
  distinct(Phylum) %>%
  arrange(Phylum) %>% #in a alphabetic order
  rbind("p__Crenarchaeota") %>%# here to get an order
  pull()

V1_Phylum_group <- V1_tax_rarefy_df %>%
  group_by(group) %>%
  slice_max(count, n = 400) %>%
  ungroup() %>%
  distinct(Phylum) %>%
  arrange(Phylum) %>% #in a alphabetic order
  pull()

V6_Phylum_group <- V6_tax_rarefy_df %>%
  group_by(group) %>%
  slice_max(count, n = 18) %>%
  ungroup() %>%
  distinct(Phylum) %>%
  arrange(Phylum) %>% #in a alphabetic order
  rbind("p__Thaumarchaeota", "p__Nanoarchaeaeota") %>%# here to get an order
  pull()

ITS1_Class_group <- ITS1_tax_rarefy_df %>%
  group_by(group) %>%
  slice_max(count, n = 50) %>%
  ungroup() %>%
  distinct(Class) %>%
  arrange(Class) %>% #in a alphabetic order
  pull()

phylum_colors_rel_mod<- c("#D9D9D9","#FFFF99","#B2DF8A","#33A02C","#1F78B4","#E31A1C", # "C51B7D"
                          "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#B15928","#FB9A99") #brewer.pal(12,"Paired")
phylum_colors_rel_mod2<- c("#D9D9D9","#FFFF99","#B2DF8A","#F781BF","#1F78B4","#E31A1C",
                           "#FF7F00","#CAB2D6","#6A3D9A","#B15928","#A6CEE3","#FB9A99") # for V6
phylum_colors_rel_mod3<- c("#D9D9D9","#FDAE61","#FFED6F","#BC80BD","#80B1D3","#D6604D",
                           "#B3DE69","#8DD3C7","#BEBADA","#CCEBC5","#FCCDE5","#DFC27D") # for ITS1

# V4 relative abundance plot (top ones are others and unclassified)
V4_bac_phyla_rel_ab <- # region number can be replaced
  V4_tax_rarefy_df %>%
  #mutate(type = case_when(type == "con" ~ "control", type == "pot" ~ "pottery", TRUE ~ "Other"))
  mutate(Phylum = ifelse(!Phylum%in% V4_Phylum_group, "others", Phylum)) %>%
  mutate(Taxonomy = ifelse(Domain == "unclassified" & Phylum == "unclassified"|
                           Domain == "k__Bacteria"& Phylum == "unclassified"|
                           Domain == "k__Archaea" & Phylum == "unclassified",
                          "unclassified", paste0(Domain,"; ",Phylum))) %>%
  mutate(Taxonomy = case_when(Taxonomy == "k__Bacteria; others" ~ "others",
                              Taxonomy == "k__Archaea; others" ~ "others", TRUE ~ as.character(Taxonomy))) %>%
  ggplot(aes(x = label, y = count, fill = Taxonomy, color = Taxonomy)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE), #"stack" shows original count; # fill
           width = 0.8, linewidth = 0.1) + # linewidth to reduce the border lines
  scale_y_continuous(labels = scales::percent, expand = c(0.01, 0.01)) +
  scale_fill_manual(values = rev(phylum_colors_rel_mod)) + #breaks = V4_Phylum_group_order
  scale_color_manual(values = rev(phylum_colors_rel_mod)) +
  guides(fill = guide_legend(title = "Taxonomy (phylum level)"), color = "none") +
  labs(x = NULL, y = "Relative Abundance") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10)) +
  facet_grid(~group, scale= "free", space = "free")

# export the df before plotting
write.csv(V6_bac_phyla_rel_ab,
          here::here("analysis", "data", "derived_data", "V6_bac_phyla_rel_ab.csv"),
                     row.names = FALSE)

library(Cairo)
CairoPNG(here::here("analysis","figures","03_rel_ab_V4.png"), width = 10, height = 4, dpi = 360, units = "in")
V4_bac_phyla_rel_ab
dev.off()

# V1 relative abundance plot (top ones are others and unclassified)
V1_bac_phyla_rel_ab <-
  V1_tax_rarefy_df %>%
  mutate(Phylum = ifelse(!Phylum%in% V1_Phylum_group, "others", Phylum)) %>%
  mutate(Taxonomy = ifelse(Domain == "unclassified" & Phylum == "unclassified"|
                           Domain == "k__Bacteria"& Phylum == "unclassified",
                           "unclassified", paste0(Domain,"; ",Phylum))) %>%
  mutate(Taxonomy = case_when(Domain == "k__Archaea" & Phylum == "unclassified" ~ "others",
                              Taxonomy == "k__Bacteria; others" ~ "others",
                              Taxonomy == "k__Archaea; others" ~ "others", TRUE ~ as.character(Taxonomy))) %>%
  ggplot(aes(x = label, y = count, fill = Taxonomy)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE), width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = rev(phylum_colors_rel_mod[c(1:3, 5:7, 9:10)])) +
  guides(fill = guide_legend(title = "Taxonomy (phylum level)")) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(~type, scale= "free", space = "free")

CairoPNG(here::here("analysis","figures", "03_rel_ab_V1.png"), width = 9, height = 4, dpi = 360, units = "in")
V1_bac_phyla_rel_ab
dev.off()

# V6 relative abundance plot (top ones are others and unclassified)
V6_bac_phyla_rel_ab <-
  V6_tax_rarefy_df %>%
  mutate(Phylum = ifelse(!Phylum%in% V6_Phylum_group, "others", Phylum)) %>%
  mutate(Taxonomy = ifelse(Domain == "unclassified" & Phylum == "unclassified"|
                           Domain == "k__Archaea" & Phylum == "unclassified"|
                           Domain == "k__Bacteria"& Phylum == "unclassified",
                           "unclassified", paste0(Domain,"; ",Phylum))) %>%
  mutate(Taxonomy = case_when(Taxonomy == "k__Bacteria; others" ~ "others",
                              Taxonomy == "k__Archaea; others" ~ "others", TRUE ~ as.character(Taxonomy))) %>%
  ggplot(aes(x = label, y = count, fill = Taxonomy)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE), width = 0.8) + #"stack" shows original count; # fill
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = rev(phylum_colors_rel_mod2)) +
  guides(fill = guide_legend(title = "Taxonomy (phylum level)")) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(~type, scale= "free", space = "free")

CairoPNG(here::here("analysis","figures","03_rel_ab_V6.png"), width = 9, height = 4, dpi = 360, units = "in")
V6_bac_phyla_rel_ab
dev.off()

# ITS1 relative abundance plot (top ones are others and unclassified)
ITS1_bac_phyla_rel_ab <-
  ITS1_tax_rarefy_df %>%
  mutate(Class = ifelse(!Class%in% ITS1_Class_group, "others", Class)) %>%
  mutate(Taxonomy = case_when(Domain == "unclassified" & Phylum == "unclassified" & Class == "unclassified" ~ "unclassified",
                              Domain == "k__Fungi" & Phylum == "unclassified" & Class == "unclassified" ~ "unclassified",
                              TRUE ~ as.character(paste0(Phylum,"; ", Class)))) %>%
  ggplot(aes(x = label, y = count, fill = Taxonomy)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE), width = 0.8) + #"stack" shows original count; # fill
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = rev(phylum_colors_rel_mod3)) +
  guides(fill = guide_legend(title = "Taxonomy (class level)")) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(~type, scale= "free", space = "free")

CairoPNG(here::here("analysis","figures","03_rel_ab_ITS1.png"), width = 9.5, height = 4, dpi = 360, units = "in")
ITS1_bac_phyla_rel_ab
dev.off()
