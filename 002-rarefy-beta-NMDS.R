library(tidyverse)
library(ggplot2)
library(vegan) #tools for descriptive community ecology. Bray-Curtis and Jaccard analyses
library(readxl)
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

# tidy meta data and get the same order as the rarefy data
context_for_join <-
  context %>%
  mutate(`adpater no` = tolower(`adpater no`),
         type = tolower(str_extract(`adpater no`, "[[:alpha:]]+(?=\\d)")),
         type = case_when(type == "con" ~ "control", TRUE ~ as.character(type)),
         area = case_when(area == "inside" ~ "interior", area == "outside" ~ "exterior"),
         group = if_else(!is.na(area)&!area%in% c("belly", "outside belly"),
                         paste0(type, "-", area), type),
         group = factor(group, levels = c("control","pot-interior","pot-exterior","soil"))) %>%
  arrange(`adpater no`)

# to get the same order for all
rarefy_V4 <- rarefy_V4 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V1 <- rarefy_V1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_V6 <- rarefy_V6 %>% mutate(Group = tolower(Group)) %>% arrange(Group)
rarefy_ITS1 <- rarefy_ITS1 %>% mutate(Group = tolower(Group)) %>% arrange(Group)

# function to have OTU reads only
reads_only <- function(df) {
  df %>%
    select(-label, -Group, -numOtus)
}

rarefy_V4_reads_only  <- reads_only(rarefy_V4)
rarefy_V1_reads_only  <- reads_only(rarefy_V1)
rarefy_V6_reads_only  <- reads_only(rarefy_V6)
rarefy_ITS1_reads_only  <- reads_only(rarefy_ITS1)

# function to join the meta data
join_meta<- function(df1, df2, region){
  df1 %>%
    left_join(filter(df2, Region == region), by =c("Group" = "adpater no")) %>%
    select(-contains("ASV")) %>%
    relocate(Region, Group, type, group, sample, area, site) %>%
    rownames_to_column()
}

rarefy_V4_meta_join <- join_meta(rarefy_V4, context_for_join , region = "V4")
rarefy_V1_meta_join <- join_meta(rarefy_V1, context_for_join, region = "V1")
rarefy_V6_meta_join <- join_meta(rarefy_V6, context_for_join, region = "V6")
rarefy_ITS1_meta_join <- join_meta(rarefy_ITS1, context_for_join, region = "ITS1")

# assessing differences in community composition using permutational MANOVA
V4_pot_perm <- adonis2(rarefy_V4_reads_only ~ group, data = rarefy_V4_meta_join) # microbes differ by group
V1_pot_perm <- adonis2(rarefy_V1_reads_only ~ group, data = rarefy_V1_meta_join)
V6_pot_perm <- adonis2(rarefy_V6_reads_only ~ group, data = rarefy_V6_meta_join)
ITS1_pot_perm <- adonis2(rarefy_ITS1_reads_only ~ group, data = rarefy_ITS1_meta_join)

# Non-metric Multidimensional Scaling (NMDS): collapse all species axes into 2 to visualize the differences between samples
V4_NMDS <- metaMDS(rarefy_V4_reads_only)
V1_NMDS <- metaMDS(rarefy_V1_reads_only)
V6_NMDS <- metaMDS(rarefy_V6_reads_only)
ITS1_NMDS <- metaMDS(rarefy_ITS1_reads_only)
stressplot(ITS1_NMDS)
plot(ITS1_NMDS)

# function to make a dataframe for NMDS
NMDS_df <- function(metaMDS, df) {
  df <- scores(metaMDS, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    full_join(df, by = c("name" = "rowname"))
}

V4_plot_NMDS_df <- NMDS_df(V4_NMDS , rarefy_V4_meta_join)
V1_plot_NMDS_df <- NMDS_df(V1_NMDS , rarefy_V1_meta_join)
V6_plot_NMDS_df <- NMDS_df(V6_NMDS , rarefy_V6_meta_join)
ITS1_plot_NMDS_df <- NMDS_df(ITS1_NMDS , rarefy_ITS1_meta_join)
NMDS_df_all <- rbind(V1_plot_NMDS_df, V4_plot_NMDS_df, V6_plot_NMDS_df, ITS1_plot_NMDS_df)

# NMDS plots
nmds_plot_all <-
  ggplot(NMDS_df_all,
         aes(x = NMDS1, y = NMDS2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(linetype = 2, size = 1) +
  labs(title = "NMDS") +
  scale_color_discrete(name = "Sample group",
                      labels = c("control", "pot-interior", "pot-exterior", "soil")) +
  scale_shape_discrete(name = "Sample group",
                       labels = c("control", "pot-interior", "pot-exterior", "soil")) +
  theme_minimal() +
  facet_wrap(~Region, scales = "free") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 12))

ggsave(here::here("analysis", "figures", "03-NMDS-plot.png"),
       width = 10, height = 8, units = "in")

# How do species contribute to the dissimilarity of communities? It takes an hour
# fit <- envfit(V4_pot_NMDS, rarefy_V4_reads_only, perm = 999) # takes the output of metaMDS() and the species matrix
# saveRDS(fit, here::here("analysis","data","derived_data","fit.rds"))
fit <- readRDS(here::here("analysis","data","derived_data","fit.rds"))
# extract p-values for each species
fit_pvals <- fit$vectors$pvals %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  dplyr::rename("pvals" = ".")

# extract coordinates for species, only keep species with p-val = 0.001
fit_spp <- fit %>%
  scores(., display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  full_join(.,fit_pvals, by = "species") %>%
  filter(pvals == 0.001)

# new plot
nmds_plot_new <- ggplot(V4_plot_NMDS_df, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() +
  geom_point(aes(color = group, shape = group), size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = group)) +
  geom_segment(data = fit_spp, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               col = "black") +
  geom_text(data = fit_spp, aes(label = species))
nmds_plot_new

# How is community structure related to specific variables?
V4_potCCA <- cca(rarefy_V4_reads_only ~ type + site + group,
              data = rarefy_V4_meta_join)
V4_potCCA
plot(V4_potCCA)

# vectors
V4_ccavectors <-
  as.matrix(scores(V4_potCCA, display = "bp", scaling = "species")*3) %>%
  as.data.frame()

# type coordinates
V4_pot_meta_data <-
  scores(V4_potCCA, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("name") %>%
  mutate(name = str_remove(name, "row")) %>%
  full_join(rarefy_V4_meta_join, by = c("name" = "rowname"))

V4_species_data <-
  scores(V4_potCCA, display = "species") %>%
  as.data.frame()

V4_plot_cca <- ggplot(V4_pot_meta_data) +
  geom_point(aes(x = CCA1, y = CCA2, color = group), shape = 19, size = 2, alpha = 0.8) +
  coord_fixed() +
  geom_segment(data = V4_ccavectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = V4_species_data, aes(x = CCA1, y = CCA2), shape = 17, size = 1, color = "gray") +
  #scale_x_continuous(limits = c(-12, 10)) +
  #scale_y_continuous(limits = c(-3, 12)) +
  geom_text(data = V4_ccavectors, aes(x = CCA1, y = CCA2, label = rownames(V4_ccavectors)), nudge_x = 0.3, nudge_y = 0.3) +
  labs(title = "Canonical Correspondence Analysis")
V4_plot_cca

#################################### Tutorial Practice ########
# run OTUs, metadata, and taxonomy data in 000-tidy-data, and rarefy data in 001-rarefy
#library(reshape)
library(tidyverse)
library(ggplot2)
library(vegan)
#library(hrbrthemes)
library(ggtext)

# follow https://www.youtube.com/watch?v=3iOTGaWQZp4
# wider dataframe with original OTUs
OTU_clean_wider_V4 <-
  OTU_clean_wider_V4 %>% # in 000 file
  mutate(Group = tolower(Group)) %>%
  select(-label, -numPhylos) %>%
  column_to_rownames("Group")

# wider dataframe with rarefied OTUs
rarefy_V4_clean_wider <-
  rarefy_V4 %>%
  mutate(Group = tolower(Group)) %>%
  select(-label, -numOtus) %>%
  column_to_rownames("Group")

# make a distance matrix for heatmap
rarefy_V4_clean_wider_tidy <-
  rarefy_V4_clean_wider %>%
  #rownames_to_column(rarefy_V4_clean_wider) %>%
  mutate(rowname = case_when(rowname == "pot1"~ "pot1\nin", rowname == "pot2"~ "pot2\nin",
                             rowname == "pot3"~ "pot3\nin", rowname == "pot4"~ "pot4\nin",
                             rowname == "pot5"~ "pot5\nin", rowname == "pot11"~ "pot6\nin",
                             rowname == "pot6"~ "pot1\nex", rowname == "pot7"~ "pot2\nex",
                             rowname == "pot8"~ "pot3\nex", rowname == "pot9"~ "pot4\nex",
                             rowname == "pot10"~ "pot5\nex", rowname == "pot12"~ "pot6\nex",
                             rowname == "soil4"~ "soil4\nbelly", rowname == "soil5"~ "soil5\nbelly\nex",
                             TRUE ~ as.character(rowname))) %>%
  column_to_rownames()

bray <-
  avgdist(rarefy_V4_clean_wider_tidy , dmethod = "jaccard", sample = 60216) %>%
  as.matrix() %>%
  as_tibble(rownames = "A") %>%
  mutate(part = str_extract(A, "(?<=pot[1-9]\n).*"),
         part = factor(part, levels = c("in", "ex")),
         order = str_extract(A, "[1-9]"),
         type = str_extract(A, "[[:alpha:]]+(?=\\d)")) %>%
  arrange(type, part, order) %>%
  select(-type, -part, -order) %>%
  rownames_to_column() %>%
  pivot_longer(c(-A, -rowname), names_to = "B", values_to = "distances")

jaccard <-
  avgdist(rarefy_V4_clean_wider_tidy , dmethod = "jaccard", sample = 60216) %>%
  as.matrix() %>%
  as_tibble(rownames = "A") %>%
  mutate(part = str_extract(A, "(?<=pot[1-9]\n).*"),
         part = factor(part, levels = c("in", "ex")),
         order = str_extract(A, "[1-9]"),
         type = str_extract(A, "[[:alpha:]]+(?=\\d)")) %>%
  arrange(type, part, order) %>%
  select(-type, -part, -order) %>%
  rownames_to_column() %>%
  pivot_longer(c(-A, -rowname), names_to = "B", values_to = "distances")

# join bray and jaccard

labels <- tibble(
  x= c(5, 16), y=c(5, 16),
  label = c("Bray-Curtis", "Jaccard")
)

bray_jaccard <-
  inner_join(bray, jaccard, by = c("A", "B")) %>%
  mutate(x = rep(c(1:4, 14, 9, 15, 5:8, 10:13,16:20),length.out = 400)) %>%
  select(rowname.x, A, x, B, bray = distances.x, jaccard = distances.y) %>%
  mutate(distances = if_else(as.numeric(rowname.x) < as.numeric(x), bray, jaccard)) %>%
  mutate(A = fct_reorder(A, as.numeric(rowname.x))) %>%
  mutate(B = fct_reorder(B, -as.numeric(x))) %>%
  ggplot(aes(x =A, y = B, fill = distances)) +
  geom_tile() +
  geom_text(data = labels, aes(x= x , y= y, label = label), inherit.aes = F, size = 8) +
  scale_fill_gradient(low="#FF0000", high = "#FFFFFF") +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(hjust = 0.5))

ggsave(here::here("analysis", "figures", "02_beta_dis_matrix.png"),
       width = 8, height = 7, units = "in")
