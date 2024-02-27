# must run independently to get the same order of OTU
library(vegan)
library(ggplot2)
library(viridis)
library(tidyverse)
library(readxl)

# read rarefied  data
rarefy_V4 <- read.csv(here::here("analysis","data","derived_data", "V4_shared_format.1.subsample.60216.csv"), sep = ",")
rarefy_V1 <- read.csv(here::here("analysis","data","derived_data","V1_shared_format.1.subsample.13658.csv"), sep = "")
rarefy_V6 <- read.csv(here::here("analysis","data","derived_data", "V6_shared_format.1.subsample.28402.csv"), sep = ",")
rarefy_ITS1 <- read.csv(here::here("analysis","data","derived_data", "ITS1_shared_format.1.subsample.387.csv"), sep = ",")
context <- read_excel(here::here("analysis","data", "raw_data", "GJB_amp_meta.xlsx"))

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
  mutate(label = case_when(`adpater no` == "pot1"~ "pot1-in", `adpater no` == "pot2"~ "pot2-in",
                           `adpater no` == "pot3"~ "pot3-in", `adpater no` == "pot4"~ "pot4-in",
                           `adpater no` == "pot5"~ "pot5-in", `adpater no` == "pot11"~ "pot6-in",
                           `adpater no` == "pot6"~ "pot1-ex", `adpater no` == "pot7"~ "pot2-ex",
                           `adpater no` == "pot8"~ "pot3-ex", `adpater no` == "pot9"~ "pot4-ex",
                           `adpater no` == "pot10"~ "pot5-ex", `adpater no` == "pot12"~ "pot6-ex",
                           `adpater no` == "soil4"~ "soil4-belly", `adpater no` == "soil5"~ "soil5-belly\nex",
                           TRUE ~ as.character(`adpater no`))) %>%
  mutate(potex_soil = case_when(group == "pot-exterior" ~ "pot-exterior", type== "soil" ~ "soil"),
         potin_soil = case_when(group == "pot-interior" ~ "pot-interior", type== "soil" ~ "soil"),
         pot_soil = case_when(type== "pot" ~ "pot", type== "soil" ~ "soil"),
         pot_con = case_when(type== "pot" ~ "pot", type== "control" ~ "control"),
         potex_con = case_when(group == "pot-exterior" ~ "pot-exterior", type== "control" ~ "control"),
         potin_con = case_when(group == "pot-interior" ~ "pot-interior", type== "control" ~ "control"))

V4_meta<- filter(context_tidy, Region == "V4")
V1_meta<- filter(context_tidy, Region == "V1")
V6_meta<- filter(context_tidy, Region == "V6")
ITS1_meta<- filter(context_tidy, Region == "ITS1")

# wider dataframe based on rarefied OTUs (use this to make all the order consistent)
rarefy_clean_wider <- function(df){
  df %>%
    mutate(Group = tolower(Group)) %>%
    select(-label, -numOtus) %>%
    arrange(Group) %>%
    column_to_rownames("Group")
}

rarefy_V4_clean_wider <- rarefy_clean_wider(rarefy_V4)
rarefy_V1_clean_wider <- rarefy_clean_wider(rarefy_V1)
rarefy_V6_clean_wider <- rarefy_clean_wider(rarefy_V6)
rarefy_ITS1_clean_wider <- rarefy_clean_wider(rarefy_ITS1)

# calculate Bray-Curtis dissimilarity - distance matrix for abundance
bray_calculate <- function(df){
  bray <- vegdist(df, dmethod = "bray")
  pcoa_bray <- cmdscale(bray, k = 2, eig = T)
  return(pcoa_bray)
}

V4_bray <- bray_calculate(rarefy_V4_clean_wider)
V1_bray <- bray_calculate(rarefy_V1_clean_wider)
V6_bray <- bray_calculate(rarefy_V6_clean_wider)
ITS1_bray <- bray_calculate(rarefy_ITS1_clean_wider)

# permutational MANOVA for assessing differences in Bray-Curtis dissimilarity
V4_bray_dist <- vegdist(rarefy_V4_clean_wider, dmethod = "bray")
V1_bray_dist <- vegdist(rarefy_V1_clean_wider, dmethod = "bray")
V6_bray_dist <- vegdist(rarefy_V6_clean_wider, dmethod = "bray")
ITS1_bray_dist <- vegdist(rarefy_ITS1_clean_wider, dmethod = "bray")

library(pairwiseAdonis)
V4_bray_perm <- adonis2(V4_bray_dist ~ group, data = V4_meta, permutations = 1000)
V4_bray_perm_pair <- pairwise.adonis(V4_bray_dist, V4_meta$group)
V1_bray_perm <- adonis2(V1_bray_dist ~ group, data = V1_meta, permutations = 1000)
V1_bray_perm_pair <- pairwise.adonis(V1_bray_dist, V1_meta$group)
V6_bray_perm <- adonis2(V6_bray_dist ~ group, data = V6_meta, permutations = 1000)
V6_bray_perm_pair <- pairwise.adonis(V6_bray_dist, V6_meta$group)
ITS1_bray_perm <- adonis2(ITS1_bray_dist ~ group, data = ITS1_meta, permutations = 1000)
ITS1_bray_perm_pair <- pairwise.adonis(ITS1_bray_dist, ITS1_meta$group)

V4_pv_bray <- c(format(round(V4_bray_perm$R2[1], 2), nsamll= 2), round(V4_bray_perm$`Pr(>F)`[1], 3))
V1_pv_bray <- c(format(round(V1_bray_perm$R2[1], 2), nsamll= 2), ceiling(V1_bray_perm$`Pr(>F)`[1]*2)/20)
V6_pv_bray <- c(format(round(V6_bray_perm$R2[1], 2), nsamll= 2), round(V6_bray_perm$`Pr(>F)`[1], 3))
ITS1_pv_bray <- c(format(round(ITS1_bray_perm$R2[1], 2), nsamll= 2), round(ITS1_bray_perm$`Pr(>F)`[1], 2))

# check pairwise beta dispersion
disp_bray_group <- betadisper(ITS1_bray_dist, ITS1_meta$group) # change data set
permutest(disp_bray_group, pairwise=TRUE, permutations=1000)

# get stress values for dissimilarity. (0.2 is suspect, 0.1 is fair, and 0.05 is good)
stress <- function(df, vec) {
  nmds <- metaMDS(df, distance= vec, k= 2, trymax= 1000)
  s_value <- ifelse(round(nmds$stress, 2) > 0.001, round(nmds$stress, 2),
                    format(round(nmds$stress, 4), scientific = T))
  return(s_value)
}
V4_stress_b <- stress(rarefy_V4_clean_wider, "bray")
V1_stress_b <- stress(rarefy_V1_clean_wider, "bray")
V6_stress_b <- stress(rarefy_V6_clean_wider, "bray")
ITS1_stress_b <- stress(rarefy_ITS1_clean_wider, "bray")
V4_pv_bray <- c(V4_pv_bray, V4_stress_b)
V1_pv_bray <- c(V1_pv_bray, V1_stress_b)
V6_pv_bray <- c(V6_pv_bray, V6_stress_b)
ITS1_pv_bray <- c(ITS1_pv_bray, ITS1_stress_b)

# extract axis positions for each sample & create a dataframe for plotting
beta_df <-function(list, df){
  pcoa_bray_df <- as.data.frame(list$points)
  colnames(pcoa_bray_df) <- c("axis_1", "axis_2")
  pcoa_bray_df$Sample <- rownames(pcoa_bray_df)
 # make a dataframe that combines metadata
  pcoa_bray_for_plot <-
    pcoa_bray_df %>%
    left_join(df, by= c("Sample" = "adpater no")) %>%
    arrange(type, area, order) %>%
    mutate(label = factor(label, levels = label))
  return(pcoa_bray_for_plot)
}

V4_pcoa_bray_for_plot <- beta_df(V4_bray, V4_meta)
V1_pcoa_bray_for_plot <- beta_df(V1_bray, V1_meta)
V6_pcoa_bray_for_plot <- beta_df(V6_bray, V6_meta)
ITS1_pcoa_bray_for_plot <- beta_df(ITS1_bray, ITS1_meta)

# calculate the proportion of variance in the data explained by the first two PCoA axes
V4_bray_eig <- round(100*V4_bray$eig/sum(V4_bray$eig), 2)
V1_bray_eig <- round(100*V1_bray$eig/sum(V1_bray$eig), 2)
V6_bray_eig <- round(100*V6_bray$eig/sum(V6_bray$eig), 2)
ITS1_bray_eig <- round(100*ITS1_bray$eig/sum(ITS1_bray$eig), 2)

# create a PCoA bray plot
library(ggrepel)
library(ggtext)
pcoa_bray_plot <-
  ITS1_pcoa_bray_for_plot %>%
  ggplot(aes(x= axis_1, y= axis_2, color= group)) +
  geom_point(size= 3) +
  geom_text_repel(aes(label = order), size = 3.5,  hjust = 2,
                  min.segment.length = 0.2, show.legend = FALSE) +
  annotate("text", Inf, Inf, label = paste(c("italic(R)^2==" , "P<", "S=="), ITS1_pv_bray), # ITS1: Inf, Inf, c(1.1, 1.1), c(1.25, 3)
           hjust = c(1.1, 1.1, 1.1), vjust = c(1.25, 3, 4.25), parse = T ) + #-Inf, Inf, c(-0.2, -0.2), c(1.25, 3) for others
  scale_colour_viridis_d(option = "magma", begin= 0.1, end= 0.8, direction = -1) + #change color by option = "magma"
  theme_bw() + # theme_bw()
  labs(x= paste0("PCoA", "(", ITS1_bray_eig[1], "%)"), y= paste0("PCoA", "(", ITS1_bray_eig[2], "%)"),
       title = "Bray-Curtis") +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 13), legend.title=element_text(size=rel(1.2)),
        legend.text=element_text(size=rel(1.2)),
        strip.text.x = element_text(size=4),strip.text.y = element_text(size=4))

ggsave(here::here("analysis", "figures", "PcoA_bray_rare.png"), width = 6.5, height = 4.5, units = "in")

# calculate Jaccard dissimilarity - distance matrix for presence/absence
jac_calculate <- function(df){
  #calculate principal coordinates analysis (Bray-Curtis)
  jac <- vegdist(df, dmethod = "jaccard", binary = T)
  pcoa_jac <- cmdscale(jac, k = 2, eig = T, add = T) # add to reduce negative values
  return(pcoa_jac)
}

V4_jac <- jac_calculate(rarefy_V4_clean_wider)
V1_jac <- jac_calculate(rarefy_V1_clean_wider)
V6_jac <- jac_calculate(rarefy_V6_clean_wider)
ITS1_jac <- jac_calculate(rarefy_ITS1_clean_wider)

# permutational MANOVA for assessing differences in Jaccard dissimilarity
V4_jac_dist <- vegdist(rarefy_V4_clean_wider, dmethod = "jaccard")
V1_jac_dist <- vegdist(rarefy_V1_clean_wider, dmethod = "jaccard")
V6_jac_dist <- vegdist(rarefy_V6_clean_wider, dmethod = "jaccard")
ITS1_jac_dist <- vegdist(rarefy_ITS1_clean_wider, dmethod = "jaccard")

V4_jac_perm <- adonis2(V4_jac_dist  ~ group, data = V4_meta, permutations = 1000)
V1_jac_perm <- adonis2(V1_jac_dist  ~ group, data = V1_meta, permutations = 1000)
V6_jac_perm <- adonis2(V6_jac_dist  ~ group, data = V6_meta, permutations = 1000)
ITS1_jac_perm <- adonis2(ITS1_jac_dist  ~ group, data = ITS1_meta, permutations = 1000)

V4_pv_jac <- c(format(round(V4_jac_perm$R2[1], 2), nsamll =2), round(V4_jac_perm$`Pr(>F)`[1], 3))
V1_pv_jac <- c(format(round(V1_jac_perm$R2[1], 2), nsamll =2), ceiling(V1_bray_perm$`Pr(>F)`[1]*2)/20)
V6_pv_jac <- c(format(round(V6_jac_perm$R2[1], 2), nsamll =2), round(V6_jac_perm$`Pr(>F)`[1], 3))
ITS1_pv_jac <- c(format(round(ITS1_jac_perm$R2[1], 2), nsamll =2), round(ITS1_jac_perm$`Pr(>F)`[1], 2))

# check pairwise beta dispersion
disp_jac_group <- betadisper(ITS1_jac_dist, ITS1_meta$group) # change data set
permutest(disp_jac_group, pairwise=TRUE, permutations=1000)

# get stress values for dissimilarity. (0.2 is suspect, 0.1 is fair, and 0.05 is good)
V4_stress_j <- stress(rarefy_V4_clean_wider, "jaccard")
V1_stress_j <- stress(rarefy_V1_clean_wider, "jaccard")
V6_stress_j <- stress(rarefy_V6_clean_wider, "jaccard")
ITS1_stress_j <- stress(rarefy_ITS1_clean_wider, "jaccard")
V4_pv_jac <- c(V4_pv_jac, V4_stress_j)
V1_pv_jac <- c(V1_pv_jac, V1_stress_j)
V6_pv_jac <- c(V6_pv_jac, V6_stress_j)
ITS1_pv_jac <- c(ITS1_pv_jac, ITS1_stress_j)

# make them as data frames (use previous functions)
V4_pcoa_jac_for_plot <- beta_df(V4_jac, V4_meta)
V1_pcoa_jac_for_plot <- beta_df(V1_jac, V1_meta)
V6_pcoa_jac_for_plot <- beta_df(V6_jac, V6_meta)
ITS1_pcoa_jac_for_plot <- beta_df(ITS1_jac, ITS1_meta)

# calculate the proportion of variance in the data explained by the first two PCoA axes
V4_jac_eig <- round(100*V4_jac$eig/sum(V4_jac$eig), 2)
V1_jac_eig <- round(100*V1_jac$eig/sum(V1_jac$eig), 2)
V6_jac_eig <- round(100*V6_jac$eig/sum(V6_jac$eig), 2)
ITS1_jac_eig <- round(100*ITS1_jac$eig/sum(ITS1_jac$eig), 2)

tibble(pe = ITS1_jac_eig,
       axis = 1:length(ITS1_jac_eig)) %>%
  ggplot(aes(x=axis, y =pe)) +
  geom_line()

# create a PCoA jaccard plot
pcoa_jac_plot <-
  ITS1_pcoa_jac_for_plot %>%
  ggplot(aes(x= axis_1, y= axis_2, color= group)) +
  geom_point(size= 3) +
  geom_text_repel(aes(label = order), size = 3.5,  hjust = 2,
                  min.segment.length = 0.2, show.legend = FALSE) +
  annotate("text", -Inf, Inf, label = paste(c("italic(R)^2==" , "P<", "S=="), ITS1_pv_jac),
           hjust = c(-0.2, -0.2, -0.2), vjust = c(1.25, 3, 4.25), parse = T ) +
  scale_colour_viridis_d(option = "magma", begin= 0.1, end= 0.8, direction = -1) + #change color by option = "magma"
  theme_bw() + # theme_bw()
  labs(x= paste0("PCoA1", "(", ITS1_jac_eig[1], "%)"), y= paste0("PCoA2", "(", ITS1_jac_eig[2], "%)"),
       title = "Jaccard") +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 13), legend.position = "none",
        strip.text.x = element_text(size=4), strip.text.y = element_text(size=4))

ggsave(here::here("analysis", "figures", "PcoA_jac_rare.png"), width = 6, height = 4.5, units = "in")

# combine two plots
library(cowplot)
beta_di <-
  plot_grid(pcoa_jac_plot, pcoa_bray_plot, rel_widths = c(0.75, 1),
            labels = c('A', 'B'), label_size = 12)

ggsave(here::here("analysis", "figures", "ITS1_PcoA_plots.png"), width = 11, height = 4, units = "in")

# get the percentage of contribution of ASV to the Bray-Curtis dissimilarity
#simper_result <- simper_modify(rarefy_V4_clean_wider, V4_meta$`adpater no`, permutations = 999)
#saveRDS(simper_result, here::here("analysis", "data", "derived_data", "simper_result.rds"))
simper_result <- readRDS(here::here("analysis", "data", "derived_data", "simper_result.rds"))

