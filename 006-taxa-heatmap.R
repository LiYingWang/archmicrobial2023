library(ggtext)
# function to get the top ASVs based on the counts by 4 groups
top_freq_ASV <- function(df, num){
  df %>% # from 004
    #filter(!Genus == "unclassified") %>%
    group_by(group) %>%
    slice_max(count, n = num) %>%
    ungroup() %>%
    distinct(name) %>%
    arrange(name) %>%
    pull(name)
}
# function to get the top Genus based on the counts by 4 groups
top_freq_genus <- function(df, num){
  df %>% # from 004
    filter(!Genus == "unclassified") %>%
    group_by(group) %>%
    slice_max(count, n =num) %>%
    pull(Genus)
}
# function for target family heatmap
target_family_heatmap <- function(df){
  df %>%
    filter(Family %in% c("f__Leuconostocaceae", "f__Enterococcaceae","f__Enterobacteriaceae","f__Phycisphaeraceae", "f__Paenibacillaceae")) %>%
    group_by(label, Family) %>%
    mutate(sum = sum(count)) %>%
    distinct(label, Family, sum) %>%
    ungroup() %>%
    mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
    ggplot(aes(x = label, y = Family, fill = z_score)) +
    geom_tile() + #color='White', linewidth=0.1
    scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
    theme(axis.text.y = element_text(size = 6)) +
    coord_equal() +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_markdown(), legend.position="left")
}

#####################  V4: get the top ASVs based on the counts by 4 groups
V4_tax_group_for_heatmap <- top_freq_ASV(V4_tax_rarefy_df, 15)
# plot heat map (ASV level)
V4_heatmap_plot <-
  V4_tax_rarefy_df %>%
  filter(name %in% c(V4_tax_group_for_heatmap)) %>%
  mutate(Taxa = ifelse(!Species == "unclassified", Species, Genus)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Family, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Order, Taxa)) %>%
  #mutate(Genus = str_remove(Genus, "g__"),
        # Genus = str_replace(Genus, "(.*)", "*\\1*")) %>% #label = str_replace(label, "([0-9]+)", "\\1\n")
  mutate(z_scores = (count-mean(count))/sd(count)) %>% #heatmap based on count table (z-score)
  ggplot(aes(x = label, y = Taxa, fill = z_scores)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# genus level and heatmap
V4_tax_Genus_for_heatmap <- top_freq_genus(V4_tax_rarefy_df, 10) # "g__Candidatus_Xiphinematobacter"
V4_tax_genus_heatmap <-
  V4_tax_rarefy_df %>%
  filter(Genus %in% V4_tax_Genus_for_heatmap) %>%
  group_by(label, Genus) %>%
  mutate(sum = sum(count)) %>%
  distinct(label, Genus, sum) %>%
  ungroup() %>%
  mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
  filter(!Genus %in% c("g__Thermoanaerobacterium","g__Deinococcus","g__Thermoanaerobacterium","g__Sphingomonas","g__Fimbriiglobus")) %>%
  ggplot(aes(x = label, y = Genus, fill = z_score)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# get the target Family based on references and taxa testing
V4_target_family_heatmap <- target_family_heatmap(V4_tax_rarefy_df)

#####################  V1: get the top ASVs based on the counts by 4 groups
V1_tax_group_for_heatmap <- top_freq_ASV(V4_tax_rarefy_df, 50)
V1_heatmap_plot <-
  V1_tax_rarefy_df %>%
  filter(name %in% c(V1_tax_group_for_heatmap)& !Phylum == "unclassified") %>%
  mutate(Taxa = ifelse(!Species == "unclassified", Species, Genus)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Family, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Order, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Class, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Phylum, Taxa)) %>%
  #mutate(Genus = str_remove(Genus, "g__"),
  # Genus = str_replace(Genus, "(.*)", "*\\1*")) %>% #label = str_replace(label, "([0-9]+)", "\\1\n")
  mutate(z_scores = (count-mean(count))/sd(count)) %>% #heatmap based on count table (z-score)
  ggplot(aes(x = label, y = Taxa, fill = z_scores)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

#Genus level and heatmap
V1_tax_Genus_for_heatmap <- top_freq_genus(V1_tax_rarefy_df, 5)
V1_tax_genus_heatmap <-
  V1_tax_rarefy_df %>%
  filter(Genus %in% V1_tax_Genus_for_heatmap) %>%
  group_by(label, Genus) %>%
  mutate(sum = sum(count)) %>%
  distinct(label, Genus, sum) %>%
  ungroup() %>%
  mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
  #filter(!Genus %in% c("g__Thermoanaerobacterium","g__Deinococcus","g__Thermoanaerobacterium","g__Sphingomonas","g__Fimbriiglobus")) %>%
  ggplot(aes(x = label, y = Genus, fill = z_score)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# get the target Family based on references and taxa testing
V1_target_family_heatmap <- target_family_heatmap(V1_tax_rarefy_df)

##################### V6: get the top ASVs based on the counts by 4 groups
V6_tax_group_for_heatmap <- top_freq_ASV(V6_tax_rarefy_df, 15)
# plot heat map (ASV level)
V6_heatmap_plot <-
  V6_tax_rarefy_df %>%
  filter(name %in% c(V6_tax_group_for_heatmap)& !Phylum == "unclassified") %>%
  mutate(Taxa = ifelse(!Species == "unclassified", Species, Genus)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Family, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Order, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Class, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Phylum, Taxa)) %>%
  #mutate(Genus = str_remove(Genus, "g__"),
  # Genus = str_replace(Genus, "(.*)", "*\\1*")) %>% #label = str_replace(label, "([0-9]+)", "\\1\n")
  mutate(z_scores = (count-mean(count))/sd(count)) %>% # heatmap based on count table (z-score)
  ggplot(aes(x = label, y = Taxa, fill = z_scores)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

#Genus level and heatmap
V6_tax_Genus_for_heatmap <- top_freq_genus(V6_tax_rarefy_df, 5)
V6_tax_genus_heatmap <-
  V6_tax_rarefy_df %>%
  filter(Genus %in% V6_tax_Genus_for_heatmap) %>%
  group_by(label, Genus) %>%
  mutate(sum = sum(count)) %>%
  distinct(label, Genus, sum) %>%
  ungroup() %>%
  mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
  #filter(!Genus %in% c("g__Thermoanaerobacterium","g__Deinococcus","g__Thermoanaerobacterium","g__Sphingomonas","g__Fimbriiglobus")) %>%
  ggplot(aes(x = label, y = Genus, fill = z_score)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# get the target Family based on references and taxa testing
V6_target_family_heatmap <- target_family_heatmap(V6_tax_rarefy_df)

#####################  ITS1: get the top ASVs based on the counts by 4 groups
ITS1_tax_group_for_heatmap <- top_freq_ASV(ITS1_tax_rarefy_df, 15)
# plot heat map (ASV level)
ITS1_heatmap_plot <-
  ITS1_tax_rarefy_df %>%
  filter(name %in% c(ITS1_tax_group_for_heatmap)& !Phylum == "unclassified") %>%
  mutate(Taxa = ifelse(!Species == "unclassified", Species, Genus)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Family, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Order, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Class, Taxa)) %>%
  mutate(Taxa = ifelse(Taxa == "unclassified", Phylum, Taxa)) %>%
  #mutate(Genus = str_remove(Genus, "g__"),
  # Genus = str_replace(Genus, "(.*)", "*\\1*")) %>% #label = str_replace(label, "([0-9]+)", "\\1\n")
  mutate(z_scores = (count-mean(count))/sd(count)) %>% #heatmap based on count table (z-score)
  ggplot(aes(x = label, y = Taxa, fill = z_scores)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

#Genus level and heatmap
ITS1_tax_Genus_for_heatmap <- top_freq_genus(ITS1_tax_rarefy_df, 5)
ITS1_tax_genus_heatmap <-
  ITS1_tax_rarefy_df %>%
  filter(Genus %in% ITS1_tax_Genus_for_heatmap) %>%
  group_by(label, Genus) %>%
  mutate(sum = sum(count)) %>%
  distinct(label, Genus, sum) %>%
  ungroup() %>%
  mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
  #filter(!Genus %in% c("g__Thermoanaerobacterium","g__Deinococcus","g__Thermoanaerobacterium","g__Sphingomonas","g__Fimbriiglobus")) %>%
  ggplot(aes(x = label, y = Genus, fill = z_score)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 6)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# get the target Family based on references and taxa testing
ITS1_tax_family_for_heatmap <-
  ITS1_tax_rarefy_df %>% # from 004
    filter(!Family == "unclassified") %>%
    group_by(group) %>%
    slice_max(count, n = 6) %>%
    pull(Family)

ITS1_target_family_heatmap <-
ITS1_tax_rarefy_df %>%
    filter(Family %in% ITS1_tax_family_for_heatmap) %>%
    group_by(label, Family) %>%
    mutate(sum = sum(count)) %>%
    distinct(label, Family, sum) %>%
    ungroup() %>%
    mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
    ggplot(aes(x = label, y = Family, fill = z_score)) +
    geom_tile() + #color='White', linewidth=0.1
    scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
    theme(axis.text.y = element_text(size = 6)) +
    coord_equal() +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_markdown(), legend.position="left")

