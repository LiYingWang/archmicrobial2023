# run 004 code first to get rarefied data
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
  theme(axis.text.y = element_text(size = 5)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "ITS1 region: ASV level") +
  geom_segment(data=sep_line, aes(x,y,xend=xend, yend=yend), linewidth=0.5, color= "white", inherit.aes=F) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

CairoPNG(here::here("analysis","figures","04-ITS1-ASV-heatmap.png"), width = 11, height = 6, dpi = 360, units = "in")
ITS1_heatmap_plot
dev.off()

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
  theme(axis.text.y = element_text(size = 5)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "ITS1 region: Genus level") +
  geom_segment(data=sep_line, aes(x,y,xend=xend, yend=yend), linewidth=0.5, color= "white", inherit.aes=F) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

CairoPNG(here::here("analysis","figures","04-ITS1-genus-heatmap.png"), width = 11, height = 6, dpi = 360, units = "in")
ITS1_tax_genus_heatmap
dev.off()

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
  theme(axis.text.y = element_text(size = 5)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "Family level") +
  geom_segment(data=sep_line, aes(x,y,xend=xend, yend=yend), linewidth=0.5, color= "white", inherit.aes=F) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")

# target specious: Aspergillus spp., Penicillium spp., Fusarium spp., Rhizopus spp.
ITS1_target_genus_heatmap <-
  ITS1_tax_rarefy_df %>%
  filter(Genus %in% c("g__Aspergillus", "g__Penicillium", "g__Fusarium")) %>%
  group_by(label, Genus) %>%
  mutate(sum = sum(count)) %>%
  distinct(label, Genus, sum) %>%
  ungroup() %>%
  mutate(z_score = (sum-mean(sum))/sd(sum)) %>%
  ggplot(aes(x = label, y = Genus, fill = z_score)) +
  geom_tile() + #color='White', linewidth=0.1
  scale_fill_viridis_c(begin= 0.1, end= 0.9) + #direction= -1; #scale_fill_distiller(palette = "RdPu") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "Genus level") +
  geom_segment(data=sep_line, aes(x,y,xend=xend, yend=yend), linewidth=0.5, color= "white", inherit.aes=F) +
  theme_minimal() +
  theme(axis.text.y = element_markdown(), legend.position="left")
