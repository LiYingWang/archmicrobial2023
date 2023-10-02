# run packages, context(metadata), and taxonomy data in 000 and objects in 004-taxa
# Venn Diagram to get the specific taxa from pottery
library(VennDiagram)
library(ggvenn)

# ASVs/species in common or difference across four groups for venn diagram
V4_con_df <- V4_tax_rarefy_df %>% filter(group == "control"&!count == 0&!Species == "unclassified") #based on species
V4_pot_in_df <- V4_tax_rarefy_df %>% filter(group == "pot-interior"&!count == 0&!Species == "unclassified")
V4_pot_ex_df <- V4_tax_rarefy_df %>% filter(group == "pot-exterior"&!count == 0&!Species == "unclassified")
V4_soil_df <- V4_tax_rarefy_df %>% filter(group =="soil"&!count == 0&!Species == "unclassified")

V4_con_list <- as.vector(V4_con_df$Species) # based on species
V4_pot_in_list <- as.vector(V4_pot_in_df$Species)
V4_pot_ex_list<- as.vector(V4_pot_ex_df$Species)
V4_soil_list <- as.vector(V4_soil_df$Species)

V1_con_df <- V1_tax_rarefy_df %>% filter(group == "control"&!count == 0&!Species == "unclassified")
V1_pot_in_df <- V1_tax_rarefy_df %>% filter(group == "pot-interior"&!count == 0&!Species == "unclassified")
V1_pot_ex_df <- V1_tax_rarefy_df %>% filter(group == "pot-exterior"&!count == 0&!Species == "unclassified")
V1_soil_df <- V1_tax_rarefy_df %>% filter(group =="soil"&!count == 0&!Species == "unclassified")

V1_con_list <- as.vector(V1_con_df$Species) # based on species
V1_pot_in_list <- as.vector(V1_pot_in_df$Species)
V1_pot_ex_list<- as.vector(V1_pot_ex_df$Species)
V1_soil_list <- as.vector(V1_soil_df$Species)

V6_con_df <- V6_tax_rarefy_df %>% filter(group == "control"&!count == 0&!Species == "unclassified")
V6_pot_in_df <- V6_tax_rarefy_df %>% filter(group == "pot-interior"&!count == 0&!Species == "unclassified")
V6_pot_ex_df <- V6_tax_rarefy_df %>% filter(group == "pot-exterior"&!count == 0&!Species == "unclassified")
V6_soil_df <- V6_tax_rarefy_df %>% filter(group =="soil"&!count == 0&!Species == "unclassified")

V6_con_list <- as.vector(V6_con_df$Species) # based on species
V6_pot_in_list <- as.vector(V6_pot_in_df$Species)
V6_pot_ex_list<- as.vector(V6_pot_ex_df$Species)
V6_soil_list <- as.vector(V6_soil_df$Species)

ITS1_con_df <- ITS1_tax_rarefy_df %>% filter(group == "control"&!count == 0&!Species == "unclassified")
ITS1_pot_in_df <- ITS1_tax_rarefy_df %>% filter(group == "pot-interior"&!count == 0&!Species == "unclassified")
ITS1_pot_ex_df <- ITS1_tax_rarefy_df %>% filter(group == "pot-exterior"&!count == 0&!Species == "unclassified")
ITS1_soil_df <- ITS1_tax_rarefy_df %>% filter(group =="soil"&!count == 0&!Species == "unclassified")

ITS1_con_list <- as.vector(ITS1_con_df$Species) # based on species
ITS1_pot_in_list <- as.vector(ITS1_pot_in_df$Species)
ITS1_pot_ex_list<- as.vector(ITS1_pot_ex_df$Species)
ITS1_soil_list <- as.vector(ITS1_soil_df$Species)

# Venn diagram
V4_venn <- list(control= V4_con_list , "pot-interior"= V4_pot_in_list,
                "pot-exterior"= V4_pot_ex_list, soil= V4_soil_list)

V1_venn <- list(control= V1_con_list , "pot-interior"= V1_pot_in_list,
                "pot-exterior"= V1_pot_ex_list, soil= V1_soil_list)

V6_venn <- list(control= V6_con_list , "pot-interior"= V6_pot_in_list,
                "pot-exterior"= V6_pot_ex_list, soil= V6_soil_list)

ITS1_venn <- list(control= ITS1_con_list , "pot-interior"= ITS1_pot_in_list,
                "pot-exterior"= ITS1_pot_ex_list, soil= ITS1_soil_list)

ggvenn(ITS1_venn, fill_color = c("#0073C2FF", "#EFC000FF", "#33A02C", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4)+
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 12))
ggsave(here::here("analysis", "figures", "ITS1_venndiag.png"), width = 5, height = 4, units = "in")

# check the unique ones for each group
unique_ASV <- function(target_list, list1, list2, list3, df) {
  setdiff(target_list, c(list1, list2, list3)) %>%
    as.data.frame() %>%
    setNames("Species") %>%
    inner_join(df, by = c("Species"= "Species"))
}

V4_con_uni_ASV_df <- unique_ASV(V4_con_list, V4_pot_in_list, V4_pot_ex_list, V4_soil_list, V4_con_df)
V4_potin_uni_ASV_df <- unique_ASV(V4_pot_in_list, V4_con_list, V4_pot_ex_list, V4_soil_list, V4_pot_in_df)
V4_potex_uni_ASV_df <- unique_ASV(V4_pot_ex_list, V4_con_list, V4_pot_in_list, V4_soil_list, V4_pot_ex_df)
V4_soil_uni_ASV_df <- unique_ASV(V4_soil_list, V4_con_list, V4_pot_in_list, V4_pot_ex_list, V4_soil_df)

V1_con_uni_ASV_df <- unique_ASV(V1_con_list, V1_pot_in_list, V1_pot_ex_list, V1_soil_list, V1_con_df)
V1_potin_uni_ASV_df <- unique_ASV(V1_pot_in_list, V1_con_list, V1_pot_ex_list, V1_soil_list, V1_pot_in_df)
V1_potex_uni_ASV_df <- unique_ASV(V1_pot_ex_list, V1_con_list, V1_pot_in_list, V1_soil_list, V1_pot_ex_df)
V1_soil_uni_ASV_df <- unique_ASV(V1_soil_list, V1_con_list, V1_pot_in_list, V1_pot_ex_list, V1_soil_df)

V6_con_uni_ASV_df <- unique_ASV(V6_con_list, V6_pot_in_list, V6_pot_ex_list, V6_soil_list, V6_con_df)
V6_potin_uni_ASV_df <- unique_ASV(V6_pot_in_list, V6_con_list, V6_pot_ex_list, V6_soil_list, V6_pot_in_df)
V6_potex_uni_ASV_df <- unique_ASV(V6_pot_ex_list, V6_con_list, V6_pot_in_list, V6_soil_list, V6_pot_ex_df)
V6_soil_uni_ASV_df <- unique_ASV(V6_soil_list, V6_con_list, V6_pot_in_list, V6_pot_ex_list, V6_soil_df)

# in both pot-interior and pot-exterior
a <- intersect(V1_pot_in_list, V1_pot_ex_list)
b <- intersect(V1_pot_ex_list, V1_soil_list)
c <- intersect(V1_pot_in_list, V1_con_list)

x <- setdiff(a, c(b, c)) %>%
  as.data.frame() %>%
  setNames("Species") %>%
  inner_join(V1_pot_in_df, by = c("Species"= "Species"))

# barplot for 16S unique ones at Species/Genus level
bar_pot_in <-
  V4_potin_uni_ASV_df %>% # replace the region here
  #distinct(ASV_ID, .keep_all = TRUE)
  filter(!Species == "unclassified") %>%
  ggplot(aes(x= label, fill= Species)) + # Species for V1, V4 & V6
  geom_bar(position="stack") + #stack
  #scale_fill_manual(values = brewer.pal(12,"Set3")) +
  scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) +
  labs(x = NULL, title = "pottery interior") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 12))

bar_pot_ex <-
  V4_potex_uni_ASV_df %>%
  filter(!Species == "s__Aspergillus_nomius"&!Species == "unclassified") %>%
  ggplot(aes(x= label, fill= Species)) +
  geom_bar(position="stack") +
 # scale_fill_manual(values = brewer.pal(12,"Paired")) +   #scale_fill_viridis_d() +
  scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) +
  labs(x = NULL, title = "pottery exterior") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 12))

library(cowplot)
unique_taxa_bar <-
  plot_grid(bar_pot_in, bar_pot_ex, ncol = 1, align = "v", labels = c('A', 'B')) #

ggsave(here::here("analysis", "figures", "V6_uni_taxa.png"), width = 14, height = 10, units = "in") # 9, 6.5

# ITS1
ITS1_con_uni_ASV_df <- unique_ASV(ITS1_con_list, ITS1_pot_in_list, ITS1_pot_ex_list, ITS1_soil_list, ITS1_con_df)
ITS1_potin_uni_ASV_df <- unique_ASV(ITS1_pot_in_list, ITS1_con_list, ITS1_pot_ex_list, ITS1_soil_list, ITS1_pot_in_df)
ITS1_potex_uni_ASV_df <- unique_ASV(ITS1_pot_ex_list, ITS1_con_list, ITS1_pot_in_list, ITS1_soil_list, ITS1_pot_ex_df)
ITS1_soil_uni_ASV_df <- unique_ASV(ITS1_soil_list, ITS1_con_list, ITS1_pot_in_list, ITS1_pot_ex_list, ITS1_soil_df)

bar_pot_in <-
  ITS1_potin_uni_ASV_df %>% # replace the region here
  #distinct(ASV_ID, .keep_all = TRUE)
  filter(!Species == "s__Aspergillus_nomius"&!Order == "unclassified") %>%
  ggplot(aes(x= label, y= Order)) + # Species for V1, V4 & V6
  geom_point() + #stack
  geom_text(aes(Order))+
  #scale_fill_manual(values = brewer.pal(12,"Set3")) +
  scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) +
  labs(x = NULL, title = "pottery interior") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size= 12))

