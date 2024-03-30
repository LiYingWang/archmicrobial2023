library(tidyverse)
#library(ape) # Phylogenetics and Evolution package. Tree calculations need to be used with phyloseq
#library(gplots) #calculate and plot Venn diagrams as well as heatmaps
#library(lme4) #linear mixed-effects models like repeated measures analysis
#library(phangorn) #read in mothur-formatted files
#library(phyloseq) #for UniFrac analyses. Organizing, linking, storing, and analyzing of phylogenetic sequencing data
library(ggplot2)
library(vegan) #tools for descriptive community ecology. Bray-Curtis and Jaccard analyses
#library(venneuler)

# read meta data to get context
library(readxl)
context <- read_excel(here::here("analysis","data", "raw_data", "GJB_amp_meta.xlsx"))
context <- context %>%
  mutate(type= tolower(str_extract(`adpater no`, "[[:alpha:]]+(?=\\d)")),
         type= case_when(type == "con" ~ "control", TRUE ~ as.character(type))) %>%
  mutate(group = if_else(!is.na(area)&!area%in% c("belly", "outside belly"),
                         paste0(type, "-", area), type),
         group = case_when(group == "pot-inside"~ "pot-interior",
                           group == "pot-outside"~ "pot-exterior", TRUE ~ as.character(group))) %>%
  mutate(label = case_when(`adpater no` == "con1"~ "C1\nblank", `adpater no` == "con2"~ "C2\nblank",
                           `adpater no` == "Pot1"~ "P1\nGJB", `adpater no` == "Pot2"~ "P2\nGJB",
                           `adpater no` == "Pot3"~ "P3\nGJB", `adpater no` == "Pot4"~ "P4\nGJB",
                           `adpater no` == "Pot5"~ "P5\nGJB", `adpater no` == "Pot11"~ "P6\nGS",
                           `adpater no` == "Pot6"~ "P1\nGJB", `adpater no` == "Pot7"~ "P2\nGJB",
                           `adpater no` == "Pot8"~ "P3\nGJB", `adpater no` == "Pot9"~ "P4\nGJB",
                           `adpater no` == "Pot10"~ "P5\nGJB", `adpater no` == "Pot12"~ "P6\nGS",
                           `adpater no` == "soil1"~ "S1\nGJB", `adpater no` == "soil2"~ "S2\nGJB",
                           `adpater no` == "soil3"~ "S3\nGJB", `adpater no` == "con3"~ "C3\nPCR",
                           `adpater no` == "soil4"~ "S4\nGS", `adpater no` == "soil5"~ "S5\nGS",
                           TRUE ~ as.character(`adpater no`))) %>%
  mutate(part = str_extract(label, "(?<=pot[1-9]\n).*"),
         part = factor(part, levels = c("in", "ex")),
         order = str_extract(label, "[1-9]"),
         type = str_extract(label, "[[:alpha:]]+(?=\\d)"),
         group = factor(group, levels = c("control", "pot-interior", "pot-exterior", "soil")))

# read original count data
OTU_V4 <- read.table(here::here("analysis","data","raw_data","V4_asv_table.txt"), header = T, sep="\t")
OTU_V1 <- read.table(here::here("analysis","data","raw_data","V1_asv_table.txt"), header = T, sep="\t")
OTU_V6 <- read.table(here::here("analysis","data","raw_data","V6_asv_table.txt"), header = T, sep="\t")
OTU_ITS1 <- read.table(here::here("analysis","data","raw_data","ITS1_asv_table.txt"), header = T, sep="\t")

check <- OTU_ITS1 %>%
  select(OTU.ID, taxonomy)

# function to get reads only with ASVs as row names
reads_only <- function(OTU) {
  OTU_df <- OTU %>%
    select(-taxonomy)
  row.names(OTU_df) <-  OTU_df$OTU.ID
  OTU_clean <- OTU_df[,-which(names(OTU_df) %in% c("OTU.ID"))]
  return(OTU_df)
}

OTU_df_V4 <- reads_only(OTU_V4)
OTU_df_V1 <- reads_only(OTU_V1)
OTU_df_V6 <- reads_only(OTU_V6)
OTU_df_ITS1 <- reads_only(OTU_ITS1)

# function to clean up the taxonomy data
tax_clean <- function(df, vec){
  df %>%
    select(-starts_with(vec)) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, "; c__; o__; f__; g__; s__", ""),
           taxonomy = str_replace(taxonomy, "; o__; f__; g__; s__", ""),
           taxonomy = str_replace(taxonomy, "; f__; g__; s__", ""),
           taxonomy = str_replace(taxonomy, "; g__; s__", ""),
           taxonomy = str_replace(taxonomy, "(; s__)([^(; s__)]*)$", ""),
           taxonomy = str_replace(taxonomy, "(; f__)([^(; f__)]*)$", ""),
           taxonomy = str_replace(taxonomy, "(; g__)([^(; g__)]*)$", ""))
  #row.names(V4_tax_clean) <-  V4_tax_clean$OTU.ID
}

V4_tax_clean <- tax_clean(OTU_V4, "V4")
V1_tax_clean <- tax_clean(OTU_V1, "V1")
V6_tax_clean <- tax_clean(OTU_V6, "V6")
ITS1_tax_clean <- tax_clean(OTU_ITS1, "ITS1")

# function to separate taxonomic rank
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

V4_tax_clean_sep <- tax_clean_sep(V4_tax_clean)
V1_tax_clean_sep <- tax_clean_sep(V1_tax_clean)
V6_tax_clean_sep <- tax_clean_sep(V6_tax_clean)
ITS1_tax_clean_sep <- tax_clean_sep(ITS1_tax_clean)

sum(is.na(V6_tax_clean_sep$Species))

# combine tidied reads and taxonomy dataframes
V4_context <- context %>% filter(Region == "V4")
OTU_df_V4_original <-
  OTU_df_V4 %>%
  left_join(V4_tax_clean_sep) %>%
  pivot_longer(!c(OTU.ID, Domain, Phylum, Class, Order, Family, Genus, Species),
               names_to = "adpater no", values_to = "count") %>%
  mutate(`adpater no` = str_remove(`adpater no`, "V4.")) %>%
  left_join(V4_context)

# function to transpose the original count table to shared file format (ASVs as columns)
transpose_OTU <- function(OTU_df) {
  OTU_clean_wider <-
    OTU_df %>%
    pivot_longer(-`OTU.ID`) %>%
    pivot_wider(names_from= `OTU.ID`, values_from= value) %>%
    separate(name, c("label", "Group")) %>%
    mutate(numPhylos = ncol(select(., starts_with("ASV")))) %>%
    relocate(label, Group, numPhylos)
}

OTU_clean_wider_V4 <- transpose_OTU(OTU_df_V4)
OTU_clean_wider_V1 <- transpose_OTU(OTU_df_V1)
OTU_clean_wider_V6 <- transpose_OTU(OTU_df_V6)
OTU_clean_wider_ITS1 <- transpose_OTU(OTU_df_ITS1)

# function to calculate the minimum sum of reads
threshold <- function(df) {
  total <- df %>%
    mutate(total = rowSums(across(contains("ASV")), na.rm = T)) %>%
    select(total)
  return(min(total))
}

threshold(OTU_clean_wider_V4)
threshold(OTU_clean_wider_V1)
threshold(OTU_clean_wider_V6)
threshold(OTU_clean_wider_ITS1)

# make shared file format and saved as csv files
write.csv(OTU_clean_wider_V4,
          file = here::here("analysis","data","derived_data","V4_shared_format.csv"))
write.csv(OTU_clean_wider_V1,
          file = here::here("analysis","data","derived_data","V1_shared_format.csv"))
write.csv(OTU_clean_wider_V6,
          file = here::here("analysis","data","derived_data","V6_shared_format.csv"))
write.csv(OTU_clean_wider_ITS1,
          file = here::here("analysis","data","derived_data","ITS1_shared_format.csv"))

# example: calculate richness and Chao1 using vegan package
data_richness <- estimateR(OTU_clean_wider)
data_evenness <- diversity(OTU_clean_wider) / log(specnumber(OTU_clean_wider))
data_shannon <- diversity(OTU_clean_wider, index = "shannon")
data_alphadiv <- cbind(t(data_richness), data_shannon, data_evenness)
