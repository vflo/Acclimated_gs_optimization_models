rm(list=ls())

library(tidyverse)
library(TrEvol)
library(V.PhyloMaker)
library(Taxonstand)

#LOAD DATASET
data_traits_total <- read_csv("data/XFT_full_database_compressed.csv", 
                        col_types = cols(P50..MPa. = col_double(),
                                         Longitude = col_double()))
data_traits_ini <- data_traits_total %>% 
  dplyr::select(-References) %>% 
  # mutate(Cleaned.binomial = case_when(is.na(Cleaned.binomial)~paste(genus,species),
  #                                     TRUE ~ Cleaned.binomial),
  #        taxon = gsub(" ", "_", Cleaned.binomial))%>% 
  filter(!is.na(MATbest),!is.na(PPTbest))
# data_traits <- data_traits[c(1:400),]
#CREATE DATASET OF TARGET SPECIES
species_new <- read_csv("data/Species_parameters_template.csv")%>% 
  dplyr::select(-References)

data_traits <- bind_rows(data_traits_ini,species_new) %>% distinct()

#CORRECT SPECIES NAMES
taxon_corr = data_traits %>%
  select(c(Genus, Species)) %>% 
  split(seq(nrow(.))) %>% 
  purrr::map_df(function(x){
    res <- Taxonstand::TPL(genus = x$Genus, species = x$Species, corr = TRUE)
    res %>% select(Family = Family, Genus = New.Genus, Species = New.Species)
  }) 

#SUMMARISE SPECIES VALUES
data_traits_cured <- data_traits %>% 
  ungroup() %>% 
  select(-c(Genus, Species)) %>% 
  cbind(taxon_corr) %>% 
  mutate(taxon = paste(Genus, Species,sep="_"),
         Family = case_when(Genus == "Capparis"~"Capparaceae",
                            Genus == "Dussia"~"Fabaceae",
                            Genus == "Euclea"~"Ebenaceae",
                            Genus == "Gonocaryum"~"Cardiopteridaceae",
                            Genus == "Hibiscus"~"Malvaceae",
                            Genus == "Horsfieldia"~"Myristicaceae",
                            Genus == "Nectandra"~"Lauraceae",
                            Genus == "Platanus"~"Platanaceae",
                            Genus == "Tinospora"~"Menispermaceae",
                            Genus == "Cornus"~"Cornaceae",
                            Genus == "Nyssa"~'Nyssaceae',
                            # Genus == "Baccaurea"~'Phyllanthaceae',
                            # Genus == "Neea"~'Nyctaginaceae',
                            # Genus == "Phyllanthus"~'Phyllanthaceae',
                            # Genus == "Saccharum"~'Poaceae',
                            TRUE ~ Family),
         Family = as_factor(Family),
         Family = fct_recode(Family,
                             Petiveriaceae = 'Phytolaccaceae',
                             Asteraceae = "Compositae",
                             Fabaceae = "Leguminosae"))  %>% 
  # group_by(taxon, Genus, Species, Family) %>%
  # summarise_all(.funs=mean, na.rm = TRUE) %>%
  mutate(across(where(is.numeric),~ifelse(is.nan(.), NA, .)),
         ID = row_number()) %>%
  ungroup() %>% 
  filter(taxon != "NA_na", !is.na(Genus)) %>% 
  arrange(taxon)


#GET TREE
df <- data_traits_cured %>% select(taxon,Genus,Family)
df_tree <- phylo.maker(df)
df_tree$scenario.3$node.label <- NULL


#TRAIT IMPUTATION
# TrEvol::initializeTrEvo()
# res <- imputeTraits(dataset = data.frame(data_traits_cured) %>% select(-Family),
#              phylogeny = df_tree$scenario.3,
#              imputationVariables = c('Height.max..m.', 'P50..MPa.', 'P12..MPa.', 'P88..MPa.', 
#                                      'Slope' 
#                                      # 'Ks..kg.m.1.MPa.1.s.1.','KL..kg.m.1.MPa.1.s.1.', 
#                                      # 'Huber.value','SLA..cm2.g.1.'
#                                      ),
#              predictors = c("MATbest","PPTbest"),
#              prodNAs = 0.0
#             )
# 
# 
# save(res, file = "results/imputation_result.RData")
load(file = "results/imputation_result.RData")

#EXPORT SPECIES

imputation_df <- data_traits_cured %>% 
  dplyr::select(ID,Binomial,Genus, Species, Height.actual..m.,
                Latitude,Longitude, Altitude..m.,MATbest,PPTbest) %>% 
  bind_cols(res$round3$ximp) %>% 
  filter(ID>=1598)

write_csv(imputation_df, file = "imputation_df.csv")



