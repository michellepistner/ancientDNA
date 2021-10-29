##Loading libraries
##Fido needs to be installed using github
##phyloseq and microbiome need to be installed using BioManager 
library(phyloseq)
library(fido)
library(driver)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringi)
library(stringr)

##Setting the seed
set.seed(2021)

#####Part 1: Reading, Filtering, and Processing Metadata/OTU data#####

###Reading in the metadata
metadata = read.delim(file.path("Data", "BritishData", "All_Metadata_2282021.txt"))

##Filtering metadata to Museum of London only with a value for Date_100
metadata_london = metadata %>%
  filter(str_detect(Museum,"MoL")) %>%
  filter(!is.na(BlackDeath_1346_1353)) 

dim(metadata_london)

###Reading in the OTU data
OTU_data = read.delim(file.path("Data", "BritishData", "AllSamples_20210212_RawAbsolute_AllTaxa.txt"))
dim(OTU_data)

##Light cleaning of the names so the taxanomic information is of the same type for every taxa
otu.names.all = str_sub(OTU_data$X.Datasets, end=-2)
otu.names.all = gsub(";[a-zA-Z /]+ group", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z /]+ incertae sedis", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z /]+ complex", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z /]+ subgroup", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z /]+ subdivisions", "", otu.names.all)
otu.names.all = gsub(";Polyangiaceae", "", otu.names.all)

otu.names.all = stri_list2matrix(str_split(otu.names.all, ";"), byrow=TRUE)
otu.names.all = as.data.frame(otu.names.all)
names(otu.names.all) = c("root", "cellular", "Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")

##Subsetting to the genus level
OTU_data$X.Datasets = otu.names.all$Genus

OTU_data = OTU_data %>%
  filter(!is.na(X.Datasets)) %>%
  group_by(X.Datasets) %>%
  summarise_each(funs(sum)) %>%
  remove_rownames %>% 
  column_to_rownames(var="X.Datasets") 


###Filtering the OTU data to include the samples in the metadata only
###Filter the OTU data
london.names = paste0("X",metadata_london$X.SampleID)
OTU.london = OTU_data[,colnames(OTU_data) %in% c(london.names,metadata_london$X.SampleID)]

##Now, reverse (filtering out metadata to contain only OTU samples) just in case
metadata_london = metadata_london[paste0("X",metadata_london$X.SampleID)%in% colnames(OTU.london), ]

##Checking dimensions
dim(OTU.london)
dim(metadata_london)

##Now, we need to filter out the OTU data.
##First, a quick plot of the rowSums
OTU.london %>% rowSums() %>% ecdf() %>% base::plot() %>% abline(v=1e6)
###Looks like 40% of the samples have counts of zero

###Use a conservative filtering rule
###Filter out all taxa that don't have a count of at least one in 30% of samples
filtered = rowSums(OTU.london > 1) < .25*ncol(OTU.london)
other.tot = colSums(OTU.london[filtered,])
filtered[122] = TRUE
otu.filtered = rbind(OTU.london[!filtered,], "Other" = other.tot + OTU.london[122,])
dim(otu.filtered)

###Reordering metadata_london to match the OTU table
metadata_london = metadata_london[match(colnames(otu.filtered), paste0("X", metadata_london$X.SampleID)),]

latLongs = read.csv(file.path("Data", "BritishData", "cemetryLocations.csv"))

metadata_london$LateDate = ifelse(is.na(metadata_london$LateDate), stri_list2matrix(str_split(metadata_london$Date_100,"_"))[,2],metadata_london$LateDate)

##Creating the final metadata with the needed variables
metadata_prep = metadata_london %>%
  select(Date_100, Date_200, Date_300, BlackDeath_PrePost, EarlyDate, LateDate, MedievalPostMedieval, Cemetry, MaxillaMandible, BuccalLingual, SubSupragingival, Tooth, Tooth_Simplified, BlackDeath_1346_1353)

