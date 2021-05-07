####Analysis of Covariation#########################
####Michelle Nixon##################################
####April 12, 2021##################################

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
metadata = read.csv(file.path("Data", "BritishData", "FullMetadata_20210428.csv"))

##Filtering metadata to Museum of London only with a value for Date_100
metadata_london = metadata %>%
  filter(str_detect(Museum,"MoL")) %>%
  filter(!is.na(BlackDeath_PrePost)) 

dim(metadata_london)

###Reading in the OTU data
OTU_data = read.delim(file.path("Data", "BritishData", "AllSamples_20210212_RawAbsolute_AllTaxa.txt"))
dim(OTU_data)

otu.names.all = str_sub(OTU_data$X.Datasets, end=-2)
otu.names.all = gsub(";[a-zA-Z /]+ group", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z /]+ incertae sedis", "", otu.names.all)
otu.names.all = gsub(";[a-zA-Z / \\.]+ Incertae Sedis", "", otu.names.all)
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
otu.filtered = rbind(OTU.london[!filtered,], "Other" = other.tot)
dim(otu.filtered)

###Reordering metadata_london to match the OTU table
metadata_london = metadata_london[match(colnames(otu.filtered), paste0("X", metadata_london$X.SampleID)),]

latLongs = read.csv(file.path("Data", "BritishData", "cemetryLocations.csv"))

metadata_london$LateDate = ifelse(is.na(metadata_london$LateDate), stri_list2matrix(str_split(metadata_london$Date_100,"_"))[,2],metadata_london$LateDate)

##Creating the final metadata with the needed variables
metadata_prep = metadata_london %>%
  select(Date_100, Date_200, Date_300, BlackDeath_PrePost, EarlyDate, LateDate, MedievalPostMedieval, Cemetry, MaxillaMandible, BuccalLingual, SubSupragingival)



#####Part 2: Selecting the Best Model for Time#####
var.explained.Y <- function(posterior, otu.closed){
  posterior = to_proportions(posterior)
  posterior = to_alr(posterior, ncategories(posterior))
  Y.pred = predict(posterior, response = "LambdaX", from_scratch = FALSE)
  Eta.inv = array(0, dim = c(dim(otu.closed), 2000))
  for(i in 1:dim(Y.pred)[3]){
    Eta.inv[,,i] = alrInv_array(Y.pred[,,i], d = nrow(Y.pred) + 1, 1)
  }
  Ypred = array(0, dim=dim(Eta.inv))
  size = as.data.frame(colSums(posterior$Y))
  for (i in 1:dim(Ypred)[3]){
    for (j in 1:dim(Ypred)[2]){
      Ypred[,j,i] <- rmultinom(1, size[j,] ,Eta.inv[,j,i])## Eta.inv[,j,i] * size[j,]
    }
  }
  exp.mat = matrix(NA, nrow = dim(Ypred)[3], ncol = 1)
  for(i in 1:dim(Y.pred)[3]){
    Y.pred.closed = miniclo_array(Ypred[,,i], parts = 1)
    var_fit = sum(apply(Y.pred.closed, MARGIN = 1, FUN = "var"))
    var_res = sum(apply(Y.pred.closed - otu.closed, MARGIN = 1, FUN = "var"))
    exp.mat[i,1] = var_fit/ (var_fit + var_res)
    
  }
  return(exp.mat)
}#end of function

gamm.select = 1 ##Baseline signal to noise ratio

###Prepping the OTU table
otu.closed = miniclo_array(otu.filtered, part = 1)
Y = otu_table(otu.filtered, taxa_are_rows=TRUE)

###Setting priors based on Y
upsilon = ntaxa(Y)+3
Omega = diag(ntaxa(Y))
G = cbind(diag(ntaxa(Y)-1), -1)
Xi = (upsilon-ntaxa(Y))*G%*%Omega%*%t(G)

###Intercept only model
X = t(model.matrix(~1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot=1)
#posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE) 
ppc_summary(posterior, from_scratch=TRUE)

var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Intercept-Only Model (Counts):")
mean(var.exp.Y)

###Date model
X = t(model.matrix(~BlackDeath_PrePost - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE)
ppc_summary(posterior, from_scratch=TRUE)


var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Date Model (Counts):")
mean(var.exp.Y)

###Cemetery model
X = t(model.matrix(~Cemetry - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.Y = var.explained.Y(posterior,otu.closed)

###Percentage of variance explained

print("Percent of Variation Explained by Cemetery Model (Counts):")
mean(var.exp.Y)


#####Part 3: Cross-Validation to Determine Signal to Noise Ratio#####

###Going to maximize the percentage of variance explained
gamm.opt = seq(1,5, by = 1)
avg.var = rep(NA, length(gamm.opt))

for(i in 1:length(gamm.opt)){
  ###Date model
  X = t(model.matrix(~BlackDeath_PrePost - 1, data = metadata_prep))
  
  ###Setting priors based on X
  Theta = matrix(0, ntaxa(Y)-1, nrow(X))
  Gamma = diag(nrow(X))*gamm.opt[i]
  
  posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
  posterior = to_clr(posterior)
  
  var.exp = var.explained.Y(posterior, otu.closed)
  ###Percentage of variance explained
  avg.var[i] =  mean(var.exp)
  print(i)
}

gamm.select = which.max(avg.var)


####Fitting optimal model
###Date_100 model
X = t(model.matrix(~BlackDeath_PrePost - 1, data = metadata_prep))

###Setting priors based on X
Theta = matrix(0, ntaxa(Y)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior = pibble(Y, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)
###Checking the posterior predictive checks
ppc(posterior) 
ppc_summary(posterior)

ppc(posterior, from_scratch=TRUE)
ppc_summary(posterior, from_scratch=TRUE)


var.exp = var.explained.Y(posterior, otu.closed)
###Percentage of variance explained

print("Percent of Variation Explained by Date_100 Model (Counts):")
mean(var.exp)
quantile(var.exp, c(0.025,0.975))

#####Part 4: Percent of Variability Explained by Other Sources#####


###First, getting the part "explained by" Date_100
Y.pred = predict(posterior, response = "Y", from_scratch = TRUE)

###Variation explained by other sources
Y.samp = apply(Y.pred, MARGIN=c(1,2), FUN = "mean")

otu.closed.Ysamp = miniclo_array(Y.samp, parts = 1)

rownames(Y.samp) = rownames(Y)
names(Y.samp) = names(Y)
Y.samp = otu_table(Y.samp, taxa_are_rows=TRUE)


##Sanity check, should be close to zero
X = t(model.matrix(~1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)


var.exp.int = var.explained.Y(posterior, otu.closed.Ysamp)
###Percentage of variance explained

print("Percent of Variation in Time Explained by Intercept-Only Model (Counts):")

mean(var.exp.int)
quantile(rowMeans(var.exp.int), c(0.025,0.975))

##Another sanity check, should be close to 100%
X = t(model.matrix(~BlackDeath_PrePost - 1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.date = var.explained.Y(posterior, otu.closed.Ysamp)
###Percentage of variance explained

print("Percent of Variation in Time Explained by Date_100 Model (Counts):")
mean(var.exp.date)
quantile((var.exp.date), c(0.025,0.975))

X = t(model.matrix(~Cemetry - 1, data = metadata_prep))

Theta = matrix(0, ntaxa(Y.samp)-1, nrow(X))
Gamma = diag(nrow(X))*gamm.select

posterior <- pibble(Y.samp, X, upsilon, Theta, Gamma, Xi, multDirichletBoot = 1)
posterior = to_clr(posterior)

var.exp.cem = var.explained.Y(posterior, otu.closed.Ysamp)

###Percentage of variance explained

print("Percent of Variation in Time Explained by Cemetry Model (Counts):")
mean(var.exp.cem)
quantile((var.exp.cem), c(0.025,0.975))

