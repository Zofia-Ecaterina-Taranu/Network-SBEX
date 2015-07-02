##  This script joints the data from four sources into one data set and writes it as a csv file##

# This script builds on "cleaning_script_rough.R", adding in "diatom_clean_script.R" 

# Need: Clean data, Network metrics, ways to fit alternative responses (linear vs. stable state)

## Loading libraries
library(reshape2)
library(data.table) 
library(plyr)

##  Data sources (2007 NATIONAL LAKE ASSESSMENT - US EPA) # Note: Can read in using either of 2 options below. 
phyto <- read.csv("../Raw data from NLA/CONTEMPPHYTO.csv", as.is=T) #Soft phytoplankton water-column data 
#phyto <- read.csv(file.choose(), as.is=T) #CONTEMPPHTO.CSV
diatom <- read.csv("../Raw data from NLA/CONTEMPDIATOMS.csv", as.is=T) #Diatom data counted from slides (water-column)
#diatom<- read.csv(file.choose(), as.is=T) #CONTEMPDIATOMS.CSV
zoo <- read.csv("../Raw data from NLA/CONTEMPZOO.csv", as.is=T) #Zooplankton water-column data 
#zoo<- read.csv(file.choose(), as.is=T) #CONTEMPZOO.CSV
lakewater = read.csv('../Raw data from NLA/LAKEWATERQUAL.csv', as.is=TRUE) #Lake water chemistry data 
#lakewater<- read.csv(file.choose(), as.is=T) #LAKEWATERQUAL.CSV
lake <- read.csv("../Raw data from NLA/LAKEINFO.csv", as.is=T) #Lake info data (site locations etc)
#lake<- read.csv(file.choose(), as.is=T) #LAKEINFO.CSV 

## Columns of final data set
col.subset <- c("SITE_ID","VISIT_NO","SAMPLE_CATEGORY","GENUS","TAXANAME","abund_ml", "T_GROUP", "TAXATYPE","BIOVOLUME")


## DATA MANIPULATION FOR PHYTOPLANKTON, DIATOM AND ZOOPLANKTON DATA 
## ZOOPLANKTON
# (note: modified June 22nd 2015 to add in nauplii and fix issue with duplicate taxanames per date)
zoo$NAUPLII[zoo$NAUPLII == "Y"] <- "Naup" #If yes for Nauplii call Naup
zoo$NAUPLII[zoo$NAUPLII == "N"]<- "Adult" #If no for Nauplii call Adult 
zoo$TAXANAME <- paste(zoo$TAXANAME, zoo$NAUPLII, sep="_") # add adult vs nauplii specification to TAXANAME
zoo$TAXANAME[zoo$TAXANAME == "_"] <- "" # we don't want to sum all _
zoo$TAXANAME[zoo$TAXANAME == "_Adult"] <- "" # we don't want to sum all Adults
zoo$TAXANAME[zoo$TAXANAME == "_Naup"] <- "" # we don't want to sum all Nauplii
zoo <- zoo[-which(zoo$TAXANAME==""),] # remove all blanks, don't want to sum
zoo$abund_samp <- (zoo$ABUND/zoo$VOL_COUNT)*zoo$INIT_VOL
zoo$d <- pi*(0.0635)^2*zoo$DEPTH_OF_TOW
zoo$abund_ml <- (zoo$abund_samp/zoo$d)
zoo$TAXATYPE <- "Zooplankton"
zoo$T_GROUP <- paste(zoo$TAXATYPE, zoo$FFG, sep="_") #Add trophic designation for Zooplankton (carnivores etc.)
zoo$T_GROUP[zoo$T_GROUP == "Zooplankton_"] <- "Zooplankton"
zoo$BIOVOLUME <- rep(NA, nrow(zoo))
#zoo <- zoo[-which(zoo$MESH_SIZE==80),] # remove 243 fromm MESH_SIZE (need to do this now because of dcast below - ie. do not want to lump different mesh sizes)
zoo <- zoo[zoo$abund_ml!=0,]
zoo <- zoo[!is.na(zoo$abund_ml),]
zoo <- zoo[,col.subset]
#NB: Originally we had removed mesh size 243 because we thought there was overlap with those organisms (see line above)
#found in mesh size 80. But there actually does not appear to be overlap.
#If we remove mesh size 243, we remove all Cladocerans and Copepods. 
#If we remove mesh size 80, we remove all rotifers, so my suggestion is to keep both mesh sizes. 

# sum abund_ml by TAXANAME to revome duplicates per site 
zoo_cast <- dcast(zoo, SITE_ID ~ TAXANAME, value.var="abund_ml", fun=sum,fill=0)
zoo_melt <- melt(zoo_cast, id.vars = c("SITE_ID"), measured.vas = "abund_ml",
                       variable.name = "TAXANAME", value.name = "abund_ml")
zoo_melt <- zoo_melt[zoo_melt$abund_ml!=0,]

# Add back unique values from remaining columns of zoo (ie. "VISIT_NO", "SAMPLE_CATEGORY","GENUS",
# "T_GROUP", "TAXATYPE","BIOVOLUME") that correspond to unique "SITE_ID" & "TAXANAME" of zoo_melt
zoo_melt$TAXA_SITE <- paste(zoo_melt$TAXANAME, zoo_melt$SITE_ID, sep="_")

zoo$TAXA_SITE <- paste(zoo$TAXANAME, zoo$SITE_ID, sep="_")
col.subset2 <- c("VISIT_NO","SAMPLE_CATEGORY","GENUS","T_GROUP","TAXATYPE","BIOVOLUME","TAXA_SITE")
zoo <- zoo[,col.subset2]
zoo <- zoo[!duplicated(zoo$TAXA_SITE),] # now zoo has the same rows as zoo_melt (14798 obs.)

# Merging the two
zoo <- data.table(zoo, key = "TAXA_SITE")
zoo_melt <- data.table(zoo_melt, key = "TAXA_SITE")
zoo_merged <- zoo_melt[zoo]
zoo_merged <- as.data.frame(zoo_merged) # to removed sorting by key = "TAXA_SITE"
zoo <- zoo_merged[,col.subset]
# Checked: No longer duplicates for taxa within each site and summing correctly. 


## PHYTOPLANKTON
## remove diatoms- because they are also counted in their own samples. 
phyto <- phyto[-grep("diatom",phyto$TAXATYPE,ignore.case=TRUE),]
phyto$T_GROUP<-"Phytoplankton"
names(phyto)[which(names(phyto)=='ABUND')]<-'abund_ml'
phyto <- phyto[,col.subset]
# Checked: all diatoms removed. 


## DIATOM
diatom$abund_ml <- diatom$COUNT/diatom$CONC_VOL
diatom$TAXATYPE <- "Diatoms"
diatom$T_GROUP<-rep("Phytoplankton",nrow(diatom))
diatom$BIOVOLUME<-rep(NA,nrow(diatom))
diatom <- diatom[,col.subset]

#Incorporate diatom biovolumes (From Julien's June 19, 2015 script "diatom_clean_script.R")- use either below to
#read in diatom biovolume file. NB: Amanda made some modifications after copying Julien's script- to make work
#with already condensed diatom data. 
diatom_biovols <- read.csv("../Clean data/NLAdiatoms_biovolumes_CD.csv", header=T, sep=",")
# diatom_biovols<- read.csv(file.choose(), as.is=T) #NLAdiatoms_biovolumes_CD.csv 

# Create diatom dataset
taxanames <- diatom$TAXANAME

## Standardize all taxanames
# Parsing functions
parse_tax_cf <- function(taxaname) {
  # Removes "cf." from taxaname
  taxaname <- gsub(" cf. ", " ", taxaname)
  taxaname <- gsub(" cff. ", " ", taxaname)
  taxaname
}
parse_tax_sp <- function(taxaname) {
  # Combines sp, spp, sp1, etc. variants
  taxaname <- gsub(" sp. [0-9]", " sp.", taxaname)
  taxaname <- gsub(" spp.", " sp.", taxaname)
  taxaname <- gsub(" sp.[0-9]", " sp.", taxaname)
  taxaname <- gsub("[?]", "", taxaname)
  taxaname
}
# Parse taxaname function
parse_tax <- function(taxaname) {
  # Combines parsing sub-functions
  taxaname <- parse_tax_cf(taxaname)
  taxaname <- parse_tax_sp(taxaname)
  taxaname
}
# Parse all taxanames
taxanames_new <- unlist(sapply(taxanames, FUN = parse_tax))

# Split and remake taxanames (removing non species epithets)
split_join_tax <- function(taxaname) {
  split_string <- strsplit(as.character(taxaname), " ")[[1]]
  taxaname <- paste(split_string[1], split_string[2], sep=".")
  taxaname
}
# Split and join taxanames
taxanames_new <- unlist(sapply(taxanames_new, FUN = split_join_tax))

# Replace diatom taxanames with new taxanames
diatom$TAXANAME <- taxanames_new

## Combine diatom row entries by uniqueness
# Subset data.frame and sum abund_ml
diatom_unique <- ddply(diatom, .(SITE_ID,VISIT_NO,SAMPLE_CATEGORY,GENUS,TAXANAME,T_GROUP,TAXATYPE,BIOVOLUME), colwise(sum, c("abund_ml")))
# Match biovolumes from NLAdiatoms_biovolumes_CD.csv
diatom_biovols_index <- match(diatom_unique$TAXANAME, diatom_biovols$species)
unit_biovols <- unlist(sapply(diatom_biovols_index, FUN=function(x) { diatom_biovols$biovolume[x] } ))
# Multiply unit-level biovolumes by abundance
BIOVOLUME <- unit_biovols * diatom_unique$abund_ml
# Add BIOVOLUME to unique diatoms
diatom_unique$BIOVOLUME <- BIOVOLUME

# Reorder diatom data columns
diatom_unique <- diatom_unique[colnames(diatom)]

# Combine data into new full dataset
diatom_new <- diatom[diatom$TAXATYPE!='Diatoms',]
diatom_new <- rbind(diatom_new, diatom_unique)
# Checked: Diatom biovolumes are now added in and names are fixed. 
#NOTE: THERE ARE STILL SOME MISSING DIATOM BIOVOLUMES. 

# Combining condensed dataframes
full <- rbind(zoo, phyto, diatom_new) #uses fixed diatom data. 

# Subsetting only first visit to each site
full <- full[full$VISIT_NO==1,]

# Using only the SAMPLE_CATEGORY=="P"
full <- full[full$SAMPLE_CATEGORY=="P",]

# Removing rows with NA in SITE_ID and abund_ml, and abund_ml==0
full <- full[!is.na(full$abund_ml),]
# full <- full[full$abund_ml>0,]  # Note (from Zo & Amanda): we should not delete abund_ml==0

## Removing empty TAXANAME and TAXATYPE
#full<-full[-which(full$TAXANAME==""),] #Zo, when I run this line, all data removed. 
# --> Zo: looks like it because there are no entries where TAXANAME == ""; see dim(full[full$TAXANAME=="",]) 
full<-full[-which(full$TAXATYPE==""),] #Works 

dim(full)# 49160 rows
# --> Zo: I had 48755 rows

##  Merging PTL and NTL
id_col = paste0(lakewater$SITE_ID, lakewater$VISIT_NO, lakewater$SAMPLE_CATEGORY)
## Matching PTL and NTL form lakewater to data
aux  = sapply(paste0(full$SITE_ID, full$VISIT_NO, full$SAMPLE_CATEGORY),
              function(r){
                aux = which(r==id_col)
                c(NTL = unique(lakewater[aux,'NTL'])[1],
                  PTL = unique(lakewater[aux,'PTL'])[1]) })

full = cbind(full, t(aux))

## new code to add lake lon/lat
aux <- sapply(full$SITE_ID, function(r)  which(r==lake$SITE_ID)[1])
full <- cbind(full, lake[aux,c('LAKE_ORIGIN', 'LON_DD', 'LAT_DD')])

## Converting abund/ml data into biomass data
full$biomass_ind <- rep(NA,nrow(full))

## Applying conversion factors for specific zooplankton taxa to adults only
# (note: THIS IS WHERE WE FILL IN OTHER CONVERSION FACTORS AS WE GET THEM FROM THE LIT)
full$biomass_ind[full$TAXANAME=="Ascomorpha_Adult"] <- 0.0177 # 6 Adults entries, 58 Nauplii
full$biomass_ind[full$TAXANAME=="Asplanchna_Adult"] <- 0.57 # 17 Adults entries, 527 Nauplii
full$biomass_ind[full$TAXANAME=="Bosmina_Adult"] <- 1.577 # 444 Adults, 8 Nauplii
full$biomass_ind[full$TAXANAME=="Brachionus angularis_Adult"] <- 0.029 # 3 Adults, 388 Nauplii
full$biomass_ind[full$TAXANAME=="Calanoida_Adult"] <- 4.2075 # 28 Adults, 1 Nauplii
full$biomass_ind[full$TAXANAME=="Calanoida_Naup"] <- 0.252 
full$biomass_ind[full$TAXANAME=="Collotheca_Adult"] <- 0.0001 # 25 Adults, 138 Nauplii
# full$biomass_ind[full$TAXANAME=="Collothecidae_Adult"] <- 0.0001 # # 0 Adults, 3 Nauplii entries
full$biomass_ind[full$TAXANAME=="Ceriodaphnia_Adult"] <- 0.7107 # 457 Adults, 4 Nauplii
full$biomass_ind[full$TAXANAME=="Chydoridae_Adult"] <- 1.0742 # 216 Adults, 0 Nauplii
full$biomass_ind[full$TAXANAME=="Cyclopidae_Adult"] <- 2.7266 # 1000 Adults, 9 Nauplii
full$biomass_ind[full$TAXANAME=="Cyclopidae_Naup"] <- 0.241 
full$biomass_ind[full$TAXANAME=="Daphnia ambigua_Adult"] <- 3.332 # 72 Adults, 2 Nauplii
full$biomass_ind[full$TAXANAME=="Daphnia mendotae complex_Adult"] <- 8.1962 # 372 Adults, 1 Nauplii
full$biomass_ind[full$TAXANAME=="Daphnia retrocurva_Adult"] <- 2.599 # 111 Adults, 0 Nauplii
full$biomass_ind[full$TAXANAME=="Diaphanosoma_Adult"] <- 2.8804 # 599 Adults, 5 Nauplii
full$biomass_ind[full$TAXANAME=="Diaptomidae_Adult"] <- 3.311 # 944 Adults, 7 Nauplii
full$biomass_ind[full$TAXANAME=="Euchlanis_Adult"] <- 0.1641 # 2 Adults, 62 Nauplii
full$biomass_ind[full$TAXANAME=="Filinia_Adult"] <- 0.0235 # 12 Adults, 339 Nauplii
full$biomass_ind[full$TAXANAME=="Gastropus_Adult"] <- 0.00956 # 20 Adults, 138 Nauplii
# full$biomass_ind[full$TAXANAME=="Kellicottia_Adult"] <- 0.0044 # 0 Adults, 1 Nauplii entry
full$biomass_ind[full$TAXANAME=="Kellicottia bostoniensis_Adult"] <- 0.0044 # 14 Adults, 286 Nauplii
full$biomass_ind[full$TAXANAME=="Kellicottia longispina_Adult"] <- 0.0044 # 28 Adults, 348 Nauplii
full$biomass_ind[full$TAXANAME=="Keratella hiemalis_Adult"] <- 0.0375 # 2 Adults, 42 Nauplii
full$biomass_ind[full$TAXANAME=="Keratella quadrata_Adult"] <- 0.074 # 14 Adults, 258 Nauplii
full$biomass_ind[full$TAXANAME=="Notholca_Adult"] <- 0.0061 # 3 Adults, 14 Nauplii
full$biomass_ind[full$TAXANAME=="Ploesoma_Adult"] <- 0.0225 # 11 Adults, 148 Nauplii
full$biomass_ind[full$TAXANAME=="Polyarthra_Adult"] <- 0.0378 # 33 Adults, 695 Nauplii
full$biomass_ind[full$TAXANAME=="Pompholyx_Adult"] <- 0.0209 # 11 Adults, 187 Nauplii
full$biomass_ind[full$TAXANAME=="Synchaeta_Adult"] <- 0.0409 # 13 Adults, 367 Nauplii

## Converted phytoplankton biovolumes to dry biomass (assuming a specific gravity of 1) and a 
# dry mass: wet mass ratio of 0.10
full$BIOMASS <- 0.1*full$BIOVOLUME # (note: only phytoplankton and diatoms have entries that are not NA, so this conversion does not affect zooplankton)

## Convert phytoplankton biovolumes in μg/L to µg/ml (ie. divide by 1000)
full$BIOMASS <-  0.001*full$BIOMASS

## Calculating Zooplankton biomass (note "Zooplankton" now has a capital "Z")
# (note: according to National Lake Assessment metadata the Phytoplankton BIOVOLUME = CELL_VOLUME * ABUND 
# where CELL_VOLUME is the taxa specific biovolume (µm^3 cell/mL water) and ABUND = COUNT/mL
# Therefore, for zooplankton we took the biomass/indiv. (µg/ indiv. count) x the abund_ml (count/ mL) 
# to obtain total biomass in µg/ mL)
full$BIOMASS[full$TAXATYPE=="Zooplankton"] <- full$biomass_ind[full$TAXATYPE=="Zooplankton"] * full$abund_ml[full$TAXATYPE=="Zooplankton"]

col.subset <- c("SITE_ID","VISIT_NO","SAMPLE_CATEGORY","GENUS","TAXANAME","abund_ml", "T_GROUP", "TAXATYPE", 
                "NTL", "PTL", "LAKE_ORIGIN", "LON_DD", "LAT_DD", "BIOMASS","BIOVOLUME") #Removed Mesh size column
full <- full[,col.subset]

## Writing the file
write.csv(full, "full_combined.csv", row.names=FALSE)

## removing all variables except 'full'
rm(list= ls()[ls()!='full'])

