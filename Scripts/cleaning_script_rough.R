##  This script joints the data from four sources into one data set and writes it as a csv file.
##

# Need: Clean data, Network metrics, ways to fit alternative responses (linear vs. stable state)

##  Data sources
phyto <- read.csv("../Raw data from NLA/CONTEMPPHYTO.csv", as.is=T)
diatom <- read.csv("../Raw data from NLA/CONTEMPDIATOMS.csv", as.is=T)
zoo <- read.csv("../Raw data from NLA/CONTEMPZOO.csv", as.is=T)
lakewater = read.csv('../Raw data from NLA/LAKEWATERQUAL.csv', as.is=TRUE)
lake <- read.csv("../Raw data from NLA/LAKEINFO.csv", as.is=T)
## Columns of final data set
col.subset <- c("SITE_ID","VISIT_NO","SAMPLE_CATEGORY","GENUS","TAXANAME","abund_ml","MESH_SIZE", "T_GROUP", "TAXATYPE","BIOVOLUME")

## ZOOPLANKTON
zoo$abund_samp <- (zoo$ABUND/zoo$VOL_COUNT)*zoo$INIT_VOL
zoo$d <- pi*(0.0635)^2*zoo$DEPTH_OF_TOW
zoo$abund_ml <- (zoo$abund_samp/zoo$d)
zoo$TAXATYPE <- "zooplankton"
zoo$T_GROUP <- rep('Zooplankton', nrow(zoo))
zoo$BIOVOLUME <- rep(NA, nrow(zoo))
zoo <- zoo[,col.subset]

## PHYTOPLANKTON
## remove diatoms
phyto <- phyto[-grep("diatom",phyto$TAXATYPE,ignore.case=TRUE),]
phyto$MESH_SIZE <- rep(NA, nrow(phyto))
phyto$T_GROUP<-"Phytoplankton"
names(phyto)[which(names(phyto)=='ABUND')]<-'abund_ml'
phyto <- phyto[,col.subset]

## DIATOM
diatom$abund_ml <- diatom$COUNT/diatom$CONC_VOL
diatom$TAXATYPE <- "Diatoms"
diatom$MESH_SIZE <- rep(NA, nrow(diatom))
diatom$T_GROUP<-rep("Phytoplankton",nrow(diatom))
diatom$BIOVOLUME<-rep(NA,nrow(diatom))
diatom <- diatom[,col.subset]

# Combining condensed dataframes
full <- rbind(zoo, phyto, diatom)

# Subsetting only first visit to each site
full <- full[full$VISIT_NO==1,]

# Using only the SAMPLE_CATEGORY=="P"
full <- full[full$SAMPLE_CATEGORY=="P",]

# Removing rows with NA in SITE_ID and abund_ml, and abund_ml==0
full <- full[!is.na(full$abund_ml),]
# full <- full[full$abund_ml>0,]  # Note (from Zo & Amanda): we should not delete abund_ml==0

## remove 243 form MESH_SIZE
full<-full[-which(full$MESH_SIZE==243),]

## Removing empty TAXANAME and TAXATYPE
full<-full[-which(full$TAXANAME==""),]
full<-full[-which(full$TAXATYPE==""),]

dim(full)# 44343 rows

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

## Applying conversion factors for specific zooplankton taxa
# (note: THIS IS WHERE WE FILL IN OTHER CONVERSION FACTORS AS WE GET THEM FROM THE LIT)
full$biomass_ind[full$TAXANAME=="Ascomorpha"] <- 0.0177
full$biomass_ind[full$TAXANAME=="Bosmina"] <- 1.577 
full$biomass_ind[full$TAXANAME=="Calanoida"] <- 4.2075
full$biomass_ind[full$TAXANAME=="Collotheca"] <- 0.0001
full$biomass_ind[full$TAXANAME=="Collothecidae"] <- 0.0001
full$biomass_ind[full$TAXANAME=="Ceriodaphnia"] <- 0.7107
full$biomass_ind[full$TAXANAME=="Chydoridae"] <- 1.0742
full$biomass_ind[full$TAXANAME=="Cyclopidae"] <- 2.7266
full$biomass_ind[full$TAXANAME=="Daphnia ambigua"] <- 3.332
full$biomass_ind[full$TAXANAME=="Daphnia mendotae complex"] <- 8.1962
full$biomass_ind[full$TAXANAME=="Daphnia retrocurva"] <- 2.599
full$biomass_ind[full$TAXANAME=="Diaphanosoma"] <- 2.8804
full$biomass_ind[full$TAXANAME=="Diaptomidae"] <- 3.311
full$biomass_ind[full$TAXANAME=="Filinia"] <- 0.0235
full$biomass_ind[full$TAXANAME=="Gastropus"] <- 0.00956
full$biomass_ind[full$TAXANAME=="Kellicottia"] <- 0.0044
full$biomass_ind[full$TAXANAME=="Kellicottia bostoniensis"] <- 0.0044
full$biomass_ind[full$TAXANAME=="Kellicottia longispina"] <- 0.0044
full$biomass_ind[full$TAXANAME=="Keratella hiemalis"] <- 0.0375
full$biomass_ind[full$TAXANAME=="Keratella quadrata"] <- 0.074
full$biomass_ind[full$TAXANAME=="Notholca"] <- 0.0061
full$biomass_ind[full$TAXANAME=="Ploesoma"] <- 0.0225
full$biomass_ind[full$TAXANAME=="Polyarthra"] <- 0.0378
full$biomass_ind[full$TAXANAME=="Pompholyx"] <- 0.0209
full$biomass_ind[full$TAXANAME=="Synchaeta"] <- 0.0409


## Converted phytoplankton biovolumes to dry biomass (assuming a specific gravity of 1) and a 
# dry mass: wet mass ratio of 0.10
full$BIOMASS <- 0.1*full$BIOVOLUME # (note: so far zooplankton data has NAs, so this conversion does not affect zoops)

## Convert phytoplankton biovolumes in μg/L to µg/ml (ie. divide by 1000)
full$BIOMASS <-  0.001*full$BIOMASS

## Calculating zooplankton biomass
# (note: according to National Lake Assessment metadata the Phytoplankton BIOVOLUME = CELL_VOLUME * ABUND 
# where CELL_VOLUME is the taxa specific biovolume (µm^3 cell/mL water) and ABUND = COUNT/mL
# Therefore, for zooplankton we took the biomass/indiv. (µg/ indiv. count) x the abund_ml (count/ mL) 
# to obtain total biomass in µg/ mL)
full$BIOMASS[full$TAXATYPE=="zooplankton"] <- full$biomass_ind[full$TAXATYPE=="zooplankton"] * full$abund_ml[full$TAXATYPE=="zooplankton"]

col.subset <- c("SITE_ID","VISIT_NO","SAMPLE_CATEGORY","GENUS","TAXANAME","abund_ml","MESH_SIZE", "T_GROUP", "TAXATYPE", 
                "NTL", "PTL", "LAKE_ORIGIN", "LON_DD", "LAT_DD", "BIOMASS","BIOVOLUME")
full <- full[,col.subset]

## Writing the file
write.csv(full, "full_rough.csv", row.names=FALSE)

## removing all variables except 'full'
rm(list= ls()[ls()!='full'])

