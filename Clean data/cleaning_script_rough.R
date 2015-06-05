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
col.subset <- c("SITE_ID","VISIT_NO","SAMPLE_CATEGORY","GENUS","TAXANAME","abund_ml","MESH_SIZE", "T_GROUP", "TAXATYPE")

## ZOOPLANKTON
zoo$abund_samp <- (zoo$ABUND/zoo$VOL_COUNT)*zoo$INIT_VOL
zoo$d <- pi*(0.0635)^2*zoo$DEPTH_OF_TOW
zoo$abund_ml <- (zoo$abund_samp/zoo$d)
zoo$TAXATYPE <- "zooplankton"
zoo$T_GROUP <- rep('Zooplankton', nrow(zoo))
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
diatom <- diatom[,col.subset]

# Combining condensed dataframes
full <- rbind(zoo, phyto, diatom)


# Subsetting only first visit to each site
full <- full[full$VISIT_NO==1,]

# Using only the SAMPLE_CATEGORY=="P"
full <- full[full$SAMPLE_CATEGORY=="P",]

# Removing rows with NA in SITE_ID and abund_ml, and abund_ml==0
full <- full[!is.na(full$abund_ml),]
full <- full[full$abund_ml>0,]

## remove 243 form MESH_SIZE
full<-full[-which(full$MESH_SIZE==243),]

## Removing empty TAXANAME and TAXATYPE
full<-full[-which(full$TAXANAME==""),]
full<-full[-which(full$TAXATYPE==""),]

dim(full)# 57177 rows


##  Merging PTL and PTL
id_col = paste0(lakewater$SITE_ID, lakewater$VISIT_NO, lakewater$SAMPLE_CATEGORY)
## Matching PTL and NTL form lakewater to data
aux  = sapply(paste0(full$SITE_ID, full$VISIT_NO, full$SAMPLE_CATEGORY),
    function(r){
        aux = which(r==id_col)
        c(NTL = unique(lakewater[aux,'NTL'])[1],
          PTL = unique(lakewater[aux,'PTL'])[1]) })

full = cbind(full, t(aux))

#### adding lake origin and lon/lat
## LAKE_ORIGIN<-LON_DD<-LAT_DD<-numeric(nrow(full))
## lake.mat<-matrix(nrow=nrow(full),ncol=4)
## lake.mat<-as.data.frame(lake.mat)
## lake.mat[,c(1:4)]<-c(full$SITE_ID,LAKE_ORIGIN,LON_DD,LAT_DD)
## colnames(lake.mat)<-c("SITE_ID","LAKE_ORIGIN","LON_DD","LAT_DD")

## lake2<-lake[,c("SITE_ID","LAKE_ORIGIN","LON_DD","LAT_DD")]
## for(i in 1:nrow(lake.mat)){
##   ind<-which(lake$SITE_ID==lake.mat[i,1])[1]
##   lake.mat[i,c(2,3,4)]<-lake2[ind,c(2,3,4)]
## }
## new code to add lake lon/lat
aux <- sapply(full$SITE_ID, function(r)  which(r==lake$SITE_ID)[1])
full <- cbind(full, lake[aux,c('LAKE_ORIGIN', 'LON_DD', 'LAT_DD')])


## Writing the file
write.csv(full, "full_rough.csv", row.names=FALSE)

## removing all variables except 'full'
rm(list= ls()[ls()!='full'])
