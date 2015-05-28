
# Diversity analyses

# SBEX NETWORKS

# Max Farrell

require(reshape2)
require(vegan)
options(help_type ='html')
full <- read.csv("../Clean data/full_rough.csv",header= TRUE, as.is=T)
full$MESH_SIZE[is.na(full$MESH_SIZE)]<-1
full = full[grep("(1|243)", full$MESH_SIZE),]
unique(full$MESH_SIZE)
str(full)

# Abundances across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
abund  <- dcast(full, SITE_ID ~ TAXANAME, value.var="abund_ml", fun=mean, fill=0)
abund = abund[,-1]


div = diversity(abund,   index = 'shannon',MARGIN=1 )
# PTL across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
index = data.frame(unique(cbind(SITE_ID = full$SITE_ID,PTL= full$PTL, NTL = full$NTL)))


plot(as.double(index$NTL), div, type='p')
plot(as.double(index$PTL), div, type='p')
library(copula)
library(QRM)
ntl_ptl <- cbind(as.double(index$PTL), as.double(index$NTL))
plot(pobs(ntl_ptl))
plot(as.double(index$PTL), as.double(index$NTL))
log(div)



##  Zooplankton
full_zoo = subset(full, TAXANAME ='zooplankton')

# Abundances across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
abund  <- dcast(full_zoo, SITE_ID ~ TAXANAME, value.var="abund_ml", fun=mean, fill=0)
abund = abund[,-1]

div = diversity(abund,   index = 'shannon',MARGIN=1 , base=10)
# PTL across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
index = data.frame(unique(cbind(SITE_ID = full_zoo$SITE_ID,PTL= full_zoo$PTL, NTL = full_zoo$NTL)))


plot(as.double(index$NTL), div, type='p')
plot(as.double(index$PTL), div, type='p')


##  Diatom
full_dia = subset(full, TAXANAME ='diatom')
# Abundances across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
abund  <- dcast(full_dia, SITE_ID ~ TAXANAME, value.var="abund_ml", fun=mean, fill=0)
abund = abund[,-1]

div = diversity(abund,   index = 'shannon',MARGIN=1 )
# PTL across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
index = data.frame(unique(cbind(SITE_ID = full_dia$SITE_ID,PTL= full_dia$PTL, NTL = full_dia$NTL)))


plot(as.double(index$NTL), div, type='p')
plot(as.double(index$PTL), div, type='p')


str(full)

#Script for diversity_analyses for Evenness (i.e., H’).  This is not calculated for each functional group.  I have to read Krebs (1989)'s description of H': some researchers calculate log(n-1).  However, Vegan writes log(specnumber). 
rowSums(abund)
log(rowSums(abund))
specnumber<-log(rowSums(abund))
Evenness<-div/specnumber
plot(log(index$PTL),Evenness)
plot(log(index$NTL),Evenness)

#When I tried to run the diversity index and evenness for phytoplankton in full_rough, I found that TAXATYPE only had zooplankton and diatom.  Hercules started working to solve this in full_rough.


#ZOOPLANKTON ONLY DIVERSITY AND EVENNESS
#Interaction matrix between site and species for zooplankton only (needs to be verified!)
full.zoo<- subset(full, TAXATYPE == "zooplankton", drop=T)
unique(full.zoo$TAXATYPE)
abund.zoo<-dcast(full.zoo, SITE_ID ~ TAXANAME, value.var="abund_ml", fun=mean, fill=0)
abund.zoo = abund.zoo[,-1]
View(abund.zoo)

div.zoo = diversity(abund.zoo,   index = 'shannon',MARGIN=1 )
# PTL across non-unique site-taxa (zooplankton only) combinations were SUMMED **** DOUBLE CHECK THIS.
index.zoo = data.frame(unique(cbind(SITE_ID = full.zoo$SITE_ID,PTL= full.zoo$PTL, NTL = full.zoo$NTL)))


par(mfrow=c(2,2))
plot(index.zoo$NTL, div.zoo, type='p')
plot(index.zoo$PTL, div.zoo, type='p')
rowSums(abund.zoo)
log(rowSums(abund.zoo))
specnumber.zoo<-log(rowSums(abund.zoo))
Evenness.zoo<-div.zoo/specnumber.zoo
plot(log(index.zoo$PTL),Evenness.zoo)
plot(log(index.zoo$NTL),Evenness.zoo)
abline(lm(Evenness.zoo~log(index.zoo$NTL)),col="red")
lm(Evenness.zoo~log(index.zoo$NTL))
summary(lm(Evenness.zoo~log(index.zoo$NTL)))#R_squared=0.07.  p-values are sig, but that is to be expected for large data sets.

simp.zoo <- diversity(abund.zoo, MARGIN=1, "simpson")
richness.zoo <- specnumber(abund.zoo, MARGIN=1)
par(mfrow=c(2,2))
plot(index.zoo$NTL, simp.zoo, type='p')
plot(index.zoo$PTL, simp.zoo, type='p')
plot(index.zoo$NTL, richness.zoo, type='p')
plot(index.zoo$PTL, richness.zoo, type='p')

#Continue with diversity indices and evenness for Phytoplankton, Phytoplankton functional groups, and diatom.




