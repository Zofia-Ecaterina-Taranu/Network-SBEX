
# Diversity analyses

# SBEX NETWORKS

# Max Farrell

require(reshape2)
require(vegan)
require(plyr)


full <- read.csv("../Clean data/full_rough.csv", as.is=T)

str(full)

# Abundances across non-unique site-taxa combinations were SUMMED **** DOUBLE CHECK THIS.
com <- dcast(full, TAXANAME ~ SITE_ID, value.var="abund_ml", fun=mean, fill=0)
rownames(com) <- com$TAXANAME
com <- as.matrix(t(com[,-1]))


# Diversity metrics

# WHOLE COMMUNITY
h <- diversity(com, MARGIN=1, "shannon")
simp <- diversity(com, MARGIN=1, "simpson")
richness <- specnumber(com, MARGIN=1)
zoosp <- sort(unique(full$TAXANAME[full$TAXATYPE=="zooplankton"]))
diatomsp <- sort(unique(full$TAXANAME[full$TAXATYPE=="diatom"]))

div <- data.frame(rownames(com),h,simp,richness)
names(div)[1] <- "SITE_ID"

div <- join(div,full[,c("SITE_ID","PTL","NTL")], type="left")

div$PTL <- log(div$PTL)
div$NTL <- log(div$NTL)

pairs(div[,-1])


# PHYTO
phytosp <- sort(unique(full$TAXANAME[full$TAXATYPE!="diatom" & full$TAXATYPE!="zooplankton"]))
com.phyto <- com[,colnames(com)%in%phytosp]

h <- diversity(com.phyto, MARGIN=1, "shannon")
simp <- diversity(com.phyto, MARGIN=1, "simpson")
richness <- specnumber(com.phyto, MARGIN=1)

div <- data.frame(rownames(com.phyto),h,simp,richness)
names(div)[1] <- "SITE_ID"

div <- join(div,full[,c("SITE_ID","PTL","NTL")], type="left")

div$PTL <- log(div$PTL)
div$NTL <- log(div$NTL)

pairs(div[,-1])
#Hi Max.  This is Sacha.  Heads up! There are no Phytoplankton in the TAXATYPE in the full_rough file.  Hercules is working on this.

