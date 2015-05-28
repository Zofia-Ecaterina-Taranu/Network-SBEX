
# Diversity analyses

# SBEX NETWORKS

# Max Farrell
options(help_type = "html") # open help in browser
require(reshape2)


full <- read.csv("../Clean data/full_rough.csv", as.is=T)
str(full)


com <- dcast(full, TAXANAME ~ SITE_ID, value.var= 'abund_ml', fun.aggregate=mean,fill=0)
range(com)
aux = paste(full$TAXANAME, full$SITE_ID)
?dcast

length(aux)
length(aux)  - length(unique(aux))
aux[1:10]
dup = duplicated(aux)

length(dup)

full[1:10,]
str(full)


ntw_oligo <- makenetwork_pa(data_oligo_sp_site)
ntw_oligo_abund <- makenetwork_abund(data_oligo_sp_site)


ntw_oligo_abund <- makenetwork_abund(data_oligo_sp_site, method='jaccard')
ntw_oligo_abund <- makenetwork_abund(data_oligo_sp_site)

function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="jaccard"){

ntw_oligo_abund <- makenetwork_pa(data_oligo_sp_site, method='bray')
ntw_oligo_abund <- makenetwork_abund(data_oligo_sp_site,method = 'bray')
