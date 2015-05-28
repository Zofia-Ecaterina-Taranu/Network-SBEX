#************************************************************************************************#
#                                    NLA network analysis                                   ####
#************************************************************************************************#
options(help_type = "html")

library(reshape2)
library(igraph)
# Loading data ####
data <- read.csv("../Clean data/full_rough.csv", as.is=TRUE)

# All sites
data_sp_site <- dcast(data,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
re_data_sp_site = t(apply(data_sp_site[,-1],1,  function(r) r /sum(r)))
groups = sapply(colnames(data_sp_site)[-1], function(r)  data$TAXATYPE[which(r == data$TAXANAME)[1]] )
groups = as.numeric(factor(groups))


makenetwork <- function(spe, threshold, groups, plot= FALSE){
    ##  Makeing network per site based on relative abundance.
    if(missing(threshold)) threshold =0

    spe = spe %*% t(spe)
    ## set diagonal to 0 (self-interaction)
    diag(spe) <- 0 
    ## spe = 1*(spe>threshold)
    ## Creating the graph
    ig = graph.adjacency(spe, mode= "undirected", weighted = TRUE)
    if(plot)
        plot(ig, layout=layout.fruchterman.reingold,
             vertex.size=0.6, vertex.label.dist=0.1,
             edge.arrow.mode="-",vertex.color="red",
             vertex.label=NA,edge.color="blue")
    ## title("Abund dist")
    list(mod = modularity(x = ig, membership=groups, weights =E(ig)$weight),
         no.edge = ecount(ig), 
         e.density = ecount(ig)/(ecount(ig)*(ecount(ig)-1)/2)*100, # realized/ total
         # avg.path = average.path.length(ig),
         no.clust = no.clusters(ig),
         max.clust.size = max(clusters(ig)$csize),
         avg.clust.size = mean(clusters(ig)$csize))
         # transit = transitivity(ig),
         # connectance = vcount(ig)/(ecount(ig))^2, # from Coll et al. (2011) Table 2
         # l.density = vcount(ig)/ecount(ig), # from Coll et al. (2011) Table 2 
         
}

makenetwork(re_data_sp_site[1,], groups=groups,plot=T)
stats = apply(re_data_sp_site, 1, function(r,g) makenetwork(r, groups=g), g= groups)
stats = do.call(rbind, stats)

NTW = unique(data[,c('SITE_ID', 'PTL', 'NTL','LAKE_ORIGIN')])
NTW$NTL_PTL = NTW$NTL/NTW$PTL

aux = sapply(NTW[,'SITE_ID'], function(r) which(r==data_sp_site[,'SITE_ID']))
NTW = cbind(NTW, mod = stats[aux, ])
colnames(NTW)
NTW$LAKE_ORIGIN <- as.factor(NTW$LAKE_ORIGIN)

par(mfrow=c(1,3))
# Modularity vs Phosphorus, Nitrogen and TN:TP ratio
# The modularity can be either positive or negative, with positive values indicating the 
# possible presence of community structure. Negative values indicate disassortative mixing
plot(log(NTW[,'PTL']), (NTW[,'mod.mod']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL']), (NTW[,'mod.mod']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL_PTL']), (NTW[,'mod.mod']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))

# Edge number vs Phosphorus, Nitrogen and TN:TP ratio
plot(log(NTW[,'PTL']), (NTW[,'mod.no.edge']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL']), (NTW[,'mod.no.edge']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL_PTL']), (NTW[,'mod.no.edge']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))

# Cluster number vs Phosphorus, Nitrogen and TN:TP ratio
plot(log(NTW[,'PTL']), (NTW[,'mod.no.clust']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL']), (NTW[,'mod.no.clust']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL_PTL']), (NTW[,'mod.no.clust']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))

# Max cluster size vs Phosphorus, Nitrogen and TN:TP ratio
plot(log(NTW[,'PTL']), (NTW[,'mod.max.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL']), (NTW[,'mod.max.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL_PTL']), (NTW[,'mod.max.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))

# Avg cluster size vs Phosphorus, Nitrogen and TN:TP ratio
plot(log(NTW[,'PTL']), (NTW[,'mod.avg.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'PTL']), (NTW[,'mod.avg.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
plot(log(NTW[,'NTL_PTL']), (NTW[,'mod.avg.clust.size']), col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
