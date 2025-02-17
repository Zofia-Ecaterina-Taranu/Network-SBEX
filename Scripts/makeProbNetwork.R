#************************************************************************************************#
#                                    NLA network analysis                                   ####
#************************************************************************************************#
options(help_type = "html")

library(reshape2)
library(igraph)
library(vegan)
library(ggplot2)

# Loading data ####
full <- read.csv("../Clean data/full_combined.csv", as.is=TRUE)
data <- full[!is.na(full$BIOMASS),] # Removing rows with NAs for BIOMASS

# Note (from Zo): Suggestion to replace TAXATYPE by T_GROUP for Zooplankton to break them down a bit more
data$TAXATYPE[data$TAXATYPE=="Zooplankton"] <- data$T_GROUP[data$TAXATYPE=="Zooplankton"]

# All sites
data_sp_site <- dcast(data,SITE_ID ~ TAXANAME, value.var="BIOMASS",fun=mean,fill=0)
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
  ig <- graph.adjacency(spe, mode= "undirected", weighted = TRUE)
  ly = make_layout(ig, groups)
  deg  = degree(ig)
  deg[deg>0] <- 6
  deg = deg + 2
  if(plot)
    plot(ig, layout=20*ly,
         vertex.size=deg, vertex.label.dist=0.5,
         edge.arrow.mode="-",vertex.color=groups,
         vertex.label=NA, edge.color="blue", edge.width= 50*E(ig)$weight)
  ## title("Abund dist")
   list(mod = modularity(x = ig, membership=groups, weights =E(ig)$weight),
        no.edge = ecount(ig),
        e.density = ecount(ig)/(ecount(ig)*(ecount(ig)-1)/2)*100, # realized/ total
        # avg.path = average.path.length(ig),
        no.clust = no.clusters(ig),
        connectance = vcount(ig)/(ecount(ig))^2, # from Coll et al. (2011) Table 2
        # transit = transitivity(ig),
        # l.density = vcount(ig)/ecount(ig), # from Coll et al. (2011) Table 2
        avg.degree = mean(degree(ig)),
        avg.closeness = mean(closeness(ig)), # vertex is ‘central’ if it is ‘close’ to many other vertices
        avg.betweenness = mean(betweenness(ig)), # extent to which a vertex is located ‘between’ other pairs of vertices; ‘importance’ relates to where a vertex is located with respect to the paths in the network graph
        avg.eigenvector= mean(evcent(ig)$vector), # 'prestige' the more central the neighbors of a vertex are, the more central that vertex itself is
        max.clust.size = max(clusters(ig)$csize),
        avg.clust.size = mean(clusters(ig)$csize))
}

make_layout<-function(ig, groups){
  ## A function that takes an igraph object and group labels and creates a circle
  ## layout.
  ly = layout.circle(ig)
  ly.new = ly[order(order(groups)),]
  ly.new
}

# Plotting by groups ####

# Creating groups
gp =list(g1 = sample(unique(data$SITE_ID[which(data$PTL<=10)]), 4),
         g2 = sample(unique(data$SITE_ID[which(data$PTL>10 & data$PTL <= 30)]), 4),
         g3 = sample(unique(data$SITE_ID[which(data$PTL>30 & data$PTL <= 60)]), 4),
         g4 = sample(unique(data$SITE_ID[which(data$PTL>60)]), 4))
k=1

gp_id = sapply(gp, function(r, g){
  ## r is the sublist g1, g2, etc.
  pdf(paste(k, '.pdf'))
  par(mfrow=c(2,2))
  aux = which(data_sp_site$SITE_ID %in% r)
  for (i in 1:length(aux)){
    makenetwork(re_data_sp_site[aux[i],],groups= g, plot=TRUE)
    tl = data$PTL[which(data$SITE_ID==r[i])[1]]
    title(paste('PTL:', tl))
  }
  dev.off()
  k<<-k+1
}, g=groups)

# Abund vs biomass ####
# (note: looking at the same lakes but using different dataset - ie. looked at all species with abund_ml as the response
# then subset of lakes with just biomass values and looking at either abund_ml or biomass as the response)
# pdf(paste(k, '.pdf'))
# par(mfrow=c(2,2))
# Oligotrophic example
# r="NLA06608-1167"
# r="NLA06608-3890"
# r="NLA06608-NH4912"
# aux = which(data_sp_site$SITE_ID %in% r)
# makenetwork(re_data_sp_site[aux,],groups= groups, plot=TRUE)
# tl <- data$PTL[which(data$SITE_ID==r[1])[1]]
# title(paste('PTL:', tl))

# Mesotrophic example
# r="NLA06608-1056"
# r="NLA06608-1206"
# r="NLA06608-1489"
# aux = which(data_sp_site$SITE_ID %in% r)
# makenetwork(re_data_sp_site[aux,],groups= groups, plot=TRUE)
# tl <- data$PTL[which(data$SITE_ID==r[1])[1]]
# title(paste('PTL:', tl))

# Eutrophic example
# r="NLA06608-2955"
# r="NLA06608-1724"
# r="NLA06608-0219"
# aux = which(data_sp_site$SITE_ID %in% r)
# makenetwork(re_data_sp_site[aux,],groups= groups, plot=TRUE)
# tl <- data$PTL[which(data$SITE_ID==r[1])[1]]
# title(paste('PTL:', tl))

# Hypereutrophic example
# r="NLA06608-1239"
# r="NLA06608-2123"
# r="NLA06608-2831"
# aux = which(data_sp_site$SITE_ID %in% r)
# makenetwork(re_data_sp_site[aux,],groups= groups, plot=TRUE)
# tl <- data$PTL[which(data$SITE_ID==r[1])[1]]
# title(paste('PTL:', tl))
# dev.off()

# Exploratory analysis - network properties ####
# (note: to run need to remove # from list in makenetwork function above)
# makenetwork(re_data_sp_site[1,], groups=groups,plot=T)
stats = apply(re_data_sp_site, 1, function(r,g) makenetwork(r, groups=g), g= groups)
stats = do.call(rbind, stats)

NTW = unique(data[,c('SITE_ID', 'PTL', 'NTL','LAKE_ORIGIN')])
aux = sapply(NTW[,'SITE_ID'], function(r) which(r==data_sp_site[,'SITE_ID']))
NTW = cbind(NTW, mod = stats[aux, ])
NTW$LAKE_ORIGIN <- as.factor(NTW$LAKE_ORIGIN)

# Modularity vs Phosphorus and Nitrogen
# The modularity can be either positive or negative, with positive values indicating the
# possible presence of community structure. Negative values indicate disassortative mixing
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.mod))) + geom_point() + labs(x="PTL",y="log10(Modularity)") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.mod))) + geom_point() + labs(x="NTL",y="log10(Modularity)") + scale_x_log10()

# Edge number vs Phosphorus and Nitrogen
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.no.edge))) + geom_point() + labs(x="log10(PTL)",y="# edges") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.no.edge))) + geom_point() + labs(x="log10(NTL)",y="# edges") + scale_x_log10()

# Cluster number vs Phosphorus and Nitrogen
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.no.clust))) + geom_point() + labs(x="log10(PTL)",y="# edges") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.no.clust))) + geom_point() + labs(x="log10(NTL)",y="# edges") + scale_x_log10()

# Max cluster size vs Phosphorus and Nitrogen
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.max.clust.size))) + geom_point() + labs(x="log10(PTL)",y="Max # clusters") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.max.clust.size))) + geom_point() + labs(x="log10(NTL)",y="Max # clusters") + scale_x_log10()

# Avg cluster size vs Phosphorus and Nitrogen
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.avg.clust.size))) + geom_point() + labs(x="log10(PTL)",y="Avg # clusters") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.avg.clust.size))) + geom_point() + labs(x="log10(NTL)",y="Avg # clusters") + scale_x_log10()

# Avg closeness vs Phosphorus and Nitrogen
ggplot(data=NTW, aes(x=PTL,y=as.numeric(mod.avg.closeness))) + geom_point() + labs(x="log10(PTL)",y="Avg # closeness") + scale_x_log10()
ggplot(data=NTW, aes(x=NTL,y=as.numeric(mod.avg.closeness))) + geom_point() + labs(x="log10(NTL)",y="Avg # closeness") + scale_x_log10()

# Does LCBD of each site relate to network prop? ####
# data_sp_site <- data_sp_site[aux,]
# data_sp_site <- data_sp_site[,-c(1)] # remove SITE_ID column
# mat.chord <- decostand(data_sp_site,"normalize")
# s <- scale(mat.chord, center=TRUE, scale=FALSE)^2
# SStotal <- sum(s)
# BD <- SStotal/(ncol(s)-1)
# LCBD <- apply(s, 1, sum)/SStotal
# SCBD <- aply(s,2,sum)/SStotal
#
# x <- seq(1:length(LCBD))
# y <- rep(1,length(LCBD))
# plot(x,y,cex=LCBD*50,pch=19,col="firebrick2")
#
# x <- seq(1:length(SCBD))
# y <- rep(1,length(SCBD))
# plot(x,y,cex=SCBD*20,pch=19,col="firebrick2")
#
# sort(SCBD*20,decreasing=T)[1:20]
# sort(LCBD*20,decreasing=T)[1:20]
#
# NTW = cbind(NTW,LCBD)
# # Avg cluster size vs LCBD
# plot(NTW[,'LCBD'], NTW[,'mod.avg.clust.size'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Modularity vs LCBD. Negative values indicate disassortative mixing
# plot(NTW[,'LCBD'], NTW[,'mod.mod'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Edge number vs LCBD
# plot(NTW[,'LCBD'], NTW[,'mod.no.edge'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Cluster number vs LCBD
# plot(NTW[,'LCBD'], NTW[,'mod.no.clust'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Max cluster size vs LCBD
# plot(NTW[,'LCBD'], NTW[,'mod.max.clust.size'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Avg cluster size vs Phosphorus, Nitrogen and TN:TP ratio
# plot(NTW[,'LCBD'], NTW[,'mod.avg.clust.size'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
#
# # Connectance vs Phosphorus, Nitrogen and TN:TP ratio
# plot(NTW[,'LCBD'], NTW[,'mod.connectance'], col=NTW[,'LAKE_ORIGIN'], pch=as.numeric(NTW[,'LAKE_ORIGIN']))
