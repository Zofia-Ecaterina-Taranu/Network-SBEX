## A script that aggregates BIOMASS across TAXATYPES
## First a network graph is constructed for each SITE_ID, where nodes are all the species
## and an edges between a pair of species exist if both do exit, with edge weight as the product of BIOMASSES.
## i.e., SITE_ID ~ TAXANAME by BIOMASS.
## Then for each pair of TAXATYPES the scripts aggregates all the edge weights connecting all species in both pairs.
## Aggregation could be done by a SUM or a PROD, or any desired function.
library(reshape2)
library(igraph)
rm(list=ls())
## Loading data
full <- read.csv("../Clean data/full_combined.csv", as.is=TRUE)
data <- full[!is.na(full$BIOMASS),] # Removing rows with NAs for BIOMASS (24056 obs. remaining out of 48755)

# Note (from Zo): Suggestion to replace TAXATYPE by T_GROUP for Zooplankton to break them down a bit more
data$TAXATYPE[data$TAXATYPE=="Zooplankton"] <- data$T_GROUP[data$TAXATYPE=="Zooplankton"]

## Constructing a matrix of SITE_ID's as rows, TAXANAME as columns filled with each species biomass.
data_sp_site <- dcast(data,SITE_ID ~ TAXANAME, value.var="BIOMASS",fun=mean,fill=0) # Now have 1154 unique sites


## Aggregating the TAXATYPE
aggregate.TAXATYPE<-function(spe, taxa.name.type, agg.fun = mean){
    ## spe is the SITE_ID ~ TAXANAME matrix, where the first column.
    ## taxa.name.type is a matrix two columns (taxaname, taxatype) and n rows,
    ## for example, taxa.name.type= unique(data[ ,c('TAXANAME', 'TAXATYPE')]).
    ## agg.fun is the aggregation function, default is the mean.

    if(colnames(spe)[1]== 'SITE_ID') spe = spe[,-1]
    if(!is.matrix(spe)) spe = as.matrix(spe)
    ## For each species (column) in spe, what TAXATYPE is it.
    spe.per.taxatype = tapply(taxa.name.type[,1], taxa.name.type[,2],
        function(r) which(colnames(spe) %in% r) )
    taxatypes= as.factor(unique(taxa.name.type[,2]))
    ## All pairs and TAXATYPES
    grid= expand.grid(levels(taxatypes), levels(taxatypes))
    ## For each site, return the probabilities for each row in taxa.name.type
    grid.P = apply(spe,1,function(r){
                           r = r %*% t(r)/(sum(r)^2)
                           apply(grid, 1, function(x){
                                           col = spe.per.taxatype[[unlist(x[1])]]
                                           row = spe.per.taxatype[[unlist(x[2])]]
                                           agg.fun(r[row,col])
                                       })
                       })
}

makenetwork.TAXATYPE<-function(groups, prob){
    ## Creates and plots an igraph object.
    ## Given a list of groups
    ## Plotting function.
    taxatypes = as.factor(unique(groups[,2]))
    M= matrix(prob,nrow = nlevels(taxatypes), ncol=nlevels(taxatypes))
    colnames(M)<-levels(groups[,1])
    rownames(M)<-levels(groups[,1])
    ig = graph.adjacency(M, mode= "undirected", weighted = TRUE)
    ly = layout.circle(ig)
    w = E(ig)$weight/mean(E(ig)$weight)
    plot(ig, layout=ly,vertex.color=taxatypes, edge.width=w, edge.label=round(100*w,2))
}


##  Running form here
taxa.name.type = unique(data[ ,c('TAXANAME', 'TAXATYPE')])
prob = aggregate.TAXATYPE(data_sp_site, taxa.name.type)

## Plotting for 4 levels of  PTL. From each level select 4 random sites.
gp =list(g1 = sample(unique(data$SITE_ID[which(data$PTL<=10)]), 4),
         g2 = sample(unique(data$SITE_ID[which(data$PTL>10 & data$PTL <= 30)]), 4),
         g3 = sample(unique(data$SITE_ID[which(data$PTL>30 & data$PTL <= 60)]), 4),
         g4 = sample(unique(data$SITE_ID[which(data$PTL>60)]), 4))


## Plotting loop
k=1
gp_id = sapply(gp, function(r){
                   ## r is the sublist g1, g2, etc.
                   pdf(paste(k, '.pdf'))
                   par(mfrow=c(2,2))
                   aux = which(data_sp_site$SITE_ID %in% r)
                   for (i in 1:length(aux)){
                       makenetwork.TAXATYPE(taxa.name.type,prob[,aux[i]])
                       tl = data$PTL[which(data$SITE_ID==r[i])[1]]
                       title(paste('PTL:', tl))
                   }
                   dev.off()
                   k<<-k+1
               })
