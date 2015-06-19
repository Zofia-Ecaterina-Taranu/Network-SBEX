##  A script that aggregates BIOMASS across TAXATYPES
rm(list=ls())
full <- read.csv("../Clean data/full_rough.csv", as.is=TRUE)
data <- full[!is.na(full$BIOMASS),] # Removing rows with NAs for BIOMASS
# All sites
data_sp_site <- dcast(data,SITE_ID ~ TAXANAME, value.var="BIOMASS",fun=mean,fill=0)

## Aggregating the TAXATYPE
uni = unique(data[ ,c('TAXANAME', 'TAXATYPE')])
g = as.factor(sapply(colnames(data_sp_site[,-1]), function(r) uni[which(r==uni[,1])[1],2]))
grid= expand.grid(levels(g), levels(g), stringsAsFactors=T)
grid.P = apply(data_sp_site[,-1],1,  function(r){
          r = r %*% t(r)/(sum(r)^2)
          #diag(r) <- 0
          apply(grid, 1, function(x){
                   col = which(g==unlist(x[1]))
                   row = which(g==unlist(x[2]))
                   mean(unlist(r[row,col]))
          })
 })

makenetwork.TAXATYPE<-function(groups, prob){
    ## Plotting function
    M= matrix(prob,nrow = nlevels(groups[,1]), ncol=nlevels(groups[,2]))
    colnames(M)<-levels(groups[,1])
    rownames(M)<-levels(groups[,1])

    ig = graph.adjacency(M, mode= "undirected", weighted = TRUE)
    ly = layout.circle(ig)
    w = E(ig)$weight/mean(E(ig)$weight)
    plot(ig, layout=ly,vertex.color=groups[,1], edge.width=w, edge.label=round(100*w,2))
}


## Plotting for 4 level of  PTL. From each level select 4 random sites.
gp =list(g1 = sample(unique(data$SITE_ID[which(data$PTL<=10)]), 4),
         g2 = sample(unique(data$SITE_ID[which(data$PTL>10 & data$PTL <= 30)]), 4),
         g3 = sample(unique(data$SITE_ID[which(data$PTL>30 & data$PTL <= 60)]), 4),
         g4 = sample(unique(data$SITE_ID[which(data$PTL>60)]), 4))


## Plotting loop
k=1
gp_id = sapply(gp, function(r, g){
  ## r is the sublist g1, g2, etc.
  pdf(paste(k, '.pdf'))
  par(mfrow=c(2,2))
  aux = which(data_sp_site$SITE_ID %in% r)
  for (i in 1:length(aux)){
    makenetwork.TAXATYPE(grid,grid.P[,aux[i]])
    tl = data$PTL[which(data$SITE_ID==r[i])[1]]
    title(paste('PTL:', tl))
}
  dev.off()
  k<<-k+1
}, g=groups)
