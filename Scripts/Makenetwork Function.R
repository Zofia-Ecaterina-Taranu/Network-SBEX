#**********************************************************************************************#
#                             makenetwork functions                                         ####
#**********************************************************************************************#

# http://web.stanford.edu/class/stats366/tables.html

# Functions that inputs an abundance table as a matrix and makes a graph adjacency
# matrix ready to plot with the igraph package.

makenetwork_pa = function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="jaccard"){

  # abundance is the abundance table with no rows that are all zero
  # plotgraph is a toggle for whether a plotted network is requested
  # if commun=TRUE the function will also output the community groups
  # threshold is the number that is fixed for when prsence occurs
  # for instance, threshold>1 means that species with 2 or more reads

  # are considered present
  require(vegan); require(igraph)

  # Keep the original row numbers as labels
  # this is for the later analysis of groups
  if (is.na(dimnames(abund)[[1]][1])) dimnames(abund)=list(1:nrow(abund),1:ncol(abund))
  # Only take the rows where there are at least one value over threshold
  abundance = abund[rowSums(abund)>threshold,]
  n=nrow(abundance)
  # Convert to 1,0 binary matrix for input to vegdist. -0 converts to numeric
  presenceAbsence = (abundance > threshold) - 0
  # Compute the Jaccard distance between the rows, this will only make points
  # closer if they actually present together
  # You could use any of the other distances in vegan or elsewhere
  jaccpa = vegdist(presenceAbsence, method)
  # Distances in R are vectors by default, we make them into matrices
  jaccm = as.matrix(jaccpa)
  coinc = matrix(0,n,n)
  ind1 = which((jaccm>0 & jaccm<(1-incommon)), arr.ind=TRUE)
  coinc[ind1] = 1
  dimnames(coinc)=list(dimnames(abundance)[[1]],dimnames(abundance)[[1]])
  # If using the network package create the graph with
  # g<-as.network.matrix(coinc,matrix.type="adjacency")
  # Here I use the igraph adjacency command
  ig = graph.adjacency(coinc)
  # What class is this?
  class(ig)
  # Take out the isolates
  isolates = V(ig)[ degree(ig)==0 ]
  ignoisol = delete.vertices(ig, V(ig)[ degree(ig) == 0])
  if (plotgraph==TRUE){
    plot(ignoisol, layout=layout.fruchterman.reingold,
         vertex.size=0.6, vertex.label.dist=0.1,
         edge.arrow.mode="-",vertex.color="red",
         vertex.label=NA,edge.color="blue")
    title("p-a dist")
  }

  if (community==TRUE){
    communitywalk = walktrap.community(ignoisol)
    nonisolates= V(ig)[ degree(ig)!=0 ]
    group0 = nonisolates[which(communitywalk$membership==0)]
    group1 = nonisolates[which(communitywalk$membership==1)]
    # You can then play around with coloring and labelling in the graph
    # For help don't forget to look up plot.igraph or plot.network
    # not just plot as it inherits the plot method appropriate to its class
    groups=list(group0,group1,ig)
    return(groups)
  }
}

makenetwork_abund = function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="bray"){

  # abundance is the abundance table with no rows that are all zero
  # plotgraph is a toggle for whether a plotted network is requested
  # if commun=TRUE the function will also output the community groups
  # threshold is the number that is fixed for when presence occurs

  # are considered present
  require(vegan); require(igraph)

  # Keep the original row numbers as labels
  # this is for the later analysis of groups
  if (is.na(dimnames(abund)[[1]][1])) dimnames(abund)=list(1:nrow(abund),1:ncol(abund))
  # Only take the rows where there are at least one value over threshold
  abundance = abund[rowSums(abund)>threshold,]
  n=nrow(abundance)

  # Compute the Bray-Curtis distance between the rows
  # You could use any of the other distances in vegan or elsewhere
  bcabund = vegdist(abundance, method)
  # Distances in R are vectors by default, we make them into matrices
  bcm = as.matrix(bcabund)
  coinc = matrix(0,n,n)
  ind1 = which((bcm>0 & bcm<(1-incommon)), arr.ind=TRUE)
  coinc[ind1] = 1
  dimnames(coinc)=list(dimnames(abundance)[[1]],dimnames(abundance)[[1]])
  # If using the network package create the graph with
  # Here I use the igraph adjacency command
  ig = graph.adjacency(coinc)
  # What class is this?
  class(ig)
  # Take out the isolates
  isolates = V(ig)[ degree(ig)==0 ]
  ignoisol = delete.vertices(ig, V(ig)[ degree(ig) == 0])
  if (plotgraph==TRUE){
    plot(ignoisol, layout=circle,
         vertex.size=0.6, vertex.label.dist=0.1,
         edge.arrow.mode="-",vertex.color="red",
         vertex.label=NA,edge.color="blue")
    title("Abund dist")
  }

  if (community==TRUE){
    communitywalk = walktrap.community(ignoisol)
    nonisolates= V(ig)[ degree(ig)!=0 ]
    group0 = nonisolates[which(communitywalk$membership==0)]
    group1 = nonisolates[which(communitywalk$membership==1)]
    # You can then play around with coloring and labelling in the graph
    # For help don't forget to look up plot.igraph or plot.network
    # not just plot as it inherits the plot method appropriate to its class
    groups=list(group0,group1,ig)
    return(groups)
  }
}

#************************************************************************************************#
#                                    NLA network analysis                                   ####
#************************************************************************************************#

# Loading data ####
setwd("../Clean data")
# data <- read.csv("../Clean data/finally_hooray.csv")
data <- read.csv("full_rough.csv")
# levels(data$TAXATYPE)

gr = c(10,30,60)
# gr = 1:12*5 # grouping by interval of 5 up to 60, then just 60 onwards
data$group = sapply(data$PTL,function(r) if(r>max(gr)) max(gr)+1 else which(gr>=r)[1])

# Subsetting by PTL ###
# Note: different groups have different nodes
# data_oligo <- data[data$PTL >= 0 & data$PTL <=10,]
data_oligo <- data[data$group == 1,]
# length(unique(data_oligo$SITE_ID))
# data_meso <- data[data$PTL > 10 & data$PTL <= 30,]
data_meso <- data[data$group == 2,]
# length(unique(data_meso$SITE_ID))
# data_eu <- data[data$PTL > 30 & data$PTL <= 60,]
data_eu <- data[data$group == 3,]
# length(unique(data_eu$SITE_ID))
# data_hyper <- data[data$PTL > 60,]
data_hyper <- data[data$group == max(gr)+1,]
# length(unique(data_hyper$SITE_ID))

# length(unique(data$SITE_ID)) # length of full dataset = 1157

# Cast matrices ####
library(reshape2)

# All sites
data_sp_site <- dcast(data,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
data_sp_taxa <- dcast(data,TAXANAME ~ SITE_ID, value.var="abund_ml",fun=mean,fill=0)

# head(data_oligo_sp_site)
# data_sp_site <- round(data_sp_site[,2:ncol(data_sp_site)])
data_sp_taxa <- round(data_sp_taxa[,2:ncol(data_sp_taxa)])

# Oligotrophic sites
# data_oligo_sp_site <- dcast(data_oligo,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
data_oligo_sp_taxa <- dcast(data_oligo,TAXANAME ~ SITE_ID, value.var="abund_ml",fun=mean,fill=0)
# head(data_oligo_sp_site)
# data_oligo_sp_site <- round(data_oligo_sp_site[,2:ncol(data_oligo_sp_site)])
data_oligo_sp_taxa <- round(data_oligo_sp_taxa[,2:ncol(data_oligo_sp_taxa)])

# Mesotrophic sites
# data_meso_sp_site <- dcast(data_meso,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
data_meso_sp_taxa <- dcast(data_meso,TAXANAME ~ SITE_ID, value.var="abund_ml",fun=mean,fill=0)
# head(data_meso_sp_site)
# data_meso_sp_site <- round(data_meso_sp_site[,2:ncol(data_meso_sp_site)])
data_meso_sp_taxa <- round(data_meso_sp_taxa[,2:ncol(data_meso_sp_taxa)])

# Eutrophic sites
# data_eu_sp_site <- dcast(data_eu,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
data_eu_sp_taxa <- dcast(data_eu,TAXANAME ~ SITE_ID, value.var="abund_ml",fun=mean,fill=0)
# head(data_eu_sp_site)
# data_eu_sp_site <- round(data_eu_sp_site[,2:ncol(data_eu_sp_site)])
data_eu_sp_taxa <- round(data_eu_sp_taxa[,2:ncol(data_eu_sp_taxa)])

# Hypereutrophic sites
# data_hyper_sp_site <- dcast(data_hyper,SITE_ID ~ TAXANAME, value.var="abund_ml",fun=mean,fill=0)
data_hyper_sp_taxa <- dcast(data_hyper,TAXANAME ~ SITE_ID, value.var="abund_ml",fun=mean,fill=0)
# head(data_hyper_sp_site)
# data_hyper_sp_site <- round(data_hyper_sp_site[,2:ncol(data_hyper_sp_site)])
data_hyper_sp_taxa <- round(data_hyper_sp_taxa[,2:ncol(data_hyper_sp_taxa)])


# Display networks - Fruchterman Reingold layout ####

# Network on all data colour coded by trophic groups
# ntw <- makenetwork_pa(data_sp_site)
# colbar <- c("chartreuse", "goldenrod","firebrick","dodgerblue")
# plot.igraph(ntw_oligo[[3]],vertex.color=colbar) # Index [[3]] calls the actual IGRAPH

# Networks split by group: comparison of p/a and abundance distance matrices
par(mfrow=c(1,2))

ntw_oligo <- makenetwork_pa(data_oligo_sp_taxa)
ntw_oligo_abund <- makenetwork_abund(data_oligo_sp_taxa)
# plot.igraph(ntw_oligo_abund[[3]],vertex.color=colbar) # Index [[3]] calls the actual IGRAPH

ntw_meso <- makenetwork_pa(data_meso_sp_taxa)
ntw_meso_abund <- makenetwork_abund(data_meso_sp_taxa)

ntw_eu <- makenetwork_pa(data_eu_sp_taxa)
ntw_eu_abund <- makenetwork_abund(data_eu_sp_taxa)

ntw_hyper <- makenetwork_pa(data_hyper_sp_taxa)
ntw_hyper_abund <- makenetwork_abund(data_hyper_sp_taxa)

# Network properties ####

# Modularity pa
modularity(walktrap.community(ntw_oligo[[3]]))
modularity(walktrap.community(ntw_meso[[3]]))
modularity(walktrap.community(ntw_eu[[3]]))
modularity(walktrap.community(ntw_hyper[[3]]))

# Modularity abundance
modularity(walktrap.community(ntw_oligo_abund[[3]]))
modularity(walktrap.community(ntw_meso_abund[[3]]))
modularity(walktrap.community(ntw_eu_abund[[3]]))
modularity(walktrap.community(ntw_hyper_abund[[3]]))
# Modularity indicates the number of separate groups and is expected to increase the stability
# of networks since disturbances are less likely to spread across different modules.

# Here, modularity increased from oligo- to hypereutropic. This means that the network has
# more and more distinct groups and is more stable (??) as the lakes become more eutrophic.
# Interpretation: Communities in hypereutrophic lakes are less similar to each other and
#become more disconnected from the rest of the network. There are more singular nodes.
#How can we get at what lakes do stay similar?
#There are more singular components (single nodes)  or components with just a few
#nodes in hypereutrophic lakes.

# Path lengths, regular graph
average.path.length(ntw_oligo[[3]])
average.path.length(ntw_meso[[3]])
average.path.length(ntw_eu[[3]])
average.path.length(ntw_hyper[[3]])

diameter(ntw_oligo[[3]]) # longest path length
diameter(ntw_meso[[3]])
diameter(ntw_eu[[3]])
diameter(ntw_hyper[[3]])
# Small world property refers to the situation wherein (a) the shortest-path distance between
# pairs of vertices is generally quite small, but (b) the clustering is relatively high.
# In our network, the average path length is less than five and even the longest of paths is
# not much bigger

# number of vertices
vcount(ntw_oligo[[3]])
vcount(ntw_meso[[3]])
vcount(ntw_eu[[3]])
vcount(ntw_hyper[[3]])

# number of edges
ecount(ntw_oligo[[3]])
ecount(ntw_meso[[3]])
ecount(ntw_eu[[3]])
ecount(ntw_hyper[[3]])

# Density
e <- ecount(ntw_oligo[[3]])
etot <- e*(e-1)/2 # total edge
e/etot # realized/ total

e <- ecount(ntw_meso[[3]])
etot <- e*(e-1)/2 # total edge
e/etot # realized/ total

e <- ecount(ntw_eu[[3]])
etot <- e*(e-1)/2 # total edge
e/etot # realized/ total

e <- ecount(ntw_hyper[[3]])
etot <- e*(e-1)/2 # total edge
e/etot # realized/ total

# number of components
no.clusters(ntw_oligo[[3]])
no.clusters(ntw_meso[[3]])
no.clusters(ntw_eu[[3]])
no.clusters(ntw_hyper[[3]])

# How big are these?
table(clusters(ntw_oligo[[3]])$csize)
table(clusters(ntw_meso[[3]])$csize)
table(clusters(ntw_eu[[3]])$csize)
table(clusters(ntw_hyper[[3]])$csize)



# igraph has a function for generating random networks of varying size and connectance.
graph.random.gnp<-erdos.renyi.game(n=15,p.or.m=.5,type="gnp",directed=T)
plot.igraph(graph.random.gnp)
# Here we have created a random directed graph with 15 species ("n") and a connectance ("p") of 0.5

# Similarly we can set the number of links that we want in the system to a value "m" that we specify.
graph.random.gnm<-erdos.renyi.game(n=15,p.or.m=112,type="gnm",directed=T)
plot.igraph(graph.random.gnm)

# In igraph we can use the barabasi.game() function (Power law distribution)
graph.barabasi.1<-barabasi.game(n=15,power=0.75)
plot.igraph(graph.barabasi.1)

# Correlation between centrality measures
cent <- list(`Degree`=degree(ntw[[3]]),
             `Closeness`=closeness(ntw[[3]]),
             `Betweenness`=betweenness(ntw[[3]]),
             `Eigenvector`=evcent(ntw[[3]])$vector,
             `PageRank`=page.rank(ntw[[3]])$vector)

# Pairs plot
pairs(cent, lower.panel=function(x,y) {
  usr <- par("usr")
  text(mean(usr[1:2]), mean(usr[3:4]), round(cor(x,y), 3), cex=2, col="blue")
} )

# Degree distribution
hist(degree.distribution(ntw[[3]])) # Power law distribution tends to occur with networks instead
# of a bell curve. A lot of people with very low degree, and long tail with with very few with large number of degrees.

is.directed(ntw[[3]]) # for directed the mode is ignored
is.loop(ntw[[3]])
# Degree in-distribution
# plot(degree.distribution(ntw[[3]], mode="in"), log="xy")
# Degree out-distribution
# plot(degree.distribution(ntw[[3]], mode="out"), log="xy")

# Transitivity of the ring graph
transitivity(ntw[[3]])

# Transitivity of a random graph of the same size
g <- erdos.renyi.game(vcount(ntw[[3]]), ecount(ntw[[3]]), type="gnm")
transitivity(g)

# Transitivity of a random graph with the same degree distribution
# g <- degree.sequence.game(degree(ntw[[3]], mode="out"), degree(ntw[[3]], mode="in"), method="simple")
# transitivity(g)

# edge.betweenness.community
edge.betweenness.community(ntw[[3]])

# leading.eigenvector.community
leading.eigenvector.community(ntw[[3]])

# Calculate communities
res <- clique.community(ntw[[3]], k=2)
# Paint them to different colors
colbar <- rainbow( length(res)+1 )
for (i in seq(along=res)) {
  V(ntw[[3]])[ res[[i]] ]$color <- colbar[i+1]
}
# Paint the vertices in multiple communities to red
V(ntw[[3]])[ unlist(res)[ duplicated(unlist(res)) ] ]$color <- "red"
plot(ntw[[3]], layout=lay, vertex.label=V(ntw[[3]])$name)

# Function to test regular graph with given size
try.ring.pl <- function(n) {
  g <- watts.strogatz.game(1, n, 3, p=0)
  average.path.length(g)
}
try.ring.pl(e)

# Test a number of regular graphs
ring.size <- seq(100, 1000, by=100)
ring.pl <- sapply(ring.size, try.ring.pl)
plot(ring.size, ring.pl, type="b")

# Path lengths, random graph
rg <- erdos.renyi.game(50, 50*3, type="gnm")
rg$layout <- layout.circle
V(rg)$size <- 3
plot(rg, vertex.label=NA, main="Random graph")
average.path.length(rg)

# Path length of random graphs
try.random.pl <- function(n) {
  g <- erdos.renyi.game(n, n*3, type="gnm")
  average.path.length(g)
}
try.random.pl(100)

# Plot network size vs. average path length
random.pl <- sapply(ring.size, try.random.pl)
plot(ring.size, random.pl, type="b")
plot(ring.size, random.pl, type="b", log="x")

# Transitivity, random graph, by definition
ecount(rg) / (vcount(rg)*(vcount(rg)-1)/2)
transitivity(rg, type="localaverage")

# Rewiring
ig2 <- watts.strogatz.game(1, 50, 3, p=0.1)
ig2$layout <- layout.circle
V(ig2)$size <- 3
plot(ig2, vertex.label=NA)

# Naming vertices
ig <- graph.ring(10)
V(ig)$name <- letters[1:10]

# Sequence of all Edges
E(ntw[[3]])
# All adjacent edges of a vertex
E(ntw[[3]])[ adj(3) ]
# Outgoing edges
E(ig)[ from(3) ]
# Incoming edges
E(ig)[ to(3) ]
# Edges along a path
E(ig, path=c(1,4,5))

# Optimalize modularity
optcom <- optimal.community(ntw[[3]])
V(ntw[[3]])$comm <- membership(optcom)
plot(optcom, ntw[[3]])

# Fit a HRG model to the network
hrg <- hrg.fit(ntw[[3]])
ihrg <- as.igraph(hrg)
ihrg$layout <- layout.reingold.tilford
plot(ihrg, vertex.size=10, edge.arrow.size=0.2)

# Modify
vn <- sub("Actor ", "", V(ihrg)$name)
colbar <- rainbow(length(optcom))
vc <- ifelse(is.na(V(ihrg)$prob), colbar[V(ig)$comm], "darkblue")
V(ihrg)$label <- ifelse(is.na(V(ihrg)$prob), vn, round(V(ihrg)$prob, 2))

par(mar=c(0,0,3,0))
plot(ihrg, vertex.size=10, edge.arrow.size=0.2,
     vertex.shape="none", vertex.label.color=vc,
     main="Hierarchical network model of the Karate Club")
dendPlot(hrg)

# More modification
ntw[[3]]$layout <- layout.circle
V(ntw[[3]])$color <- "white"
V(ig)[name=="A"]$color <- "orange"
V(ntw[[3]])$size <- 10
V(ntw[[3]])$label.cex <- 1
V(ntw[[3]])$label <- V(ntw[[3]])$name
E(ntw[[3]])$color <- "black"
E(ntw[[3]])$width <- 3
# Plot 'ig' and A's transitivity
tr <- transitivity(ntw[[3]], type="local")[1]
plot(ntw[[3]], main=paste("Transitivity of 'A':", tr))

# Make a very hierarchical graph
g1 <- graph.full(5)
g2 <- graph.ring(5)
g <- g1 + g2
g <- g + edge(1, vcount(g1)+1)
plot(g)
ghrg <- hrg.fit(g)
dendPlot(ghrg)

# Using tkplot()
E(ntw[[3]])$color <- "grey"
V(ntw[[3]])[ degree(ntw[[3]])>=5 ]$color <- "yellow"
tkplot(ntw[[3]])
scc <- clusters(ntw[[3]],"strong")
# Find strongly connected vertices:
v_idx <- which(scc$csize>1)
V(ntw[[3]])[scc$membership==v_idx[2]]
unique(as.vector(ntw[[3]]))[scc$memberships==v_idx[2]]
V(ntw[[3]])[degree(ntw[[3]])>5]


#************************************************************************************************#
#                             load the igraph and rgexf packages                                 #
#************************************************************************************************#

library(rgexf)

# construct the nodes and edges data for gexf conversion
nodes <- data.frame(cbind(V(ntw[[3]]), as.character(V(ntw[[3]]))))
edges <- t(Vectorize(get.edge, vectorize.args='id')(ntw[[3]], 1:ecount(ntw[[3]])))

# do the conversion
write.gexf(nodes, edges)

#************************************************************************************************#

# Or use the following function
saveAsGEXF = function(g, filepath="converted_graph.gexf")
{

  require(igraph)

  require(rgexf)

  # gexf nodes require two column data frame (id, label)

  # check if the input vertices has label already present

  # if not, just have the ids themselves as the label

  if(is.null(V(g)$label))

    V(g)$label <- as.character(V(g)$name)

  # similarily if edges does not have weight, add default 1 weight

  if(is.null(E(g)$weight))

    E(g)$weight <- rep.int(1, ecount(g))

  nodes <- data.frame(cbind(V(g), V(g)$label))

  edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))

  # combine all node attributes into a matrix (and take care of & for xml)

  vAttrNames <- setdiff(list.vertex.attributes(g), "label")

  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))),

                         stringsAsFactors = FALSE)

  # combine all edge attributes into a matrix (and take care of & for xml)

  eAttrNames <- setdiff(list.edge.attributes(g), "weight")

  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))),

                         stringsAsFactors = FALSE)

  # generate the gexf object

  output <- write.gexf(nodes, edges,

                       edgesWeight=E(g)$weight,

                       edgesAtt = edgesAtt,

                       nodesAtt = nodesAtt)

  print(output, filepath, replace=T)

}
saveAsGEXF(ntw[[3]])

#************************************************************************************************#

# Or: http://www.vesnam.com/Rblog/viznets2/

# Plotting networks in R
# An example how to use R and rgexf package to create a .gexf file for network visualization in Gephi

library("plyr")

# Read a data set.
# Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
dataSet <- read.table("lesmis.txt", header = FALSE, sep = "\t")

# Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
gD <- simplify(graph.data.frame(dataSet, directed=FALSE))

# Print number of nodes and edges
vcount(gD)
ecount(gD)

# Try with ntw[[3]]
gD <- ntw[[3]]


# Calculate some node properties and node similarities that will be used to illustrate
# different plotting abilities

# Calculate degree for all nodes
degAll <- degree(gD, v = V(gD), mode = "all")

# Calculate betweenness for all nodes
betAll <- betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)

# Calculate Dice similarities between all pairs of nodes
dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")

# Add new node/edge attributes based on the calculated node properties/similarities

gD <- set.vertex.attribute(gD, "degree", index = V(gD), value = degAll)
gD <- set.vertex.attribute(gD, "betweenness", index = V(gD), value = betAll.norm)

# Check the attributes
# summary(gD)

F1 <- function(x) {data.frame(V4 = dsAll[which(V(gD)$name == as.character(x$V1)), which(V(gD)$name == as.character(x$V2))])}
dataSet.ext <- ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))

gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0)
gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)

# The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
# and for that reason these values cannot be assigned directly

E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)

# Check the attributes
summary(gD)

# Create a dataframe nodes: 1st column - node ID, 2nd column -node name
nodes_df <- data.frame(ID = c(1:vcount(gD)), NAME = V(gD)$name)
# Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
edges_df <- as.data.frame(get.edges(gD, c(1:ecount(gD))))

# Define node and edge attributes - these attributes won't be directly used for network visualization, but they
# may be useful for other network manipulations in Gephi
#
# Create a dataframe with node attributes: 1st column - attribute 1 (degree), 2nd column - attribute 2 (betweenness)
nodes_att <- data.frame(DEG = V(gD)$degree, BET = V(gD)$betweenness)
#
# Create a dataframe with edge attributes: 1st column - attribute 1 (weight), 2nd column - attribute 2 (similarity)
edges_att <- data.frame(WGH = E(gD)$weight, SIM = E(gD)$similarity)

# Define node/edge visual attributes - these attributes are the ones used for network visualization
#
# Calculate node coordinate - needs to be 3D
#nodes_coord <- as.data.frame(layout.fruchterman.reingold(gD, weights = E(gD)$similarity, dim = 3, niter = 10000))
# We'll cheat here, as 2D coordinates result in a better (2D) plot than 3D coordinates
nodes_coord <- as.data.frame(layout.fruchterman.reingold(gD, weights = E(gD)$similarity, dim = 2, niter = 10000))
nodes_coord <- cbind(nodes_coord, rep(0, times = nrow(nodes_coord)))
#
# Calculate node size
# We'll interpolate node size based on the node betweenness centrality, using the "approx" function
approxVals <- approx(c(1, 5), n = length(unique(V(gD)$betweenness)))
# And we will assign a node size for each node based on its betweenness centrality
nodes_size <- sapply(V(gD)$betweenness, function(x) approxVals$y[which(sort(unique(V(gD)$betweenness)) == x)])
#
# Define node color
# We'll interpolate node colors based on the node degree using the "colorRampPalette" function from the "grDevices" library
library("grDevices")
# This function returns a function corresponding to a collor palete of "bias" number of elements
F2 <- colorRampPalette(c("#F5DEB3", "#FF0000"), bias = length(unique(V(gD)$degree)), space = "rgb", interpolate = "linear")
# Now we'll create a color for each degree
colCodes <- F2(length(unique(V(gD)$degree)))
# And we will assign a color for each node based on its degree
nodes_col <- sapply(V(gD)$degree, function(x) colCodes[which(sort(unique(V(gD)$degree)) == x)])
# Transform it into a data frame (we have to transpose it first)
nodes_col_df <- as.data.frame(t(col2rgb(nodes_col, alpha = FALSE)))
# And add alpha (between 0 and 1). The alpha from "col2rgb" function takes values from 0-255, so we cannot use it
nodes_col_df <- cbind(nodes_col_df, alpha = rep(1, times = nrow(nodes_col_df)))
# Assign visual attributes to nodes (colors have to be 4dimensional - RGBA)
nodes_att_viz <- list(color = nodes_col_df, position = nodes_coord, size = nodes_size)

# Assign visual attributes to edges using the same approach as we did for nodes
F2 <- colorRampPalette(c("#FFFF00", "#006400"), bias = length(unique(E(gD)$weight)), space = "rgb", interpolate = "linear")
colCodes <- F2(length(unique(E(gD)$weight)))
edges_col <- sapply(E(gD)$weight, function(x) colCodes[which(sort(unique(E(gD)$weight)) == x)])
edges_col_df <- as.data.frame(t(col2rgb(edges_col, alpha = FALSE)))
edges_col_df <- cbind(edges_col_df, alpha = rep(1, times = nrow(edges_col_df)))
edges_att_viz <-list(color = edges_col_df)

# Write the network into a gexf (Gephi) file
#write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesWeight = E(gD)$weight, edgesAtt = edges_att, nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, defaultedgetype = "undirected", output = "lesmis.gexf")
# And without edge weights
write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesAtt = edges_att, nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, defaultedgetype = "undirected", output = "lesmis.gexf")
write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesAtt = edges_att, nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, defaultedgetype = "undirected", output = "mites.gexf")
