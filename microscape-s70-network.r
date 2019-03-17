# NETWORK ANALYSIS
#install.packages(c("SpiecEasi","parallelDist","Matrix","Rtsne","TSP","gmodels","igraph","parallel"))
library(SpiecEasi) # for sparcc, spieceasi
library(parallelDist) #for parDist
library(Matrix) #for tril
library(Rtsne) #for Rtsne
library(TSP) #for TSP
library(gmodels) #for fast.prcomp
library(igraph) #for components
library(parallel) #for detectCores
source("FastPCA.R")
library(reshape2)

source("microscape-s00-setup.r")

options(menu.graphics=FALSE)

#output.path
outfolder <- "out_dada"

# read from disk
normtab <- readRDS(file.path(outfolder,"normtab_filt.rds")) # 

# run parallel sparcc
# sparcc.out <- parsparcc(proptab, iter=detectCores()-1, inner_iter=5, th=0.1)
# sparcc.out <- parsparcc(proptab, iter=7, inner_iter=15, th=0.1) #crashed
# for computational sake, only compute correlations on sequences present in at least 5 samples (9472/19073 = 50%)
# crashed after n=21 with cores=13
# n=9 cores=5 took 10 hours
# n=5 cores=5 ran out of memory at i=5
# bonk should use a consistent number of iterations
# n=5 i=5 cores=1 took 43 hours
# new merged samples:
# crashed after n=19 with cores=13
# approximate time in days to complete: y = 121.3x^-2.209
#     #sparcc.out <- sparcc(proptab, iter=15, inner_iter=15, th=0.1)

# do subset based on TSNE clustering ?


#i is minimum number of samples
#crashed at i=45, iter=20, inner=20, th=0.1, cores=13
for(i in seq(75,5,-10)) {
    normtab_small <- normtab[,colSums(normtab>0)>i]  ### THIS MEANS I NEED TO SAVE THE NAMES WITH THE NETWORK_MELT BELOW
	small.cols <- match(colnames(normtab_small),colnames(normtab))
	small.esvs <- paste0("ESV_",small.cols)
    print(date())
    print(i)
    print(dim(normtab_small))
    sparcc.out <- parsparcc(normtab_small, iter=20, inner_iter=20, th=0.1, cores=13)
    saveRDS(sparcc.out,file=file.path(outfolder,paste0("sparcc_out_",i,".rds")))

	# Cluster sequences using sparcc correlations
	seq.sparcc.dist <- sparcc.out$Cor
	rownames(seq.sparcc.dist) <- small.esvs
	colnames(seq.sparcc.dist) <- small.esvs
	seq.sparcc.dist[is.na(seq.sparcc.dist)] <- 1
#	seq.sparcc.pca <- fast.prcomp(seq.sparcc.dist)
	seq.sparcc.pca <- FastPCA(seq.sparcc.dist,50)
	seq.sparcc.pca.uniq <- unique(seq.sparcc.pca$x)
	seq.tsne.sparcc.list <- Rtsne(seq.sparcc.pca.uniq,verbose=TRUE, theta=0.5, max_iter=2000, pca=FALSE)
	seq.tsne.sparcc <- seq.tsne.sparcc.list$Y
	if(!identical(seq.sparcc.pca$x,seq.sparcc.pca.uniq)) {
		seq.tsne.sparcc.full <- addDups(seq.sparcc.dist, seq.tsne.sparcc)
	} else {
		seq.tsne.sparcc.full <- seq.tsne.sparcc
	}
	# reorder based on clustering using traveling salesperson algorithm
	seq.tsne.sparcc.dist <- dist(seq.tsne.sparcc.full)
	seq.tsp.sparcc <- solve_TSP(TSP(seq.tsne.sparcc.dist))
	seq.tsp.sparcc.ind <- c(seq.tsp.sparcc, seq.tsp.sparcc[1])
	plot(seq.tsne.sparcc.full, pch='.', main="Sequence SparCC")
	lines(seq.tsne.sparcc.full[seq.tsp.sparcc.ind,])

	rownames(seq.tsne.sparcc.full) <- small.esvs
	
	saveRDS(seq.sparcc.dist,file.path(outfolder,"seq_sparcc_dist.rds"))
	saveRDS(seq.sparcc.pca,file.path(outfolder,"seq_sparcc_pca.rds"))
	saveRDS(seq.tsne.sparcc.full,file.path(outfolder,"seq_sparcc_tsne.rds"))
	saveRDS(seq.tsp.sparcc,file.path(outfolder,"seq_sparcc_tsp.rds"))

	# vastly simplified by just saving vectors of weights and colors instead of the whole graph
	# four columns sorted by correlation: node1, node2, correlation, weight

	network_correlations <- sparcc.out$Cor
	network_matrix <- as.matrix(Matrix::tril(network_correlations) + t(Matrix::tril(network_correlations)))/2 #to make it symmetric so igraph doesn't crash

    network_melt <- melt(network_matrix)
    colnames(network_melt) <- c("node1","node2","corr")
    network_melt <- network_melt[network_melt$node1 > network_melt$node2, ]
    network_melt <- network_melt[abs(network_melt$corr)>0.1,]
    network_melt$weight <- network_melt$corr^3/max(abs(network_melt$corr)^3)
    network_melt$color[network_melt$corr<0] <- rgb(abs(network_melt$weight[network_melt$corr<0]),0,0,abs(network_melt$weight[network_melt$corr<0]))
    network_melt$color[network_melt$corr>=0] <- rgb(0,0,abs(network_melt$weight[network_melt$corr>=0]),abs(network_melt$weight[network_melt$corr>=0]))
	network_melt$node1 <- small.esvs[network_melt$node1]
	network_melt$node2 <- small.esvs[network_melt$node2]

	saveRDS(network_melt, file.path(outfolder, "network_melt.rds"))
   
   }

# at i=75
# > dim(network_melt)
# [1] 18724     5

# at old i=19
# >     print(dim(normtab_small))
# [1] 2191 3123

# 	# save a graph at each cutoff
# 	save.graphs <- list()
# 	for(i in 100:1) {
# 		cutoff <- i/100
# 		print(cutoff)
# 
# 		# define graph
# 		network_matrix <- abs(network_correlations) > cutoff
# 		diag(network_matrix) <- 0
# 		if(sum(network_matrix)==0) { save.graphs[[as.character(i)]] <- NULL; print("done"); next }
# 
# 		network_matrix <- as.matrix(Matrix::tril(network_matrix) + t(Matrix::tril(network_matrix))) #to make it symmetric so igraph doesn't crash
# 
# 		network_graph <- adj2igraph(network_matrix)
# 
# 		save.graphs[[paste0("cor_",as.character(i))]] <- network_graph
# 	}
# 	saveRDS(save.graphs, file.path(outfolder,"network_graphs.rds")) #730Mb each... must be a better way but for now...
# 
# 	# save the edge weights at each cutoff
# 	save.weights <- list()
# 	for(i in 100:1) {
# 		print(i)
# 		network_graph <- save.graphs[[paste0("cor_",as.character(i))]]
# 		if(is.null(network_graph)) next
# 		#define weights and colors
# 		# network_matrix inherits from i=1 above
# 		network_weights <- (network_matrix & upper.tri(network_matrix))*network_correlations
# 
# 		#scale the weights
# 		E(network_graph)$weight <- network_weights[abs(network_weights) > cutoff]^3 # cubed preserves sign
# 		#recenter the weights
# 		E(network_graph)$weight <- E(network_graph)$weight/max(abs(E(network_graph)$weight))
# 		#color the rescaled weights
# 		E(network_graph)$color <- rgb(0,0,0,0)
# 		network_neg <- E(network_graph)$weight < 0
# 		network_pos <- E(network_graph)$weight > 0
# 		E(network_graph)$color[network_neg] <- rgb(abs(E(network_graph)$weight[network_neg]),0,1,abs(E(network_graph)$weight[network_neg]))
# 		E(network_graph)$color[network_pos] <- rgb(0,0,abs(E(network_graph)$weight[network_pos]),abs(E(network_graph)$weight[network_pos]))
# 		network_weights[network_weights==1] <- 0
# 
# 		save.weights[[paste0("cor_",as.character(i))]] <- network_weights
# 	}
# 	saveRDS(save.weights, file.path(outfolder,"network_weights.rds"))
# 

# }
