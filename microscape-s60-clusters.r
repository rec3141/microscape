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

options(menu.graphics=FALSE)

# setup functions and locations
source("microscape-s00-setup.r")

# read from disk
# these are the metadata-filtered versions
proptab <- readRDS(file.path(outfolder,"proptab_filt.rds")) # 

# Clustering
# first do PCA / the difference between PCA on x and t(x) is dependent only on which side centering is done

# Cluster samples using Bray-Curtis distances
sample.bray.dist <- parDist(proptab,method="bray",threads=14)
sample.bray.dist[is.na(sample.bray.dist)] <- 1
# PCA
sample.bray.pca <- fast.prcomp(sample.bray.dist)
sample.bray.pca.uniq <- unique(sample.bray.pca$x)
# bh-tsne
sample.tsne.bray.list <- Rtsne(sample.bray.pca.uniq,verbose=TRUE, theta=0.5, pca=FALSE)
sample.tsne.bray <- sample.tsne.bray.list$Y
# add duplicates back in
sample.tsne.bray.full <- addDups(sample.bray.dist, sample.tsne.bray)
# traveling salesperson algorithm
sample.tsne.bray.dist <- dist(sample.tsne.bray.full)
sample.tsp.bray <- solve_TSP(TSP(sample.tsne.bray.dist))
sample.tsp.bray.ind <- c(sample.tsp.bray, sample.tsp.bray[1])
plot(sample.tsne.bray.full, pch='.', main="Sample Bray")
lines(sample.tsne.bray.full[sample.tsp.bray.ind,])

# Cluster sequences using Bray-Curtis distances
seq.bray.dist <- parDist(t(proptab),method="bray",threads=14) # REDO THIS PATH
seq.bray.dist[is.na(seq.bray.dist)] <- 1
seq.bray.pca <- fast.prcomp(seq.bray.dist)
seq.bray.pca.uniq <- unique(seq.bray.pca$x)
seq.tsne.bray.list <- Rtsne(seq.bray.pca.uniq,verbose=TRUE, theta=0.5, pca=FALSE)
seq.tsne.bray <- seq.tsne.bray.list$Y
seq.tsne.bray.full <- addDups(seq.bray.dist, seq.tsne.bray)
# reorder based on clustering using traveling salesperson algorithm
seq.tsne.bray.dist <- dist(seq.tsne.bray.full)
seq.tsp.bray <- solve_TSP(TSP(seq.tsne.bray.dist))
seq.tsp.bray.ind <- c(seq.tsp.bray, seq.tsp.bray[1])
plot(seq.tsne.bray.full, pch='.', main="Sequence Bray")
lines(seq.tsne.bray.full[seq.tsp.bray.ind,])

# 
# # Cluster samples using Euclidean distances
# sample.euclidean.dist <- parDist(proptab,method="euclidean",threads=14)
# sample.euclidean.dist[is.na(sample.euclidean.dist)] <- 1
# sample.euclidean.pca <- fast.prcomp(sample.euclidean.dist)
# sample.euclidean.pca.uniq <- unique(sample.euclidean.pca$x)
# sample.tsne.euclidean.list <- Rtsne(sample.euclidean.pca.uniq,verbose=TRUE, theta=0.5, pca=FALSE)
# sample.tsne.euclidean <- sample.tsne.euclidean.list$Y
# sample.tsne.euclidean.full <- addDups(sample.euclidean.dist, sample.tsne.euclidean)
# # reorder based on clustering using traveling salesperson algorithm
# sample.tsne.euclidean.dist <- dist(sample.tsne.euclidean.full)
# sample.tsp.euclidean <- solve_TSP(TSP(sample.tsne.euclidean.dist))
# sample.tsp.euclidean.ind <- c(sample.tsp.euclidean, sample.tsp.euclidean[1])
# plot(sample.tsne.euclidean.full, pch='.', main="Sample Euclidean")
# lines(sample.tsne.euclidean.full[sample.tsp.euclidean.ind,])
# 
# # Cluster sequences using Euclidean distances
# seq.euclidean.dist <- parDist(t(proptab),method="euclidean",threads=14)
# seq.euclidean.dist[is.na(seq.euclidean.dist)] <- 1
# seq.euclidean.pca <- fast.prcomp(seq.euclidean.dist)
# seq.euclidean.pca.uniq <- unique(seq.euclidean.pca$x)
# seq.tsne.euclidean.list <- Rtsne(seq.euclidean.pca.uniq,verbose=TRUE, theta=0.5, pca=FALSE)
# seq.tsne.euclidean <- seq.tsne.euclidean.list$Y
# seq.tsne.euclidean.full <- addDups(seq.euclidean.dist, seq.tsne.euclidean)
# # reorder based on clustering using traveling salesperson algorithm
# seq.tsne.euclidean.dist <- dist(seq.tsne.euclidean.full)
# seq.tsp.euclidean <- solve_TSP(TSP(seq.tsne.euclidean.dist))
# seq.tsp.euclidean.ind <- c(seq.tsp.euclidean, seq.tsp.euclidean[1])
# plot(seq.tsne.euclidean.full, pch='.', main="Sequence Euclidean")
# lines(seq.tsne.euclidean.full[seq.tsp.euclidean.ind,])

dev.off()

saveRDS(sample.bray.dist,file.path(outfolder,"sample_bray_dist.rds"))
saveRDS(seq.bray.dist,file.path(outfolder,"seq_bray_dist.rds"))
# saveRDS(sample.euclidean.dist,file.path(outfolder,"sample_euclidean_dist.rds"))
# saveRDS(seq.euclidean.dist,file.path(outfolder,"seq_euclidean_dist.rds"))

saveRDS(sample.bray.pca,file.path(outfolder,"sample_bray_pca.rds"))
saveRDS(seq.bray.pca,file.path(outfolder,"seq_bray_pca.rds"))
# saveRDS(sample.euclidean.pca,file.path(outfolder,"sample_euclidean_pca.rds"))
# saveRDS(seq.euclidean.pca,file.path(outfolder,"seq_euclidean_pca.rds"))

saveRDS(sample.tsne.bray.full,file.path(outfolder,"sample_bray_tsne.rds"))
saveRDS(seq.tsne.bray.full,file.path(outfolder,"seq_bray_tsne.rds"))
# saveRDS(sample.tsne.euclidean.full,file.path(outfolder,"sample_euclidean_tsne.rds"))
# saveRDS(seq.tsne.euclidean.full,file.path(outfolder,"seq_euclidean_tsne.rds"))

saveRDS(sample.tsp.bray,file.path(outfolder,"sample_bray_tsp.rds"))
saveRDS(seq.tsp.bray,file.path(outfolder,"seq_bray_tsp.rds"))
# saveRDS(sample.tsp.euclidean,file.path(outfolder,"sample_euclidean_tsp.rds"))
# saveRDS(seq.tsp.euclidean,file.path(outfolder,"seq_euclidean_tsp.rds"))

