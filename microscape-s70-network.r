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

#output.path
outfolder <- "out_dada"


# read from disk
normtab <- readRDS(file.path(outfolder,"normtab_final.rds")) # 

# SPARCC correlation networks for sequences
# algorithm requires integer values
# replace sparcc function with parallel sparcc
parsparcc <- function (data, iter = 20, inner_iter = 10, th = 0.1, cores=5)
{
    sparccs <- mclapply(1:iter, function(i) parsparccinner(t(apply(data, 1, norm_diric)), iter = inner_iter, th = th), mc.cores=cores)
    cors <- array(unlist(lapply(sparccs, function(x) x$Cor)),
        c(ncol(data), ncol(data), iter))
    corMed <- apply(cors, 1:2, median)
    covs <- array(unlist(lapply(sparccs, function(x) x$Cov)),
        c(ncol(data), ncol(data), iter))
    covMed <- apply(cors, 1:2, median)
    covMed <- cor2cov(corMed, sqrt(diag(covMed)))
    list(Cov = covMed, Cor = corMed)
}
environment(parsparcc) <- asNamespace('SpiecEasi')

parsparccinner <- function (data.f, T = NULL, iter = 10, th = 0.1)
{
    if (is.null(T))
        T <- av(data.f)
    res.bv <- basis_var(T)
    Vbase <- res.bv$Vbase
    M <- res.bv$M
    cbase <- C_from_V(T, Vbase)
    Cov <- cbase$Cov
    Cor <- cbase$Cor
    excluded <- NULL
#    print(iter)
    for (i in 1:iter) {    
        print(i)
        res.excl <- exclude_pairs(Cor, M, th, excluded)
        M <- res.excl$M
        excluded <- res.excl$excluded
#        print(res.excl)
        if (res.excl$break_flag)
            break
        res.bv <- basis_var(T, M = M, excluded = excluded)
        Vbase <- res.bv$Vbase
        M <- res.bv$M
        K <- M
        diag(K) <- 1
        cbase <- C_from_V(T, Vbase)
        Cov <- cbase$Cov
        Cor <- cbase$Cor
    }
    gc()
    list(Cov = Cov, Cor = Cor, i = i, M = M, excluded = excluded)
}
environment(parsparccinner) <- asNamespace('SpiecEasi')


# function to add duplicates back in
# finds zero distances from distance matrix df.dist, adds de-duped rows from df.xy back into df.xy
addDups <- function(df.dist, df.xy) {
    #stackoverflow solution to insert row: https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended/11562428#11562428
    insertRow2 <- function(existingDF, newrow, r) {
      existingDF <- rbind(existingDF,newrow)
      existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
      row.names(existingDF) <- 1:nrow(existingDF)
      return(existingDF)  
    }

    mat <- as.matrix(df.dist)
    zeros <- which(mat==0,arr.ind=T)
    zeros <- zeros[zeros[,1]!=zeros[,2],]
    zg <- graph.data.frame(zeros,directed=FALSE)
    zc <- data.frame(components(zg)$membership)
    zc[,2] <- as.numeric(rownames(zc))    
    zc <- zc[order(zc[,1],zc[,2]),]
    rownames(zc) <- NULL
    colnames(zc) <- c("group","row")

    dups <- NULL
    for(k in unique(zc$group)) {
        duped <- sort(zc$row[which(zc$group==k)])
        kept <- duped[1]
        duped <- duped[-1]
        dups <- rbind(dups, cbind(rep(kept,length(duped)),duped) )
    }
    dups <- as.matrix(dups[order(dups[,2,drop=F]),])

    out.xy <- df.xy
    for(j in 1:nrow(dups)) {
        out.xy <- insertRow2(out.xy, out.xy[dups[j,1],,drop=F], unname(dups[j,2,drop=F]))
    }

    out.xy
}


# run parallel sparcc
# sparcc.out <- parsparcc(proptab, iter=detectCores()-1, inner_iter=5, th=0.1)
# sparcc.out <- parsparcc(proptab, iter=7, inner_iter=15, th=0.1) #crashed
# for computational sake, only compute correlations on sequences present in at least 5 samples (9472/19073 = 50%)
# crashed after n=21 with cores=13
# n=9 cores=5 took 10 hours
# n=5 cores=5 ran out of memory at i=5
# bonk should use a consistent number of iterations
# n=5 i=5 cores=1 took 43 hours
for(i in seq(11,1,-2)) {
    normtab_small <- normtab[,colSums(normtab>0)>i]
    print(date())
    print(i)
    print(dim(normtab_small))
    sparcc.out <- parsparcc(normtab_small, iter=i, inner_iter=20, th=0.1, cores=11) #crashed
    #sparcc.out <- sparcc(proptab, iter=15, inner_iter=15, th=0.1)
    saveRDS(sparcc.out,file=file.path(outfolder,paste0("sparcc_out_",i,".rds")))






# sparcc.tsne <- Rtsne(props.sparcc$Cor,verbose=TRUE) #,min_cost=1)
# tsne.coord <- tsne.coord.list$Y

# sparcc.out <- readRDS(file=file.path(outfolder,paste0("sparcc_out_",i,".rds")))
# Cluster sequences using sparcc correlations
seq.sparcc.dist <- sparcc.out$Cor
seq.sparcc.dist[is.na(seq.sparcc.dist)] <- 1
seq.sparcc.pca <- fast.prcomp(seq.sparcc.dist)
seq.sparcc.pca.uniq <- unique(seq.sparcc.pca$x)
seq.tsne.sparcc.list <- Rtsne(seq.sparcc.pca.uniq,verbose=TRUE, theta=0.5, pca=FALSE)
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

saveRDS(seq.sparcc.dist,file.path(outfolder,"seq_sparcc_dist.rds"))
saveRDS(seq.sparcc.pca,file.path(outfolder,"seq_sparcc_pca.rds"))
saveRDS(seq.tsne.sparcc.full,file.path(outfolder,"seq_sparcc_tsne.rds"))
saveRDS(seq.tsp.sparcc,file.path(outfolder,"seq_sparcc_tsp.rds"))



#can this be vastly simplified by just saving vectors of weights and colors instead of the whole graph?
network_correlations <- sparcc.out$Cor

# save a graph at each cutoff
save.graphs <- list()
for(i in 100:1) {
    cutoff <- i/100
    print(cutoff)

    # define graph
    network_matrix <- abs(network_correlations) > cutoff
    diag(network_matrix) <- 0
    if(sum(network_matrix)==0) { save.graphs[[as.character(i)]] <- NULL; print("done"); next }

    network_matrix <- as.matrix(tril(network_matrix) + t(tril(network_matrix))) #to make it symmetric so igraph doesn't crash

    network_graph <- adj2igraph(network_matrix)

    save.graphs[[as.character(i)]] <- network_graph
}
saveRDS(save.graphs, file.path(outfolder,"network_graphs.rds")) #730Mb each... must be a better way but for now...

# save the edge weights at each cutoff
save.weights <- list()
for(i in 100:1) {
    print(i)
    network_graph <- save.graphs[[as.character(i)]]
    if(is.null(network_graph)) next
    #define weights and colors
    network_weights <- (network_matrix & upper.tri(network_matrix))*network_correlations

    #scale the weights
    E(network_graph)$weight <- network_weights[abs(network_weights) > cutoff]^3 # cubed preserves sign
    #recenter the weights
    E(network_graph)$weight <- E(network_graph)$weight/max(abs(E(network_graph)$weight))
    #color the rescaled weights
    E(network_graph)$color <- rgb(0,0,0,0)
    network_neg <- E(network_graph)$weight < 0
    network_pos <- E(network_graph)$weight > 0
    E(network_graph)$color[network_neg] <- rgb(abs(E(network_graph)$weight[network_neg]),0,1,abs(E(network_graph)$weight[network_neg]))
    E(network_graph)$color[network_pos] <- rgb(0,0,abs(E(network_graph)$weight[network_pos]),abs(E(network_graph)$weight[network_pos]))
    network_weights[network_weights==1] <- 0

    save.weights[[as.character(i)]] <- network_weights
}
saveRDS(save.weights, file.path(outfolder,"network_weights.rds"))


}
