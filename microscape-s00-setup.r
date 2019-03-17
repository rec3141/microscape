# set up locations and reused functions

options(menu.graphics=FALSE)

#output.path
outfolder <- "out_dada"


# function to normalize subtables and re-merge
norm_esv <- function(st) {
    pt <- prop.table(st,margin=1)
    pt[is.nan(pt)] <- 0
    pt <- pt[!is.na(rownames(pt)),]
    pt
}



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
    dups <- as.matrix(dups[order(dups[,2]),])

    out.xy <- df.xy
    for(j in 1:nrow(dups)) {
        out.xy <- insertRow2(out.xy, out.xy[dups[j,1],], unname(dups[j,2]))
    }

    out.xy
}




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

# 
# # function to add duplicates back in
# # finds zero distances from distance matrix df.dist, adds de-duped rows from df.xy back into df.xy
# addDups <- function(df.dist, df.xy) {
#     #stackoverflow solution to insert row: https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended/11562428#11562428
#     insertRow2 <- function(existingDF, newrow, r) {
#       existingDF <- rbind(existingDF,newrow)
#       existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
#       row.names(existingDF) <- 1:nrow(existingDF)
#       return(existingDF)  
#     }
# 
#     mat <- as.matrix(df.dist)
#     zeros <- which(mat==0,arr.ind=T)
#     zeros <- zeros[zeros[,1]!=zeros[,2],]
#     zg <- graph.data.frame(zeros,directed=FALSE)
#     zc <- data.frame(components(zg)$membership)
#     zc[,2] <- as.numeric(rownames(zc))    
#     zc <- zc[order(zc[,1],zc[,2]),]
#     rownames(zc) <- NULL
#     colnames(zc) <- c("group","row")
# 
#     dups <- NULL
#     for(k in unique(zc$group)) {
#         duped <- sort(zc$row[which(zc$group==k)])
#         kept <- duped[1]
#         duped <- duped[-1]
#         dups <- rbind(dups, cbind(rep(kept,length(duped)),duped) )
#     }
#     dups <- as.matrix(dups[order(dups[,2,drop=F]),])
# 
#     out.xy <- df.xy
#     for(j in 1:nrow(dups)) {
#         out.xy <- insertRow2(out.xy, out.xy[dups[j,1],,drop=F], unname(dups[j,2,drop=F]))
#     }
# 
#     out.xy
# }


#write to fasta format
dada2fasta <- function(esvs,taxanames,taxon) {
    which.taxa <- grepl(taxon,taxanames)
    saved <- paste0(">",taxanames[which.taxa],"\n",colnames(esvs)[which.taxa],"\n")
    cat(saved,sep="")
}
