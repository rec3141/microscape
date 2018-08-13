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


