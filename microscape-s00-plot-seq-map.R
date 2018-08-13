# plot a heatmap of sequence counts on the plate for QC
library(ggplot2)

infile <- read.table("count-demux.txt",sep="\t")

sp1 <- as.data.frame(do.call(rbind,strsplit(x=as.character(infile[,1]),split="_")))
sp2 <- cbind(infile,sp1[,1],substr(sp1[,2], 1, 1),substr(sp1[,2], 2, 3))

for (plate in unique(sp2[,4])) {

print(plate)
mydata <- sp2[which(sp2[,4]==plate),c(5,6,3)]
colnames(mydata) <- c("x","y","z")

#we actually counted lines in FASTQ file
#to get sequences divide by 4
mydata$z <- mydata$z/4

# base graphics
# plate.m <- xtabs(z~x+y, data=mydata)[1:8,]

#try and fail to reverse x factors to match plate
#mydata$x <- factor(mydata$x, levels=(mydata$x)[rev(order(mydata$x))])

#plot absolute
ggplot(mydata, aes(y, x)) +
	ggtitle(plate, subtitle = NULL) +
    geom_tile(aes(fill = z)) + 
    geom_text(aes(label = round(z, 1))) +
    scale_fill_gradient(low = "white", high = "red") 

ggsave(filename=paste0(plate,"-map.pdf"),width=12,height=8)
    
#plot proportion
mydata$z <- 100*mydata$z/sum(mydata$z)

ggplot(mydata, aes(y, x)) +
	ggtitle(plate, subtitle = NULL) +
    geom_tile(aes(fill = z)) + 
    geom_text(aes(label = round(z, 2))) +
    scale_fill_gradient(low = "white", high = "red") 

ggsave(filename=paste0(plate,"-map-prop.pdf"),width=12,height=8)
    
}

