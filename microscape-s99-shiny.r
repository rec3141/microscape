library(shiny)
library(filehash)
library(igraph) # for E
library(ape) #phylogenetics
#library(phytools) #
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())
library(msa)
library(phangorn)
library(gdata) # for read.xls
library(marmap) # for fortify.bathy and plot.bathy
library(ncdf4) # for netCDF
library(vegan) #for metaMDS
library(Rtsne) #for Rtsne
library(spdep) #for Rotation
library(tidyr) 
library(scales)


# run app from /work/cryomics/seq_data/Arctic_Amplicon_clean
# library(shiny)
# runApp(appDir="microscape")

# or work interactively from /work/cryomics/seq_data/Arctic_Amplicon_clean/microscape

# old way
# source("microscape-s99-shiny.r"); shinyApp(ui = ui, server = server)


# deploy
# link needed files into ./microscape dir
# library(rsconnect)

# tags<- rsconnect::appDependencies()
# pkgDep(tags$package[tags$source=="CRAN"], suggests =FALSE, enhances = FALSE)
# dg <- makeDepGraph(tags$package[tags$source=="CRAN"], suggests=F, enhances = F)
# pdf(file="file.pdf",width=72,height=72)
# plot(dg, legendPosition = c(-1, -1), vertex.size = 10, cex = 0.7)
# dev.off()

# options(rsconnect.http = "curl") #got an error using default Rcurl
# options(rsconnect.max.bundle.size=3145728000)
# configureApp("microscape",logLevel="verbose")
# options(repos = BiocInstaller::biocinstallRepos())
# rsconnect::deployApp(logLevel="verbose")

# need a microscape-s90-shiny-setup.r that does all the softlinking and making the shiny directory

# each tab filters the next? 
# no... should they affect each other?
# add a "filter by viewport" checkbox to allow choice
# Step 1) Select samples
# Step 2) Select taxa
# Step 3) Select geography
# Step 4) Select environmental

# allow sparcc correlations on selected subsets?

#### load required data

source("microscape-s00-setup.r")
outfolder <- "."
options(expressions=50000)
system("ulimit -s 16384")

# options(shiny.error = browser)

taxout <- readRDS(file.path(outfolder,"taxout_edit.rds")) #changed from taxout.rds
bootout <- readRDS(file.path(outfolder,"bootout_edit.rds")) #changed from bootout.rds
proptab <- readRDS(file.path(outfolder,"proptab_filt.rds"))
proptab.ord <- t(apply(proptab,1,sort.list,dec=T))

num.taxa <- ncol(proptab)
num.samples <- nrow(proptab)
seqtab <- readRDS(file.path(outfolder,"seqtab_filt.rds"))
normtab <- readRDS(file.path(outfolder,"normtab_filt.rds"))
rs.seqtab <- rowSums(seqtab)
table_list <- readRDS(file.path(outfolder,"table_list.rds"))
primer_list <- unique(table_list)
primer_list <- sort(primer_list[!is.na(primer_list)])

sample.cex <- 1*log10(rs.seqtab)/max(log10(rs.seqtab))
bio.cex.mult <- 10
bio.cex.div <- 1 #inverse
bio.max.taxa <- 10 #num.taxa #maximum overlaid dots, decrease for faster plotting (downsample)

net.cex.mult <- 100 #scaling for network points
net.cex.div <- 1 #inverse
edge.min <- 3
edge.max <- 10

biospatial.px <- 800
network.px <- 800

#i had a note that says this wouldn't work, don't remember why
taxa.list <- readRDS(file.path(outfolder,"taxa_list.rds"))
names.list <- readRDS(file.path(outfolder,"names_list.rds"))


# select which mappings to use for samples
sample.tsne <- readRDS(file.path(outfolder,"sample_bray_tsne.rds"))

# rotate for better view
sample.tsne.rot <- Rotation(sample.tsne,235*pi/180)
tsne.here <- data.frame(sample.tsne.rot)
colnames(tsne.here) <- c("x","y")

# this sets up the sample list from the metadata sheet
metadata.in <- readRDS(file.path(outfolder,"metadata.rds"))
metadata.in <- metadata.in[,c("cellid","desc")]
rownames(metadata.in) <- NULL
metadata <- metadata.in[!duplicated(metadata.in$cellid),]
rownames(metadata) <- metadata$cellid
metadata <- metadata[rownames(proptab),"desc",drop=F]

sample.list <- c("",metadata$desc)

# select which mappings to use for sequences
#choices are differential abundance, or correlation network, but corr net much smaller
seq.tsne <- readRDS(file.path(outfolder,"seq_bray_tsne.rds")) 
seq.tsne <- Rotation(seq.tsne,90*pi/180)

ref <- "ref_dada2_silva_nr_v132_train_set.fasta"
taxlevels <- colnames(taxout[[ref]])

#make default colors
{
		#update when choosing taxon
		leg.taxa <- unlist(lapply(names.list[[ref]],function(x) paste0(strsplit(x,";")[[1]][6:7],collapse=" ")))

		#colors for biospatial plot
		# color.biospatial <- rainbow(4*((num.taxa %/% 4) +1)) #did this so that nearby taxa had different colors to easily distinguish them
		# color.biospatial <- as.character(t(matrix(color.biospatial,nrow=4)))
		# color.biospatial <- color.biospatial[order(newnames)]
		color.biospatial <- rainbow(num.taxa)[order(names.list[[ref]])]

		#default colors for network plot
		color.network <- rep("#427df4",times=num.taxa) #blue default Bacteria
		color.network[grep("Alphaproteobacteria",names.list[[ref]])] <- "#0033ff" #blue
		color.network[grep("Betaproteobacteria",names.list[[ref]])] <- "#1897ff" #blue
		color.network[grep("Gammaproteobacteria",names.list[[ref]])] <- "#1dcdf9" #blue
		color.network[grep("Bacteroidetes",names.list[[ref]])] <- "#0ee0d2" #blue
		color.network[grep("Archaea",names.list[[ref]])] <- "#f9f03e" #yellow
		color.network[grep("Eukaryota",names.list[[ref]])] <- "#f44265" #red
		color.network[grep("Metazoa",names.list[[ref]])] <- "#f463f9" #purple
		color.network[grep("Chloroplast",names.list[[ref]])] <- "#0fe047" #green
		color.network[grep("Ochrophyta",names.list[[ref]])] <- "#0ee047" #green
		color.network[grep("Chlorophyta",names.list[[ref]])] <- "#0de047" #green
		color.network[grep("Syndiniales",names.list[[ref]])] <- "#f91bf5" #purple
		color.network[grep("Mitochondria",names.list[[ref]])] <- "#42e5f4" #orange


}

# if(!any(grepl("save.weights",ls())) & !any(grepl("save.graphs",ls()))) {
# 	save.w <- save.weights[as.character(seq(1,100,10))]
# 	saveRDS(save.w,file.path(outfolder,"network_weights_reduced.rds"))

# 	save.weights <- readRDS(file.path(outfolder,"network_weights_reduced.rds"))
#	save.graphs <- readRDS(file.path(outfolder,"network_graphs.rds"))
# 	save.graphs <- save.graphs[as.character(seq(1,100,10))]

# }

	#newer simpler version
save.weights <- readRDS(file.path(outfolder,"network_melt.rds")) # esv1, esv2, weight, color

# if networking doesn't finish completely, need to limit options
netcutoff <- as.numeric(strsplit(strsplit(list.files(outfolder,pattern="sparcc_out*")[1],"_")[[1]][3],".",fixed=T)[[1]][1])
taxa.in.network <-  unname(colSums(normtab>0)>netcutoff)
taxa.in.network.esvs <- paste0("ESV_",which(taxa.in.network))

# rework to use PhyloSeq
#load(file="seqs.physeq.Rdata")

#        seqs.meta <- sample_data()
#        seqs.physeq <- merge_phyloseq(seqs.physeq, sampledata, seqs.tree)


# for mapping
# 
# # open the netCDF file
# nc <- nc_open("/scratch/SKQ201813S/underway/ETOPO2v2g_f4.nc")
# etopo.in <- ncvar_get(nc, "z") #get depths
# 
# # do some stuff I copied from the web
# colnames(etopo.in) <- ncvar_get(nc, "y")
# rownames(etopo.in) <- ncvar_get(nc, "x")
# class(etopo.in) <- "bathy"



########### START WEBSITE
# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Microscape v0.00003"),

    fluidRow(
    column(
     tabsetPanel(type = "tabs", id= "options", 
      tabPanel("Samples", value="sample-options",
		  selectInput(inputId = "sample", label = "Show single sample", choices = sample.list),
		  textAreaInput(inputId = "plot.sample.regex", label="Show all matching samples (regex)", value=NULL),
		  sliderInput("min.reads", "Minimum sample reads (log2)", min = 0.0, max = 20.0, value = 10.0),
		  sliderInput("min.richness", "Minimum sample richness (log2)", min = 0.0, max = 10.0, value = 1.0),
		  sliderInput("point.biospatial", "Point size (log)", min = -20.0, max = 20.0, value = 0),
		  actionButton(inputId="reset.zoom", label="Reset Zoom"),
		  sliderInput("plot.all.names", label="Show some random names", min=0, max=25, value = 0),
		  selectInput("tax.levels", "P: Taxonomic levels", taxlevels, multiple=TRUE, selected="Genus"),
		  textAreaInput(inputId = "plot.names", label="Show these names (regex)", value=NULL)
      ),
      tabPanel("Taxa Network", value="taxa-options",
		  selectInput(inputId = "ref", label = "Taxonomy Database", choices = names(taxa.list), selected="ref_dada2_silva_nr_v132_train_set.fasta"),
		  selectInput(inputId = "primers", label = "Primer/Target Pairs", choices = primer_list, multiple=TRUE, selected=c("16S_prokaryote", "18S_protist", "ITS_all")),
		  radioButtons(inputId = "primers.log", label = "Primer/Target required", choices = c("AND","OR"), selected = "AND", inline = TRUE),
		  selectInput(inputId = "taxon", label = "Show single taxon", choices = c("",taxa.list[[ref]])),
		  textAreaInput(inputId = "plot.taxa.regex", label="Show all matching taxa (regex)", value=NULL),
		  sliderInput("min.taxa.hits", "Minimum samples per ESV (log2)", min = 0.0, max = 12.0, value = 6),
		  sliderInput("point.network", "Point size (log)", min = -20.0, max = 20.0, value = 0),
			#slider value set to one higher than maximum reported network correlation
	#       sliderInput("network.cor", "N: Minimum SparCC correlation", min = 0.11, max = 1, value = (0.01+length(save.weights)/100), step=0.01),
		  sliderInput("network.cor", "N: Minimum SparCC correlation", min = 0, max = 1, value = 0.71, step=0.01, round=-2),
		  actionButton(inputId="reset.zoom", label="Reset Zoom"),
          sliderInput("plot.all.names", label="Show some random names", min=0, max=25, value = 0),
		  selectInput("tax.levels", "P: Taxonomic levels", taxlevels, multiple=TRUE, selected="Genus"),
		  textAreaInput(inputId = "plot.names", label="Show these names (regex)", value=NULL)
	  ),
      tabPanel("Phylogenetic", value="phylogenetic-options",
		  sliderInput("tree.height", "P: Tree height limits", min = -1000, max = 2000.0, value = c(0,20.0)),
		  sliderInput("tree.width", "P: Tree width limits", min = -500, max = 2000, value = c(0,100)),
		  sliderInput("tree.font.size", "P: Tree font size", min = 0, max = 20, value = 10),
		  sliderInput("tree.colors", "P: Rotate tree colors", min = 0, max = 6, value = c(0,0), step=0.1),
		  selectInput("tax.levels", "P: Taxonomic levels", taxlevels, multiple=TRUE, selected="Genus"),
		  checkboxInput("tax.boots", "P: Taxonomic bootstraps", FALSE)
      )
#       tabPanel("Geographic", value="geographic-options",
# 		  sliderInput("longitude", "M: Longitude", min = -200, max = 200, value = c(-170,-150)),
# 		  sliderInput("latitude", "M: Latitude", min = -90, max = 90, value = c(55,75))
#       )
	), width=3, class="well"
    ),

    column(
     tabsetPanel(type = "tabs", id = "tabs",
        tabPanel("Samples", value="biospatial",
            span(textOutput("samplename"), style="font-size:120%; text-align:center"),
            div(style = "height:800px;", plotOutput("biospatial", hover = "biospatial_hover", click = "biospatial_click", dblclick="biospatial_dblclick"))
        ),
        tabPanel("Taxa Network", value="network",
            span(textOutput("taxaname"), style="font-size:95%; text-align:center"),
            div(style = "height:700px;", plotOutput("network", hover = "network_hover", click = "network_click", dblclick="network_dblclick"))
        ),
        tabPanel("Phylogenetic", value="phylogenetic", 
             div(style = "height:1000px;", plotOutput("phylogenetic", hover = "phylogenetic_hover", click = "phylogenetic_click", dblclick="phylogenetic_dblclick"))
        )
#         tabPanel("Functional", plotOutput("functional")),
#         tabPanel("Geographic", plotOutput("geographic")),
#         tabPanel("Environmental", plotOutput("environmental"))
        

     ),width=6, height=12, class="well"
    ),
	column(
		div(tableOutput("overunder"), style = "font-size:80%"),
		div(tableOutput("rhtable"), style = "font-size:80%")	
	,width=3, class="well")
  )
)



# Define server logic
server <- function(input, output, session) {

    # set up reactive values
    ranges <- reactiveValues(biospatial_xlim=NULL, biospatial_ylim=NULL, network_xlim=c(-60,60), network_ylim=c(-60,60))
    oldxy <- reactiveValues(biospatial_x=0, network_x=0)
    oldzoom <- reactiveValues(z=0)
    picks <- reactiveValues(taxa=NULL,samples=NULL)
    point.size <- reactiveValues(biospatial=0, network=0, input=NULL)
    color <- reactiveValues(phylo=rep(NA,num.taxa))
    last <- reactiveValues(seqs=NULL); isolate(last$seqs)
	taxa <- reactiveValues(list=NULL,names=NULL)

	# keep things up to date
	observe({
		taxa$list <- taxa.list[[input$ref]]
		taxa$names <- names.list[[input$ref]]
	    updateSelectInput(session = session, inputId = "taxon", label = "Show single taxon", choices = c("",taxa$list))
		updateSelectInput(session = session, inputId = "tax.levels", label = "P: Taxonomic levels", 
			choices = colnames(taxout[[input$ref]]))

        # if input$network.cor is too high for selected taxa, lower it
        # shouldn't affect which vertices are plotted, just which edges
        # bring in saved weights
#        print(head(save.weights[[ as.character(input$network.cor*100) ]]))
        network_cor <- input$network.cor

# old        while(is.null(save.weights[[ as.character(network_cor*100) ]] ) & network_cor > 0) {
# can skip this now that we're doing it the network_melt way
#         while(is.null(save.weights[[ as.character(network_cor*100) ]] ) & network_cor > 0) {
#             network_cor <- network_cor - 0.01
# 	        print(network_cor)
#         }
#         if(network_cor == 0) network_cor <- 0.71
#         isolate({updateSliderInput(session, "network.cor", value = network_cor)})
		print(input$network.cor)
	})
	
    # some functions
   xy_biospatial <- function(e) {
        if(is.null(e)) return(NULL)
        new.pos <- c(e$x, e$y) # New position
        # Compute distance to points and select nearest index
        tmp.dist <- colSums((t(sample.tsne.rot) - new.pos)^2)
        tmp.dist[!picks$samples] <- max(tmp.dist)
        nearest.idx <- which.min(tmp.dist)
        nearest.idx
    }
    
   xy_network <- function(e) {
        if(is.null(e)) return(NULL)
        new.pos <- c(e$x, e$y) # New position
        # Compute distance to points and select nearest index
        tmp.dist <- colSums((t(seq.tsne) - new.pos)^2)
        tmp.dist[!picks$taxa] <- max(tmp.dist)
        nearest.idx <- which.min(tmp.dist)
        nearest.idx
    }

    
    ###############
    # HEADERS
    ###############

    output$samplename <- renderText({
        if(is.null(input$biospatial_hover)) {
            as.character("​") #had to use zero width space (U+200B) to avoid wiping
        } else {
            #require it to be picked
            if(picks$samples[xy_biospatial(input$biospatial_hover)]) {
                metadata$desc[xy_biospatial(input$biospatial_hover)]
            } else {
                as.character("​") #had to use zero width space (U+200B) to avoid wiping
            }
        }
    })

    output$taxaname <- renderText({
        if(is.null(input$network_hover)) {
            as.character("​") #had to use zero width space (U+200B) to avoid wiping
        } else {
            if(picks$taxa[xy_network(input$network_hover)]) {
                taxa$names[xy_network(input$network_hover)]
            } else {
                as.character("​") #had to use zero width space (U+200B) to avoid wiping
            }
        }
    })

    ###############
    # biospatial plot
    ###############

      output$biospatial <- renderPlot({

      # MECHANICS AND ZOOM
        if(is.null(ranges$biospatial_xlim)) ranges$biospatial_xlim <- c(-50,50)
        if(is.null(ranges$biospatial_ylim)) ranges$biospatial_ylim <- c(-50,50)
        point.size$biospatial <- input$point.biospatial
        
        if(input$reset.zoom > oldzoom$z & input$tabs=="biospatial") {
            ranges$biospatial_xlim=c(-50,50)
            ranges$biospatial_ylim=c(-50,50)
            oldzoom$z <- input$reset.zoom
        }

        # what to do if user clicks (zoom out) or double clicks (zoom in)
        if(!is.null(input$biospatial_dblclick$x)) {

            if(input$biospatial_dblclick$x != oldxy$biospatial_x) {

            e <- input$biospatial_dblclick
            xlim.diff = 0.4*(ranges$biospatial_xlim[2]-ranges$biospatial_xlim[1])
            ranges$biospatial_xlim <- c(e$x - xlim.diff/2, e$x + xlim.diff/2)
    
            ylim.diff = 0.4*(ranges$biospatial_ylim[2]-ranges$biospatial_ylim[1])
            ranges$biospatial_ylim <- c(e$y - ylim.diff/2, e$y + ylim.diff/2)
    
            oldxy$biospatial_x <- input$biospatial_dblclick$x
            }
        } else if (!is.null(input$biospatial_click$x)) {

            if(input$biospatial_click$x != oldxy$biospatial_x) {
            e <- input$biospatial_click
            xlim.diff = 1.25*(ranges$biospatial_xlim[2]-ranges$biospatial_xlim[1])
            ranges$biospatial_xlim <- c(e$x - xlim.diff/2, e$x + xlim.diff/2)
    
            ylim.diff = 1.25*(ranges$biospatial_ylim[2]-ranges$biospatial_ylim[1])
            ranges$biospatial_ylim <- c(e$y - ylim.diff/2, e$y + ylim.diff/2)
    
            oldxy$biospatial_x <- input$biospatial_click$x
            }
        }


        # PROCESS USER INPUT

        #### SAMPLES FIRST
        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,num.samples)
        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,num.samples)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("biospatial samples: how did you get here?")
        }


		# select subset of samples that have matching primer pairs
		samples.with.primer <- rep(TRUE, num.samples)
		if(input$primers.log == "OR") {
			# only select samples that have ANY user selected primers
			samples.with.primer <- rowSums(proptab[,table_list %in% input$primers])>0
		} else if (input$primers.log == "AND") {
			# only select samples that have ALL user selected primers
			for(primer in input$primers) {
				print(primer)
				samples.with.primer <- samples.with.primer & rowSums(proptab[,table_list %in% primer])>0
			}
		} else {
			print("go away you should not be here!")
		}

        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

        # put all the sample filters together
        # QC and holdup


        map.samples <- min.reads.samples & min.richness.samples & samples.in.box & samples.with.primer

        picks$samples <- map.samples & sample.user.regex

        # print(paste(metadata$desc[picks$samples],sep="\n"))
        pss <- sum(picks$samples)

        if(pss==0) {print("biospatial stopped; no samples selected"); return(0)}

		# PLOT 4: Sample map with circles proportional to relative abundance in n
        # plot the basic outline of selected samples
        plot(sample.tsne.rot[map.samples,], pch=19, cex=sample.cex[map.samples], col="grey", axes=FALSE, ann=FALSE, xlim=ranges$biospatial_xlim, ylim=ranges$biospatial_ylim)

# why?        plot(sample.tsne.rot[map.samples,], pch=19, cex=sample.cex[picks$samples], col="grey", axes=FALSE, ann=FALSE, xlim=ranges$biospatial_xlim, ylim=ranges$biospatial_ylim)

        #### THEN TAXA

        # initialize
        taxa.user.regex <- rep(FALSE,num.taxa)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,num.taxa)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=taxa$names,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=taxa$names, pattern=patt)
            }
        } else {
            print("biospatial taxa: how did you get here?")
        }

		# select subset of taxa that match user selected primer targets
		taxa.with.primer <- table_list %in% input$primers

        # filter by what's shown in taxa box
        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]

        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 
#         picks$taxa <- picks$taxa & taxa_mult_samples #infinite loop?

        # enforce slider limits over taxa choices
#         print(dim(seq.tsne))
#         print(c(length(taxa.user.regex) , length(taxa_mult_samples) , length(taxa.in.box)))
        picks$taxa <- taxa.user.regex & taxa_mult_samples & taxa.in.box & taxa.with.primer
        # print(paste(taxa$list[picks$taxa],sep="\n"))

        pts <- sum(picks$taxa)
        # print(sum(picks$taxa))

        if(pts==0) { print("biospatial stopped; no taxa selected");return(0)} #no valid taxa selected

        # --> diverge

        #only color selected taxa -- not sure if this is necessary
        color.taxa <- rep(NA,num.taxa)
        color.taxa[which(picks$taxa)] <- color.biospatial[which(picks$taxa)]

        # if defined by phylogeny, add some color for plotting points
        color.taxa[which(!is.na(color$phylo))] <- color$phylo[which(!is.na(color$phylo))]

        # for each sample, sort pick.taxa from greatest to least, plot them
        vsize.def <- rep(NA,num.taxa)

		# can i just sort these once and re-use them? yes. proptab.ord
		points.save <- NULL
		vsize.save <- NULL
		color.save <- NULL
        for(n in which(picks$samples)) {
            if(n %% 100 == 0) print("working")
#             seq.ra <- unname(proptab[n,,drop=F])
# #             print(which(picks$taxa))
#             seq.ra[!picks$taxa] <- 0
#             seq.ra[seq.ra==0] <- NA
#             seq.order <- order(seq.ra,decreasing=T,na.last=NA)
# #             print(seq.order)

			#scale point sizes
            vsize <- vsize.def
            vsize <- bio.cex.mult*unname(sqrt(sqrt(proptab[n,]/bio.cex.div)))
            vsize <- vsize*1.1^point.size$biospatial

			#filter
			vsize[vsize==0] <- NA
			vsize[!picks$taxa] <- NA

			#reorder
			seq.order <- proptab.ord[n,]
			
			#trim and plot/ in this order
			vsize.plot <- vsize[seq.order]

			#coords
			num.vals <- sum(!is.na(vsize))
            tsne.plot <- cbind("x"=rep(tsne.here[n,1],num.vals),"y"=rep(tsne.here[n,2],num.vals))
			
			color.plot <- color.taxa[seq.order]
			color.plot <- color.plot[!is.na(vsize.plot)]

			vsize.plot <- vsize.plot[!is.na(vsize.plot)]
			
            points(tsne.plot,pch=19,cex=vsize.plot,col=color.plot)    

#             tsne.plot <- data.frame(tsne.plot)
#             colnames(tsne.plot) <- c("x","y")
    #        tsne.plot.taxa <- rep("",length(pick.taxa))
    #        tsne.plot.taxa[pick.taxa] <- taxa$list[which(pick.taxa)]
#             print(sort(vsize))

# this is much slower
# what if you did just orders and not data
# 			points.save <- rbind(points.save,tsne.plot)
# 			vsize.save <- c(vsize.save, vsize)
# 			color.save <- c(color.save, color.taxa)
        }

#             points(points.save,pch=19,cex=vsize.save,col=color.save)    

#         print(input$plot.all.names)
#         print(input$plot.names)
        
         if(input$plot.all.names > 0 & input$tabs=="biospatial") {
                small.sample <- sample(which(picks$samples),input$plot.all.names, replace=T)
                if(length(small.sample)>0) text(sample.tsne.rot[small.sample,,drop=F],labels=metadata$desc[small.sample])
         }
#                   if(input$plot.all.names & input$tabs=="biospatial") text(sample.tsne.rot[picks$samples,],labels=metadata$desc[picks$samples])

         if(input$plot.names!="") {
            # print(input$plot.sample.names)
            patterns <- strsplit(x=input$plot.names,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                pick.user.samples <- grepl(pattern=patt, x=metadata$desc)
                if(sum(pick.user.samples)>0) text(sample.tsne.rot[pick.user.samples,,drop=F],labels=metadata$desc[pick.user.samples])
            }
        }
        print("done biospatial")
    }, width=biospatial.px,height=biospatial.px)


    ###############
    # network plot
    ###############
    output$network <- renderPlot({
    
    # pre-calculate weighted networks and load
    # to do? click on a taxon to show that taxon's network
    
        ### MECHANICS AND ZOOM
        if(input$reset.zoom > oldzoom$z & input$tabs=="network") {
            ranges$network_xlim=c(-60,60)
            ranges$network_ylim=c(-60,60)
            oldzoom$z <- input$reset.zoom
        }

        point.size$network <- input$point.network
        
        if(!is.null(input$network_dblclick$x)) {

            if(input$network_dblclick$x != oldxy$network_x) {

            e <- input$network_dblclick
            network_xlim.diff = 0.4*(ranges$network_xlim[2]-ranges$network_xlim[1])
            ranges$network_xlim <- c(e$x - network_xlim.diff/2, e$x + network_xlim.diff/2)

            network_ylim.diff = 0.4*(ranges$network_ylim[2]-ranges$network_ylim[1])
            ranges$network_ylim <- c(e$y - network_ylim.diff/2, e$y + network_ylim.diff/2)

            oldxy$network_x <- input$network_dblclick$x
            }
        } else if (!is.null(input$network_click$x)) {

            if(input$network_click$x != oldxy$network_x) {
            e <- input$network_click
            network_xlim.diff = 1.25*(ranges$network_xlim[2]-ranges$network_xlim[1])
            ranges$network_xlim <- c(e$x - network_xlim.diff/2, e$x + network_xlim.diff/2)

            network_ylim.diff = 1.25*(ranges$network_ylim[2]-ranges$network_ylim[1])
            ranges$network_ylim <- c(e$y - network_ylim.diff/2, e$y + network_ylim.diff/2)

            oldxy$network_x <- input$network_click$x
            }
        }

   
        ### PROCESS USER INPUT
        
        ### TAXA FIRST
      
        # initialize
        taxa.user.regex <- rep(FALSE,num.taxa)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,num.taxa)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=taxa$names,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=taxa$names, pattern=patt)
            }
        } else {
            print("network taxa: how did you get here?")
        }
        
		# only use taxa included in primer
		taxa.with.primer <- table_list %in% input$primers

        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 

       ## Plot the basic chart using slider limits
        plot(seq.tsne[taxa_mult_samples,], pch=19, col="grey", cex=0.3, xlim=ranges$network_xlim, ylim=ranges$network_ylim, axes = FALSE, xlab = "", ylab = "")

        # enforce slider limits over taxa choices
        picks$taxa <- taxa.user.regex & taxa_mult_samples & taxa.with.primer #& taxa.in.network 
#         print(paste(taxa$list[picks$taxa],sep="\n"))

        pts <- sum(picks$taxa)
        if(pts==0) {print("network stopped; no taxa selected"); return(0)} #no valid taxa selected



        ### SAMPLES NEXT

        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,num.samples)

        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,num.samples)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("network samples: how did you get here?")
        }

		# select subset of samples that have matching primer pairs
		samples.with.primer <- rep(TRUE, num.samples)
		if(input$primers.log == "OR") {
			# only select samples that have ANY user selected primers
			samples.with.primer <- rowSums(proptab[,table_list %in% input$primers])>0
		} else if (input$primers.log == "AND") {
			# only select samples that have ALL user selected primers
			for(primer in input$primers) {
				print(primer)
				samples.with.primer <- samples.with.primer & rowSums(proptab[,table_list %in% primer])>0
			}
		} else {
			print("go away you should not be here!")
		}

        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

       # put all the sample filters together
        # QC and holdup

        picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex & samples.in.box & samples.with.primer
        # print(paste(metadata$desc[picks$samples],sep="\n"))
        pss <- sum(picks$samples)

         if(pss==0) {print("network stopped; no samples selected");return(0)}

        # add some color for plotting points
        color.taxa <- rep(NA,num.taxa)
        color.taxa[which(picks$taxa)] <- color.network[which(picks$taxa)]
        color.taxa[which(!is.na(color$phylo))] <- color$phylo[which(!is.na(color$phylo))]
#         print(color$phylo[which(!is.na(color$phylo))])
        # --> diverge

        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]

        # scale the point and edge sizes
        network_asize <- unname(proptab[picks$samples,,drop=F])
        if(pss>1) network_asize <- colMeans(network_asize) #collapse
        network_asize <- net.cex.mult*(sqrt(network_asize * 1.8^point.size$network/net.cex.div)) # janky sizing
        network_asize[network_asize==0] <- NA #don't plot zeros
#        network_asize[!networked_taxa] <- NA #don't plot unnetworked
        network_asize[!taxa.in.box] <- NA #don't plot if not in box


# network edge setup

        # select taxa in requested network space: only plot if SparCC correlation passes threshold
        network_save_weights <- save.weights[abs(save.weights$corr) > input$network.cor,,drop=F]


		# only plot edges if there's something to plot
		if(nrow(network_save_weights)>0) {
			# select positive and negative correlations
			# here we use the reduced network size if sparcc isn't run on n=1
		
			networked_esvs <- as.numeric(do.call(rbind,strsplit(x=unique(c(network_save_weights$node1,network_save_weights$node2)),split="_"))[,2])
			networked_taxa <- 1:num.taxa %in% networked_esvs

			# re-enforce limits and merge here 
			networked_taxa <- networked_taxa & taxa_mult_samples
			networked_taxa <- networked_taxa | picks$taxa

			#change weight function
		
			network_save_weights$weight <- network_save_weights$corr/max(abs(network_save_weights$corr))
# 			network_save_weights$color[network_save_weights$corr<0] <- rgb(abs(network_save_weights$weight[network_save_weights$corr<0]),0,0,1)
# 			network_save_weights$color[network_save_weights$corr>=0] <- rgb(0,0,abs(network_save_weights$weight[network_save_weights$corr>=0]),1)
			network_save_weights$color[network_save_weights$corr<0] <- rgb(1,0,0,abs(network_save_weights$weight[network_save_weights$corr<0]))
			network_save_weights$color[network_save_weights$corr>=0] <- rgb(0,0,1,abs(network_save_weights$weight[network_save_weights$corr>=0]))

	#     	network_save_graph$color[network_save_graph$corr<0] <- rgb(abs(network_save_graph$weight[network_save_graph$corr<0]),0,0,abs(network_save_graph$weight[network_save_graph$corr<0]))
	#     	network_save_graph$color[network_save_graph$corr>=0] <- rgb(0,0,abs(network_save_graph$weight[network_save_graph$corr>=0]),abs(network_save_graph$weight[network_save_graph$corr>=0]))

			# get the network from the saved graphs file
	#        network_save_graph <- save.graphs[[as.character(input$network.cor*100)]]

			network_save_graph <- graph_from_data_frame(network_save_weights, directed=F)

			# resize the edges
			network_edgesize <- scales::rescale(E(network_save_graph)$weight,to=c(edge.min,edge.max))

			# add edges and vertices to plot
			if(!all(is.na(network_asize))) {
		
	# old one assumed network_save_graph has all vertices, now it just has those with edges
	#            plot(network_save_graph, layout=seq.tsne[networked_esvs,], xlim=ranges$network_xlim, ylim=ranges$network_ylim, vertex.shape="circle", vertex.size=network_asize, vertex.color=color.taxa, vertex.label=NA, vertex.label.cex=1.5, vertex.label.color="black", edge.label=NA, edge.width=network_edgesize, main=NULL, rescale=F, add=T)
				plot(network_save_graph, layout=seq.tsne[networked_esvs,], xlim=ranges$network_xlim, ylim=ranges$network_ylim, vertex.shape="circle", vertex.size=network_asize[networked_esvs], vertex.color=color.taxa[networked_esvs], vertex.label=NA, vertex.label.cex=1.5, vertex.label.color="black", edge.label=NA, edge.width=network_edgesize, main=NULL, rescale=F, add=T)
			}
		
		}
		
		# plot vertices
        points(seq.tsne[picks$taxa,], pch=21, col="black", bg=color.taxa[picks$taxa], cex=network_asize[picks$taxa], xlim=ranges$network_xlim, ylim=ranges$network_ylim, xlab = "", ylab = "")

        # add random text annotations if requested
         if(input$plot.all.names > 0 & input$tabs=="network") {
            shown_taxa <- taxa.in.box #networked_taxa & 
                small.sample <- sample(which(shown_taxa),input$plot.all.names, replace=T)
                if(length(small.sample)>0) text(seq.tsne[small.sample,,drop=F],labels=taxa$list[small.sample])
         }

        # add specific text annotations if requested 
         if(input$plot.names!="") {
            patterns <- strsplit(x=input$plot.names,split="\n")[[1]]
            for(patt in patterns) {
                if(patt=="\n") next
                pick.user.taxa <- grepl(pattern=patt, x=taxa$list)
                if(sum(pick.user.taxa)>0) text(seq.tsne[pick.user.taxa,,drop=F],labels=sapply(which(pick.user.taxa), function(x) paste(strsplit(taxa$names[x],split=";")[[1]][c(6,7)],collapse=" ")))
            }
        }
        print("done network")
    }, width=network.px,height=network.px)







    ###############
    # geographic plot
    ###############
    output$geographic <- renderPlot({
    
    	h2("Coming Soon!")
    	
        # set up the map range for ggplot
        range.lat <- input$latitude
        range.lon <- input$longitude

        print(range.lat)
        print(range.lon)
        # you can change the limits of the base map here
        lat.s <- which(as.numeric(colnames(etopo.in))==range.lat[1])
        lat.n <- which(as.numeric(colnames(etopo.in))==range.lat[2])

        lon.e <- which(as.numeric(rownames(etopo.in))==range.lon[1])
        lon.w <- which(as.numeric(rownames(etopo.in))==range.lon[2])

        etopo.chukchi <- etopo.in[lon.e:lon.w, lat.s:lat.n]
        class(etopo.chukchi) <- "bathy"

        # this reformats the data into a ggplot-friendly format
        chukchi.df <- fortify(etopo.chukchi)

        # we'll plot the map using ggplot2
        # ggplot is a layer-based plotting system, so we'll set up the layers one by one and then combine them at the end

        # prepare the base topographic map

        # each line adds a new layer of data or information, make sure to include the '+' at the end of each line before the final
        # the n=10 tells ggplot how many levels to use to interpolate the colors
        m.base <- ggplot(chukchi.df, aes(x=x, y=y)) +
            geom_raster(aes(fill=z), data=chukchi.df) + scale_fill_etopo(n=10) + # the raster image of the bathymetry
            scale_x_continuous(limits=c(range.lon)) + scale_y_continuous(limits=c(range.lat)) # limit the axes to the map range given previously

        # to more specifically focus on the ocean, we'll download a nice bathymetric color scheme from the shared drive, "scale_fill_etopo_bathy"
        source("http://share.sikuliaq.alaska.edu/public/Cruises/SKQ201813S/ObservationalOceanography/Eric/etopo_colors.r")

        # the limits tell ggplot which data to plot, here we're greying out everything below 3000 m and above 0 m
        m.base <- ggplot(chukchi.df, aes(x=x, y=y)) +
            geom_raster(aes(fill=z), data=chukchi.df) + scale_fill_etopo_bathy(n=10,limits=c(-3000,0)) + # the raster image of the bathymetry
            scale_x_continuous(limits=c(range.lon)) + scale_y_continuous(limits=c(range.lat)) # limit the axes to the map range given previously

        # you can add more layers as necessary, here we'll add some contour lines
        m.fancy <- m.base +
            geom_contour(aes(z=z),breaks=c(-10, -20, -30, -40, -50, -60, -70, -80, -90, -100),colour="white", size=0.2) + # contour lines for shallow depths
            geom_contour(aes(z=z),breaks=c(-200, -250, -300, -350, -400, -450),colour="grey", size=0.2) + # contour lines for medium depths
            geom_contour(aes(z=z),breaks=c(-500, -1000, -1500, -2000, -2500, -3000, -3500, -4000),colour="blue", size=0.2) + # contour lines for deep depths
            geom_contour(aes(z=z),breaks=c(0),colour="black", size=0.4) # contour lines for the shoreline

        # we can also read in some waypoints from previous cruises
#         waypoint.in <- read.xls(xls="http://share.sikuliaq.alaska.edu/public/Cruises/SKQ201813S/ObservationalOceanography/Eric/Station_Plan.xls",sheet="Plan")
        waypoint.data <- waypoint.in[,1:4]
        waypoint.data <- waypoint.data[complete.cases(waypoint.data),]

        # give them some styling
        # for points: small black circles
        l.waypoint.point <- geom_point(data=waypoint.data, aes(x=longitude, y=latitude), size=0.2, color="black", fill="white")

        # and add them to the map
        m.fancy <- m.fancy + l.waypoint.point

        # you can even add labels, but it gets a bit messy
        # styling for labels: small text, avoiding overlaps, opaque, and nudged to the right of the points 
        l.waypoint.label <-	geom_text(data=waypoint.data, aes(x=longitude, y=latitude, label=site), size=2, check_overlap=T, alpha="1", nudge_x=0.01*(abs(range.lat[1]-range.lat[2])))

        # and add them to the map
        m.fancy <- m.fancy + l.waypoint.label

        # if we want, we can remove the depth legend as follows
        m.fancy <- m.fancy + theme(legend.position="none")

        # and add a title and labels
        m.fancy <- m.fancy + ggtitle("Temperature") + xlab("Longitude") + ylab("Latitude")

        # finally, read in underway data from the Sikuliaq LDS by its URL
        data.in <- read.table("http://data.sikuliaq.alaska.edu/archive/SKQ201812T/lds/proc/20180601Z/SKQ201812T_underway_20180601Z.txt",comment="%")

        # name the columns
        data.names <- c(
            "Year",
            "Month",
            "Day",
            "Hour",
            "Minute",
            "Second",
            "ins_seapath_position Latitude",
            "ins_seapath_position Longitude",
            "ins_seapath_position Speed over ground",
            "ins_seapath_position Course over Ground",
            "ins_seapath_position Heading",
            "ins_seapath_position Roll",
            "ins_seapath_position Pitch",
            "ins_seapath_position Heading",
            "ins_seapath_position Heave",
            "mb_em302_centerbeam Latitude",
            "mb_em302_centerbeam Longitude",
            "mb_em302_centerbeam Depth",
            "met_ptu307 Atmospheric pressure",
            "met_ptu307 Air Temperature",
            "met_ptu307 Relative Humidity",
            "fluoro_turner-c6 Phycoeryth",
            "fluoro_turner-c6 CDOM",
            "fluoro_turner-c6 Chlorophyll_a",
            "fluoro_turner-c6 Crude Oil",
            "fluoro_turner-c6 Turbidity",
            "fluoro_turner-c6 Depth",
            "fluoro_turner-c6 Temp",
            "wind_gill_fwdmast_true Wind direction",
            "wind_gill_fwdmast_true Wind speed N",
            "rad_psp-pir LW",
            "rad_psp-pir SW",
            "rad_qsr2150a PAR",
            "flow_krohne_fwd flow speed",
            "flow_krohne_fwd volume flow",
            "flow_krohne_fwd coil temperature",
            "flow_krohne_fwd conductivity",
            "tsg_sbe45_fwd Temperature",
            "tsg_sbe45_fwd Conductivity",
            "tsg_sbe45_fwd Salinity",
            "tsg_sbe45_fwd Speed of sound",
            "thermo_sbe38_fwd Intake Temperature",
            "oxygen_sbe43_fwd_calc Dissolved oxygen",
            "oxygen_sbe43_fwd_calc Temperature",
            "oxygen_sbe43_fwd_calc Salinity",
            "oxygen_sbe43_fwd_calc Oxygen Solubility",
            "oxygen_sbe43_fwd_calc Dissolved Oxygen"
        )

        #make unique names with no spaces
        colnames(data.in) <- make.names(data.names)
        colnames(data.in)

        # add temperature track to fancy map
        # the way we'll do this is to plot the lat/long of our cruise track, and use colors to show the temperature
        # we'll get the lat/lon from the EM302
        # first we need to get some colors, let's do a quick and dirty version
        # since we know all of the temperature data is going to fall between 0 and 20 C, we'll use a built-in color ramp from 1:20
        twenty.colors <- heat.colors(20) #or topo.colors(), terrain.colors(), rainbow(), etc.
        simple.colors <- twenty.colors[round(data.in$tsg_sbe45_fwd.Temperature)]

        m.temperature <- m.fancy +
            geom_point(data=data.in, aes(x=mb_em302_centerbeam.Longitude, y=mb_em302_centerbeam.Latitude), color=simple.colors, size=2)

        # here's one way to make a nicer color gradient
        col.data <- data.in$tsg_sbe45_fwd.Temperature
        range.temperature <- range(col.data)
        col.row <- round(scales::rescale(col.data,to=c(1,nrow(data.in)), from=range.temperature))
        col.row[col.data < range.temperature[1]] <- 1
        col.row[col.data > range.temperature[2]] <- nrow(data.in)
        mapcolors <- rev(topo.colors(nrow(data.in)))[col.row]

        m.temperature <- m.fancy +
            geom_point(data=data.in, aes(x=mb_em302_centerbeam.Longitude, y=mb_em302_centerbeam.Latitude), color=mapcolors, size=2)

        print(m.temperature)

        # and if you want to plot the temperature against time in ggplot
        # first convert the time
#         data.in$gmt_time <- as.POSIXct(strptime(paste(data.in$Year,data.in$Month,data.in$Day,data.in$Hour,data.in$Minute,data.in$Second),format="%Y %m %d %H %M %S"))
#         data.in$local_time <- data.in$gmt_time -8*60*60
# 
#         p.temperature <- ggplot() +
#             geom_point(data=data.in, aes(x=local_time, y=tsg_sbe45_fwd.Temperature), col="grey") + #biochem lab
#             geom_point(data=data.in, aes(x=local_time, y=tsg_sbe45_fwd.Temperature), col=mapcolors) + #portside 
#             xlab("local time") + ylab("Temperature") + ylim(range.temperature)
    
print("done geographic")
    
    }, width=700, height=700)

    ###############
    # environmental plot
    ###############
    output$environmental <- renderPlot({
    
    h2("Coming Soon!")
    
    })

    ###############
    # phylogenetic plot
    ###############
    output$phylogenetic <- renderPlot({
    
        ### tree annotated with sample heatmap
        ### PROCESS USER INPUT
        
        #### TAXA FIRST

		#filter for those in network

        taxa.user.regex <- rep(FALSE,num.taxa)

        # if user doesn't provide a taxon or a regex, select all
        # print(input$taxon)
        if(input$taxon=="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- rep(TRUE,num.taxa)
        # else if user provides a lone taxon and no regex, grep it
        } else if (input$taxon!="" & input$plot.taxa.regex=="") {
            taxa.user.regex <- grepl(x=taxa$names,pattern=input$taxon)
        } else if(input$plot.taxa.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
            patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                taxa.user.regex <- taxa.user.regex | grepl(x=taxa$names, pattern=patt)
            }
        } else {
            print("phylogenetic taxa: how did you get here?")
        }

		# only use taxa included in primer
		taxa.with.primer <- table_list %in% input$primers

        # filter by what's shown in taxa box
        taxa.in.box <- seq.tsne[,1] > ranges$network_xlim[1] & 
        seq.tsne[,1] < ranges$network_xlim[2] & 
        seq.tsne[,2] > ranges$network_ylim[1] & 
        seq.tsne[,2] < ranges$network_ylim[2]

        #select taxa in multiple samples
        taxa_mult_samples <- colSums(proptab>0) >= 2^input$min.taxa.hits 
#         picks$taxa <- picks$taxa & taxa_mult_samples #infinite loop?

        # enforce slider limits over taxa choices
        print(c(length(taxa.user.regex) , length(taxa_mult_samples) , length(taxa.in.box)))

        picks$taxa <- taxa.user.regex & taxa_mult_samples & taxa.in.box & taxa.with.primer
        # print(paste(taxa$list[picks$taxa],sep="\n"))

        pts <- sum(picks$taxa)
        if(pts<3) {print("phylogenetic stopped"); plot(0,0,main="not enough sequences",axes = FALSE, xlab = "", ylab = "", pch=""); return(0)} #can't do phylogeny on <3 sequences
        


        #### THEN SAMPLES
        # Sample filters
        min.reads.samples <- rs.seqtab > 2^input$min.reads
        min.richness.samples <- rowSums(proptab>0) >= 2^input$min.richness

        # initialize 
        sample.user.regex <- rep(FALSE,num.samples)
        # print(input$sample)
        # if user doesn't provide a lone sample or a regex, select all
        if(input$sample=="" & input$plot.sample.regex=="") {
            sample.user.regex <- rep(TRUE,num.samples)
        # if user provides a lone sample and no regex, grep for it
        } else if (input$sample!="" & input$plot.sample.regex=="") {
            sample.user.regex <- metadata$desc==input$sample
        } else if(input$plot.sample.regex!="") {
        # if user provides a list of regex, split and search each one, ORing
        # override input$sample
            patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
            for(patt in patterns) {
                # print(patt)
                if(patt=="\n") next
                sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
            }
        } else {
            print("phylogenetic samples: how did you get here?")
        }

		# select subset of samples that have matching primer pairs
		samples.with.primer <- rep(TRUE, num.samples)
		if(input$primers.log == "OR") {
			# only select samples that have ANY user selected primers
			samples.with.primer <- rowSums(proptab[,table_list %in% input$primers])>0
		} else if (input$primers.log == "AND") {
			# only select samples that have ALL user selected primers
			for(primer in input$primers) {
				print(primer)
				samples.with.primer <- samples.with.primer & rowSums(proptab[,table_list %in% primer])>0
			}
		} else {
			print("go away you should not be here!")
		}

        # don't waste time calculating and plotting things that aren't in the view box
        samples.in.box <- sample.tsne.rot[,1] > ranges$biospatial_xlim[1] & 
        sample.tsne.rot[,1] < ranges$biospatial_xlim[2] & 
        sample.tsne.rot[,2] > ranges$biospatial_ylim[1] & 
        sample.tsne.rot[,2] < ranges$biospatial_ylim[2]

        # put all the sample filters together
        # QC and holdup

        picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex & samples.in.box & samples.with.primer
        # print(paste(metadata$desc[picks$samples],sep="\n"))
        pss <- sum(picks$samples)
        if(pss==0) {print("phylogenetic stopped"); return(0)}


        # keep it small for now
        if(pts > 60) {
            #want selected taxa in the top 100 in abundance for the selected 
            seqs.taxa <- rep(FALSE,num.taxa) # F F F F
            order.taxa <- order(colSums(proptab[picks$samples,picks$taxa]),decreasing=T) # 3 1 2 4
            seqs.taxa[picks$taxa] <- order.taxa #probably FALSE got converted to zeros
          
            seqs.taxa <- seqs.taxa < 100 & seqs.taxa > 0
#             order.taxa.pick <- order.taxa[picks$taxa] # 1 4 
#             order.taxa.pick.order <- order(order.taxa.pick)  # 1 2 
#             order.taxa.pick.order.which <- which(order.taxa.pick.order <=100) # 1 2
#             order.taxa.pick.order.which.which <- which(
#             order.max <- max(order.taxa[picks$taxa])         
#             seqs.taxa[order.taxa <= order.max & picks$taxa] <- TRUE
#             seqs.taxa[sample(which(picks$taxa),100)] <- TRUE
            print(paste0(pts-sum(seqs.taxa)," sequences not shown"))
            
        } else {
            seqs.taxa <- picks$taxa
        }
#         sub.physeq <- subset_taxa(seqs.physeq, seqs.taxa)
#         sub.physeq <- subset_samples(sub.physeq, picks$samples)
#        sub.bin <- seqs.bin[seqs.taxa]
#        checkAlignment(sub.bin)

        seqs.ss <- DNAStringSet(colnames(proptab))
        sub.ss <- seqs.ss[seqs.taxa,]
        # align sequences / will fail miserably with ITS
#        seqs.aln <- muscle(sub.bin, quiet=FALSE, MoreArgs="-maxiters 1 -diags")
		sub.aln <- msa(sub.ss,type="dna","ClustalOmega",order="input")
        sub.aln.c <- as.character(sub.aln@unmasked)
        sub.aln.raw <- t(sapply(strsplit(sub.aln.c,""), tolower))
        rownames(sub.aln.raw) <- colnames(proptab)[seqs.taxa]
        sub.bin <- as.DNAbin(sub.aln.raw)
        # get distances
        sub.dist <- dist.dna(sub.bin,as.matrix=T)
        sub.dist[is.na(sub.dist)] <- 1
        sub.dist[is.infinite(sub.dist)] <- 1
        sub.dist[is.nan(sub.dist)] <- 1
        # make tree
        sub.tree <- bionjs(sub.dist)
        sub.tree <- midpoint(sub.tree)
        
        #assign taxonomic labels
        user.levels <- input$tax.levels
        print(user.levels)
		print(colnames(taxout[[input$ref]]))
        if(length(user.levels)==0) { 
            sub.tree$tip.label <- paste0("ESV_",which(seqs.taxa))
        } else {
            user.taxlevels <- unname(taxout[[input$ref]][colnames(proptab)[seqs.taxa],user.levels])
            user.taxlevels <- cbind(user.taxlevels,paste0("ESV_",which(seqs.taxa)))
            sub.tree$tip.label <- apply(user.taxlevels,1,paste,collapse=";")
        }
        
        if(input$tax.boots) {
        	sub.tree$tip.label <- paste0(sub.tree$tip.label," (",bootout[[input$ref]][colnames(proptab)[seqs.taxa],user.levels[length(user.levels)]],")")        
        }

#         sub.physeq <- subset_taxa(seqs.physeq, seqs.taxa)
#         sub.physeq <- subset_samples(sub.physeq, picks$samples)
#         sub.otu <- otu_table(t(proptab[,seqs.taxa]), taxa_are_rows = TRUE)
#         sub.tax <- tax_table(taxout$ref_dada2_silva_nr_v132_train_set.fasta$tax[colnames(proptab)[seqs.taxa],])
#         sub.physeq <- phyloseq(sub.otu, sub.tax, sub.tree)
    
        #scale to relative abundances
        if(pss>1) {
            sub.edge.weight <- unname(colSums(proptab[picks$samples,seqs.taxa,drop=F]))
        } else {
            sub.edge.weight <- unname(proptab[picks$samples,seqs.taxa,drop=F])
        }

        # now need to calculate for internal nodes
        sub.node.weight <- rep(0,nrow(sub.tree$edge)+1)
        for(tip in 1:length(sub.tree$tip.label)) {
            tip.anc <- c(tip,Ancestors(sub.tree,type="all",node=tip))
#             sub.node.weight[tip] <- sub.edge.weight[tip]
#             sub.node.weight[tip.anc] <- sub.node.weight[tip.anc] + sub.edge.weight[tip]
            for(t in 2:length(tip.anc)) {
                tt <- which(sub.tree$edge[,1]==tip.anc[t] & sub.tree$edge[,2]==tip.anc[t-1])
#                 print(tt)
                sub.node.weight[tt] <- sub.node.weight[tt] + sub.edge.weight[tip]
            }
        }
        sub.node.weight <- scales::rescale(sqrt(sqrt(sub.node.weight)),to=c(1,10))
#        sub.node.weight <- sub.node.weight[match(1:length(sub.node.weight),sub.tree$edge[,1])]
        # new distance matrix from branch lengths
        sub.dist.tree <- cophenetic.phylo(sub.tree)
        # jitter to avoid Rtsne complaining
        sub.dist.jit <- apply(sub.dist.tree,1,jitter)
        # sometimes Rtsne fails "perplexity too high"
        tsne.complete <- 0
        perp <- round(ncol(sub.dist.jit)/3)
        counter <- 0
        while(tsne.complete==0) {
            counter <- counter+1
            if(counter>100) return(NULL)
            try({
                perp <- perp-1
               set.seed(42) #I don't actually understand why this works to prevent infinite loop
                sub.tsne <- unname(metaMDS(comm=sub.dist, k=3, distance="euclidean")$points)
#                sub.tsne <- Rtsne(X=sub.dist.jit,dims=3,perplexity=perp)$Y
                tsne.complete <- 1
            })
        }

#         print(color$phylo[seqs.taxa])

        # if colors for any member of this set of taxa are not set, and no rotation is requested, make default colors
        if(any(is.na(color$phylo[seqs.taxa])) | all(input$tree.colors==0)) {
            sub.colors <- cbind(sub.tsne[,1],sub.tsne[,2],sub.tsne[,3])
            sub.colors.rot <- rgb(scales::rescale(sub.colors,to=c(0,0.8)))

#             print("default colors")

        # allow reset to default colors
        
        # else use rotated color

        } else if (any(input$tree.colors!=0)) {
            sub.colors <- cbind(sub.tsne[,1],sub.tsne[,2],sub.tsne[,3])
            sub.colors.rot <- cbind(Rotation(sub.colors[,1:2],input$tree.colors[1]),Rotation(sub.colors[,c(2,3)],input$tree.colors[2]))
            sub.colors.rot <- rgb(scales::rescale(sub.colors.rot,to=c(0,0.8)))
#             print("new colors")
        }

        #kill infinite loop
        isolate({color$phylo[seqs.taxa] <- sub.colors.rot})

        # plot tree, auto-calculate scaling only if first time
        if(is.null(last$seqs) | !identical(last$seqs,picks$taxa)) { #if defaults
            sub.plot <- plot(sub.tree, no.margin=TRUE, plot=FALSE, edge.width=1, label.offset=input$tree.width[2]/50000)
            updateSliderInput(session, "tree.height", value = sub.plot$y.lim*10)
            updateSliderInput(session, "tree.width", value = sub.plot$x.lim*1000)
        }
       
       plot(sub.tree, no.margin=TRUE, x.lim=input$tree.width/1000,y.lim=input$tree.height/10, cex=input$tree.font.size/10, tip.color=sub.colors.rot, edge.width=sub.node.weight, show.node.label = FALSE, label.offset=input$tree.width[2]/50000)
        #tiplabels() for relative abundance data
       isolate({last$seqs <- picks$taxa})

        print("done phylogenetic")
    }, width=1000, height=900)


    ###############
    # functional plot
    ###############
    output$functional <- renderPlot({
    
    h2("Coming Soon!")
    
    })

    ###############
    # taxa table
    ###############
      output$rhtable <- renderTable({
        if(input$tabs == "biospatial") {
            #we want a table of selected taxa in the sample being hovered over

            xy <- xy_biospatial(input$biospatial_hover)
            if(is.null(xy)) return()
            # print(xy)
            if(picks$samples[xy]) {
                sample.count <- unname(seqtab[xy,picks$taxa,drop=F])
                sample.sum <- sum(unname(seqtab[xy,,drop=F]))
                rhtable.num <- round(1000*sample.count/sample.sum,1)[order(sample.count,decreasing=T)] #bc of sprintf
                rhtable.ppt <- sprintf('%.2f',rhtable.num)
                rhtable.taxa <- taxa$names[picks$taxa][order(sample.count,decreasing=T)]
                  data.frame(
                    cbind(
                        "per mille"=rhtable.ppt[rhtable.num>0],                  
                        "taxon"=rhtable.taxa[rhtable.num>0]
                        )
                    )
            }
        } else if(input$tabs == "network") {
            #we want a table of networked taxa connected to the taxa being hovered over
            # or do we want a table of samples that the hovered taxa is in
            xy <- xy_network(input$network_hover) #gives closest point
            if(is.null(xy)) return()
            if(picks$taxa[xy]) {
                taxa.count <- unname(proptab[picks$samples,xy,drop=F]) #proportional reads of taxon in selected samples
            	taxa.order <- order(taxa.count,decreasing=T)
                # taxa.sums <- rowSums(unname(proptab[picks$samples,,drop=F])) #proportional reads of all taxa in selected samples
                rhtable.num <- round(taxa.count*1000,1)[taxa.order] #bc of sprintf
                rhtable.ppt <- sprintf('%.2f',rhtable.num)
                rhtable.taxa <- metadata$desc[picks$samples][taxa.order]
                  data.frame(
                    cbind(
                        "per mille"=rhtable.ppt[rhtable.num>0],                  
                        "sample"=rhtable.taxa[rhtable.num>0]
                        )
                    )
            }
        } else {
            print("nothing to see here")
#            data.frame()
        }
      }, align='rl',digits=2)

# 
# 
#     ###############
#     # overunder table
#     ###############
# 
#       output$overunder <- renderTable({
#         if(input$tabs == "biospatial") {
#             #we want a table of over and underrepresented taxa in the sample being hovered over
# 
#             xy <- xy_biospatial(input$biospatial_hover)
#             # if not hovering, then table is for whole viewport
#             if(is.null(xy)) {
#                 
#             
#             } else {
#             # if hovering then table is just for that sample
#                 if(picks$samples[xy]) {
#                     sample.count <- unname(seqtab[xy,picks$taxa,drop=F])
#                     sample.sum <- sum(unname(seqtab[xy,,drop=F]))
#                     overunder.num <- round(1000*sample.count/sample.sum,1)[order(sample.count,decreasing=T)] #bc of sprintf
#                     overunder.ppt <- sprintf('%.2f',overunder.num)
#                     overunder.taxa <- taxa$list[picks$taxa][order(sample.count,decreasing=T)]
#                       data.frame(
#                         cbind(
#                             "per mille"=overunder.ppt[overunder.num>0],                  
#                             "taxa"=overunder.taxa[overunder.num>0]
#                             )
#                         )
#                 }
#             }
#         } else if(input$tabs == "network") {
#             # we want a table of samples which have enriched abundances of the taxa in the FOV
#             xy <- xy_network(input$network_hover) #gives closest point
#             if(is.null(xy)) return()
#             if(picks$taxa[xy]) {
#                 taxa.count <- unname(proptab[picks$samples,xy,drop=F]) #raw reads of taxon in selected samples
#                 taxa.sums <- rowSums(unname(proptab[picks$samples,,drop=F])) #raw reads of all taxa in selected samples
#                 overunder.num <- round(taxa.count/1000,1)[order(taxa.count,decreasing=T)] #bc of sprintf
#                 overunder.ppt <- sprintf('%.2f',overunder.num)
#                 overunder.taxa <- metadata$desc[picks$samples][order(taxa.count,decreasing=T)]
#                   data.frame(
#                     cbind(
#                         "per mille"=overunder.ppt[overunder.num>0],                  
#                         "samples"=overunder.taxa[overunder.num>0]
#                         )
#                     )
#             }
#         } else {
#             print("nothing to see here")
#             data.frame()
#         }
#       }, align='rl',digits=2)
#     
#     


}

# Create Shiny app ----
shinyApp(ui = ui, server = server)



#network
# do this once outside of shiny app and import it

# can this be vastly simplified by just saving vectors of weights and colors instead of the whole graph?
#network_correlations <- props.sparcc$Cor #does it make sense to try to choose networked taxa here?
#save.graphs <- list()
#save.weights <- list()

# for(i in 100:1) {
#     cutoff <- i/100
#     print(cutoff)
#     # define graph
#     network_matrix <- abs(network_correlations) > cutoff
#     diag(network_matrix) <- 0
#     if(sum(network_matrix)==0) { save.graphs[[i]] <- NULL; save.weights[[i]] <- NULL; next }
#     network_matrix <- as.matrix(tril(network_matrix) + t(tril(network_matrix))) #to make it symmetricso igraph doesn't crash
# #    network_matrix <- Matrix(network_matrix, sparse=TRUE) #does this actually speed anything up?
#     network_graph <- adj2igraph(network_matrix)
#     #define weights and colors
#     network_weights <- (network_matrix & upper.tri(network_matrix))*network_correlations
#     #scale the weights to x^4
#     E(network_graph)$weight <- network_weights[abs(network_weights) > cutoff]^3 #preserve sign
#     #recenter the weights
#     E(network_graph)$weight <- E(network_graph)$weight/max(abs(E(network_graph)$weight))
#     #color the rescaled weights
#     E(network_graph)$color <- rgb(0,0,0,0)
#     network_neg <- E(network_graph)$weight < 0
#     network_pos <- E(network_graph)$weight > 0
#     E(network_graph)$color[network_neg] <- rgb(abs(E(network_graph)$weight[network_neg]),0,1,abs(E(network_graph)$weight[network_neg]))
#     E(network_graph)$color[network_pos] <- rgb(0,0,abs(E(network_graph)$weight[network_pos]),abs(E(network_graph)$weight[network_pos]))
#     network_weights[network_weights==1] <- 0
#     save.graphs[[i]] <- network_graph
#     save.weights[[i]] <- network_weights
# }
# save(file="save.graphs.Rdata",save.graphs) #730Mb each... must be a better way but for now...
# save(file="save.weights.Rdata",save.weights)



#         print(input$point.scalar)
#         print(point.size$network)
#         if(input$point.scalar != 0 & input$point.scalar != point.size$network) {
#             point.size$network <- point.size$network + input$point.scalar
#             updateSliderInput(session, "point.scalar", value = 0)
#         }

#  old
# 
#         min.reads.samples <- rs.seqtab > 2^input$min.reads
#         min.richness.samples <- rowSums(proptab>0) > 2^input$min.richness
# 
#         sample.user.regex <- metadata$desc==input$sample
# #        sample.user.regex <- rep(FALSE,num.samples)
# 
#         if(input$sample=="") {
#             sample.user.regex <- rep(TRUE,num.samples)
#         } else {
#         sample.user.regex <- metadata$desc==input$sample
# 
#             if(input$plot.sample.regex != "") {
#                 patterns <- strsplit(x=input$plot.sample.regex,split="\n")[[1]]
#                 for(patt in patterns) {
#                     # print(patt)
#                     if(patt=="\n") next
#                     sample.user.regex <- sample.user.regex | grepl(x=metadata$desc, pattern=patt)
#                 }
#             }
#         }
# 
# 
#         picks$samples <- min.reads.samples & min.richness.samples & sample.user.regex
#         pss <- sum(picks$samples)
#         
#         

#
# old trash if works 
#     #choose networked taxa to plot
#     taxon <- input$taxon
#     
#     picks$taxa <- rep(FALSE,num.taxa)
#     taxa.user.regex <- rep(FALSE,num.taxa)
#     if(input$plot.taxa.regex != "") {
#         patterns <- strsplit(x=input$plot.taxa.regex,split="\n")[[1]]
#         for(patt in patterns) {
#             # print(patt)
#             if(patt=="\n") next
#             taxa.user.regex <- taxa.user.regex | grepl(x=taxa$list, pattern=patt)
#         }
#         picks$taxa <- taxa.user.regex #logical
#     } else {
#         picks$taxa <- grepl(pattern=taxon,x=taxa$list,fixed=T) #logical
#     }
#     


# can't figure out how to not make it scale both graphs
#         print(input$point.scalar)
#         print(input$tabs)
#         if(input$point.scalar != 0 & input$point.scalar != point.size$biospatial) {
#             point.size$biospatial <- point.size$biospatial + input$point.scalar
#             updateSliderInput(session, "point.scalar", value = 0)
#         }

# forget autoscaling again
# 
#         print(networked_taxa)
#         if(is.null(ranges$network_xlim)) {
# #             xlim <- c(-60,60)
#             print(networked_taxa)
#             print(seq.tsne[networked_taxa,])
#             xlim2 <- range(seq.tsne[networked_taxa,1,drop=F],na.rm=T,finite=T)
# #             ox <- mean(seq.tsne[networked_taxa,1,drop=F],na.rm=T,finite=T)
# #             xd <- max(abs(xlim2-ox))
# #             zd <- c(xd,yd)[which.max(abs(c(xd,yd)))]
# #             xlim2[1] <- ox - xd*1.1
# #             xlim2[2] <- ox + xd*1.1
# #             if(xlim2[1] < xlim[1]) xlim2[1] <- xlim[1]
# #             if(xlim2[2] > xlim[2]) xlim2[2] <- xlim[2]
#             ranges$network_xlim <- xlim2
#         }
# 
#         if(is.null(ranges$network_ylim)) {
# #            ylim <- c(-60,60)
#             ylim2 <- range(seq.tsne[networked_taxa,2,drop=F])
# #             oy <- mean(seq.tsne[networked_taxa,2,drop=F])
# #             yd <- max(abs(ylim2-oy))
# #             zd <- c(xd,yd)[which.max(abs(c(xd,yd)))]
# #             ylim2[1] <- oy - yd*1.1
# #             ylim2[2] <- oy + yd*1.1
# #             if(ylim2[1] < ylim[1]) ylim2[1] <- ylim[1]
# #             if(ylim2[2] > ylim[2]) ylim2[2] <- ylim[2]
#             ranges$network_ylim <- ylim2
#         }
# 
#         print(paste0("ranges:", ranges$network_xlim))
# 
# 
# library(shiny)
# runApp( list(ui = bootstrapPage(
#   verbatimTextOutput("results"),
#   tags$script('
#     $(document).on("keypress", function (e) {
#        Shiny.onInputChange("mydata", e.which);
#     });
#   ') 
# )
# , server = function(input, output, session) {
# 
#   output$results = renderPrint({
#     input$mydata
#   })
# }
# ))
# for keydown events you can substitute:
# 
#   tags$script('
#     $(document).on("keydown", function (e) {
#        Shiny.onInputChange("mydata", e.which);
#     });
#   ') 

