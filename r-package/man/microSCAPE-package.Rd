\name{microSCAPE-package}
\alias{microSCAPE-package}
\alias{microSCAPE}
\docType{package}
\title{
  Downstream Analysis Tools for Amplicon Sequencing Data
}
\description{
  Provides downstream analysis tools for amplicon sequencing data produced
  by DADA2. Includes sequence table quality filtering, MIMARKS metadata
  loading, taxonomic renormalization, phylogenetic tree construction,
  Bray-Curtis ordination with t-SNE, and SparCC-style correlation network
  analysis. Designed as a companion to the dada2 Bioconductor package.
}
\details{
  The main functions in this package are:

  \itemize{
    \item \code{\link{filterSeqtab}} -- Cascade quality-control filtering
    \item \code{\link{loadMetadata}} -- Load MIMARKS-compliant metadata
    \item \code{\link{renormalize}} -- Separate and normalize by taxonomic group
    \item \code{\link{buildPhylogeny}} -- Multiple alignment and NJ tree
    \item \code{\link{ordinateSamples}} -- Bray-Curtis ordination (PCA / t-SNE)
    \item \code{\link{sparccNetwork}} -- CLR-based correlation network
  }

  See \code{vignette("microSCAPE")} for a usage overview.
}
\author{
  R. Eric Collins \email{recollins@alaska.edu}

  Maintainer: R. Eric Collins \email{recollins@alaska.edu}
}
\keyword{package}
