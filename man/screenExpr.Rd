\name{screenExpr}
\alias{screenExpr}
\title{
Screens data from the original expression files
}
\description{
Screens for genes with sample variance above user-specified threshold. Consolidates data for duplicate genes by averaging
across duplicates.
}
\usage{
screenExpr(Yexpr, sdCutoff)
}
\arguments{
  \item{Yexpr}{
Expression data.
}
  \item{sdCutoff}{
std dev threshold.
}
}
\details{
Screens for genes with sample variance above user-specified threshold. Consolidates data for duplicate genes by averaging
across duplicates.
}
\value{
Returns expression data with unique gene names and variance greater than psecified threshold
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Cancer Discovery.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Yexpr, sdCutoff) 
{
    data.all <- Yexpr
    print(paste("    Read expression data for ", dim(data.all)[1], 
        " genes.", sep = ""))
    data.mapped <- data.all[toupper(data.all[, 1]) != "UNMAPPED" & 
        !is.na(data.all[, 1]), ]
    print(paste("    Mapped data to ", dim(data.mapped)[1], " genes IDs.", 
        sep = ""))
    nLastCol <- dim(data.mapped)[2]
    print("Removing genes with low sample variance...")
    sampleVar <- apply(data.mapped[, 2:nLastCol], 1, var)
    idxHighVar <- sampleVar > sdCutoff^2
    data.highVar <- data.mapped[idxHighVar, ]
    sampleVar.highVar <- sampleVar[idxHighVar]
    print(paste("    Found ", dim(data.highVar)[1], " genes exceeding SD threshold of ", 
        sdCutoff, ".", sep = ""))
    nGenesHighVar <- dim(data.highVar)[1]
    genesUnique <- as.vector(unique(data.highVar[, 1]))
    nGenesUnique <- length(genesUnique)
    nSamples <- dim(data.highVar)[2] - 1
    data <- array(dim = c(nGenesUnique, nLastCol))
    data[, 1] <- genesUnique
    colnames(data) <- colnames(data.highVar)
    print("Removing duplicate genes (selecting for max standard deviation)...")
    for (gene in genesUnique) {
        idxGenes <- seq(along = 1:nGenesHighVar)[data.highVar[, 
            1] == gene]
        data.slice <- data.highVar[idxGenes, 2:nLastCol]
        if (length(idxGenes) > 1) {
            idxMaxVar <- which.max(sampleVar.highVar[idxGenes])
            data.slice <- data.slice[idxMaxVar, ]
        }
        data[data[, 1] == gene, 2:nLastCol] <- as.matrix(data.slice)
    }
    print(paste("    ", nGenesUnique, " unique genes IDs remain.", 
        sep = ""))
    print("Screened data is the output:... ")
    print("")
    rm(data.all, data.mapped, data.highVar)
    sdata <- data
    sdata
  }
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
