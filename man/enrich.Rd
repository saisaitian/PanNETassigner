\name{enrich}
\alias{enrich}
\title{ Performs hypergeometric test of association between miR subtypes and mRNA subtypes.
}
\description{
This function use Hypergemotric test to determine association between miR and mRNA subtypes.
}
\usage{
enrich(labels1, labels2, outfile)
}
\arguments{
  \item{labels1}{
 EnrichmentAnalysis.labels1 contain the miR subtype labels and the first column contain the sample names.
}
  \item{labels2}{
 EnrichmentAnalysis.labels2 contain the mRNA subtype labels and the first column contain the sample names.
}
  \item{outfile}{
Provide the name of the output file, i.e. "miR_subtypes-mRNA_subtypes.txt".
}
}
\details{
Hypergemotric test is used to determine association between miR in and mRNA subtypes, in EnrichmentAnalysis.labels1 and EnrichmentAnalysis.labels2 respectively.
}
\value{
The output of this function is a heatmap showing the relationships between the miR and mRNA subtypes. Red colour denotes strong association and other colors indicate that there is no assocation. This is saved in your working directory under PanNetResults/Enrichment folder. 

}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Accepted.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Two datasets are required for this function and please provide the output filename. 
data(EnrichmentAnalysis.labels1)
data(EnrichmentAnalysis.labels2)
outfile<-"miR_subtypes-mRNA_subtypes.txt"
## Run the function
enrich(EnrichmentAnalysis.labels1, EnrichmentAnalysis.labels2, outfile)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
