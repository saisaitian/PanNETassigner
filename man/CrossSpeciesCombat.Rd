\name{CrossSpeciesCombat}
\alias{CrossSpeciesCombat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function for merging expression data from human and mouse PanNETs
}
\description{
Combat is used to remove batch effect and merge two expression datasets from two species, human and mouse.
}
\usage{
CrossSpeciesCombat(human, mouse, ordered)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{human}{
 is human.data dataset which contains the PanNET core clinical mRNA expression data.
}
  \item{mouse}{
 is mouse.data dataset which contains mouse RT2 mRNA expression data
}
  \item{ordered}{
 is the ordered samples across species.
}
}
\details{
This function selects the most variable genes (standard deviation greater than 0.8) in both datasets. Combat is used to remove batch effect and merge two gene expression datasets from different species.
}
\value{
The output of this function is the merged human and mouse expression data saved in your working directory under PanNetResults/CrossSpeciesCombat folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Accepted.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Three datasets are required for this function 
data(human.data)
data(mouse.data)
data(CrossSpeciesCombat.ordered)
## Run the function
CrossSpeciesCombat(human.data,mouse.data,CrossSpeciesCombat.ordered)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
