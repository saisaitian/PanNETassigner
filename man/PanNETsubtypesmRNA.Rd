\name{PanNETsubtypesmRNA}
\alias{PanNETsubtypesmRNA}
\title{
PanNET subtyping based on mRNA expression data.
}
\description{
This function selects genes with standard deviation greater than 0.8, employs NMF to discover subtypes and SAM
to identify genes associated with the subtypes.
}
\usage{
PanNETsubtypesmRNA(Indata)
}
\arguments{
  \item{Indata}{
 is PanNETsubtypesmRNA.Indata dataset which is mRNA expression data of tumors and islets. Genes are on the rows and the first column contain the gene names.
}
}
\details{
The dataset used here, PanNETsubtypesmRNA.Indata, is core clinical mRNA/gene expression PanNET data of tumors and islets. Genes with standard deviation greater than 0.8, employs NMF to discover subtypes and SAM
to identify genes associated with the subtypes.
}
\value{
The output of this function include NMF subtypes and SAM genes saved in your working directory under PanNetResults/mRNA folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Cancer Discovery.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Only one dataset required for this function 
data(PanNETsubtypesmRNA.Indata)
## Run the function
PanNETsubtypesmRNA(PanNETsubtypesmRNA.Indata)
}
