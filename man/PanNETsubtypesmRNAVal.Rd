\name{PanNETsubtypesmRNAVal}
\alias{PanNETsubtypesmRNAVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validation of gene expression PanNET subtypes using gene expression validation dataset.
}
\description{
This function selects PanNETassinger-mRNA signature from the gene expression validation dataset. NMF is employed to indentify PanNET subtypes using this dataset.
}
\usage{
PanNETsubtypesmRNAVal(Indata, IndSAMdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Indata}{
 PanNETsubtypesmRNAVal.Indata is the  gene expression validation dataset. Genes are on the rows and the first column contain the gene names.
}
  \item{IndSAMdata}{
PanNETsubtypesmRNAVal.IndSAMdata contain 221 PanNETassigner-mRNA signature for finding the PanNET subtypes.
}
}
\details{
The 221 PanNETassigner-mRNA signature in PanNETsubtypesmRNAVal.IndSAMdata is selected from PanNETsubtypesmRNAVal.Indata, gene expression validation dataset, and NMF is used to discover PanNET mRNA expression subtypes in this dataset. 
}
\value{
The output of this function include NMF mRNA expression subtypes saved in your working directory under PanNetResults/mRNAVal folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Cancer Discovery.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
### Two datasets are required for this function 
data(PanNETsubtypesmRNAVal.Indata)
data(PanNETsubtypesmRNAVal.IndSAMdata)
### Run the function
PanNETsubtypesmRNAVal(PanNETsubtypesmRNAVal.Indata,PanNETsubtypesmRNAVal.IndSAMdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
