\name{PanNETsubtypesmicRNAVal}
\alias{PanNETsubtypesmicRNAVal}
\title{
Validation of micRNA expression PanNET subtypes using miR expression validation dataset.
}
\description{
This function selects SAM gene signature for PanNET miR expression subtypes in the validation dataset. NMF is employed to
indentify PanNET subtypes in this dataset.
}
\usage{
PanNETsubtypesmicRNAVal(Indata, IndSAMdata)
}
\arguments{
  \item{Indata}{
 PanNETsubtypesmicRNAVal.Indata is the validation dataset. microRNAs are on the rows and the first column contain the microRNA names.
 }
  \item{IndSAMdata}{
PanNETsubtypesmicRNAVal.IndSAMdata contain SAM micRNAs signature for finding the PanNET subtypes.
}
}
\details{
The SAM gene signature in PanNETsubtypesmicRNAVal.IndSAMdata is selected from the validation dataset PanNETsubtypesmicRNAVal.Indata and NMF is employed to discover PanNET miR expression subtypes in this dataset. 
}
\value{
The output of this function include NMF micRNA expression subtypes saved in your working directory under PanNetResults/micRNAVal folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics.  Cancer Discovery.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Two datasets are required for this function 
data(PanNETsubtypesmicRNAVal.Indata)
data(PanNETsubtypesmicRNAVal.IndSAMdata)
## Run the function
PanNETsubtypesmicRNAVal(PanNETsubtypesmicRNAVal.Indata,PanNETsubtypesmicRNAVal.IndSAMdata)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
