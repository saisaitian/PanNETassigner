\name{PanNETsubtypesmicRNA}
\alias{PanNETsubtypesmicRNA}
\title{
PanNET subtyping based on core clinical miR expression data.
}
\description{
This function selects microRNAs with standard deviation greater than 0.5, employs NMF to discover subtypes and SAM
to identify microRNAs associated with the subtypes.
}
\usage{
PanNETsubtypesmicRNA(Indata)
}
\arguments{
  \item{Indata}{
 is the PanNETsubtypesmicRNA.Indata dataset which is micRNA expression data of tumors. micRNAs are on the rows and the first column contain the micRNA names.
}
}
\details{
The dataset used here, PanNETsubtypesmicRNA.Indata, is microRNA expression data of tumors. micRNAs with standard deviation greater than 0.5 are selected, NMF is employed to discover subtypes and SAM to identify microRNAs associated with the subtypes.
}
\value{
The output of this function include NMF subtypes and associated SAM microRNAs saved in your working directory under PanNetResults/micRNA folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Cancer Discovery.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Only one dataset required for this function 
data(PanNETsubtypesmicRNA.Indata)
## Run the function
PanNETsubtypesmicRNA(PanNETsubtypesmicRNA.Indata)
}