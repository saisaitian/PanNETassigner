\name{MLPvsother}
\alias{MLPvsother}
\title{
Find gene signature associated with combined MLP-1 and MLP-2 vs other PanNET subtypes.
}
\description{
Find gene signature associated with combined MLP-1 and MLP-2 vs other PanNET subtypes using SAM.
}
\usage{
MLPvsother(Indata, Indlabdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Indata}{
expression data
}
  \item{Indlabdata}{
Data with NMF labels 
}
}
\details{
Find gene signature associated with combined MLP-1 and MLP-2 vs other PanNET subtypes using SAM.
}
\value{
The output of this function include SAM genes saved in your working directory under PanNetResults/MLPvsother folder. 
}
\references{
Anguraj Sadanandam, et all (2015). A cross-species analysis in pancreatic neuroendocrine tumors reveals molecular subtypes with distinctive clinical, metastatic, developmental, and metabolic characteristics. Accepted.
}
\author{
Anguraj Sadanandam anguraj.sadanandam@icr.ac.uk 
}
\examples{
## Two datasets are required for this function 
data(MLPvsother.Indata)
data(MLPvsother.Indlabdata)
## Run the function
MLPvsother(MLPvsother.Indata,MLPvsother.Indlabdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
