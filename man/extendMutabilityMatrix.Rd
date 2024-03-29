% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TargetingModels.R
\name{extendMutabilityMatrix}
\alias{extendMutabilityMatrix}
\title{Extends a mutability model to include Ns.}
\usage{
extendMutabilityMatrix(mutabilityModel)
}
\arguments{
\item{mutabilityModel}{vector of 5-mer mutability rates built by 
\link{createMutabilityMatrix}.}
}
\value{
A \code{MutabilityModel} containing a 3125 vector of normalized 
          mutability rates for each 5-mer motif with names defining the 5-mer 
          nucleotide sequence. Note that "normalized" means that the mutability 
          rates for the 1024 5-mers that contain no "N" at any position sums up 
          to 1 (as opposed to the entire vector summing up to 1). 
          
          If the input \code{mutabilityModel} is of class \code{MutabilityModel}, 
          then the output \code{MutabilityModel} will carry over the input 
          \code{numMutS} and \code{numMutR} slots.
}
\description{
\code{extendMutabilityMatrix} extends a 5-mer nucleotide mutability model 
with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
}
\examples{
\donttest{
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
set.seed(112)
db <- dplyr::slice_sample(db, n=75)

# Create model using only silent mutations and ignore multiple mutations
sub_model <- createSubstitutionMatrix(db, model="s", sequenceColumn="sequence_alignment",
                                      germlineColumn="germline_alignment_d_mask",
                                      vCallColumn="v_call")
mut_model <- createMutabilityMatrix(db, sub_model, model="s", 
                                    sequenceColumn="sequence_alignment",
                                    germlineColumn="germline_alignment_d_mask",
                                    vCallColumn="v_call")
ext_model <- extendMutabilityMatrix(mut_model)
}

}
\seealso{
\link{createMutabilityMatrix}, \link{extendSubstitutionMatrix}, 
          \link{MutabilityModel}
}
