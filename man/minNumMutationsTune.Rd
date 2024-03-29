% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TargetingModels.R
\name{minNumMutationsTune}
\alias{minNumMutationsTune}
\title{Parameter tuning for minNumMutations}
\usage{
minNumMutationsTune(subCount, minNumMutationsRange)
}
\arguments{
\item{subCount}{\code{data.frame} returned by \link{createSubstitutionMatrix}
with \code{numMutationsOnly=TRUE}.}

\item{minNumMutationsRange}{a number or a vector indicating the value or range of values
of \code{minNumMutations} to try.}
}
\value{
A 3xn \code{matrix}, where n is the number of trial values of \code{minNumMutations}
             supplied in \code{minNumMutationsRange}. Each column corresponds to a value
             in \code{minNumMutationsRange}. The rows correspond to the number of 5-mers
             for which substitution rates would be computed directly using the 5-mer itself 
             (\code{"5mer"}), using its inner 3-mer (\code{"3mer"}), and using the central 
             1-mer (\code{"1mer"}), respectively.
}
\description{
\code{minNumMutationsTune} helps with picking a threshold value for \code{minNumMutations}
in \link{createSubstitutionMatrix} by tabulating the number of 5-mers for which 
substitution rates would be computed directly or inferred at various threshold values.
}
\details{
At a given threshold value of \code{minNumMutations}, for a given 5-mer,
             if the total number of mutations is greater than the threshold and there
             are mutations to every other base, substitution rates are computed directly
             for the 5-mer using its mutations. Otherwise, mutations from 5-mers with 
             the same inner 3-mer as the 5-mer of interest are aggregated. If the number 
             of such mutations is greater than the threshold and there are mutations to 
             every other base, these mutations are used for inferring the substitution 
             rates for the 5-mer of interest; if not, mutations from all 5-mers with the 
             same center nucleotide are aggregated and used for inferring the substitution
             rates for the 5-mer of interest (i.e. the 1-mer model).
}
\examples{
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")

# Count the number of mutations per 5-mer
subCount <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
                                     germlineColumn="germline_alignment_d_mask",
                                     vCallColumn="v_call",
                                     model="s", multipleMutation="independent",
                                     returnModel="5mer", numMutationsOnly=TRUE)

# Tune minNumMutations
minNumMutationsTune(subCount, seq(from=10, to=80, by=10))
                                      
}
\references{
\enumerate{
  \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
           on synonymous mutations from high-throughput immunoglobulin sequencing data. 
           Front Immunol. 2013 4(November):358.
 }
}
\seealso{
See argument \code{numMutationsOnly} in \link{createSubstitutionMatrix} 
             for generating the required input \code{data.frame} \code{subCount}. 
             See argument \code{minNumMutations} in \link{createSubstitutionMatrix}
             for what it does.
}
