% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Shmulate.R
\name{shmulateTree}
\alias{shmulateTree}
\title{Simulate mutations in a lineage tree}
\usage{
shmulateTree(
  sequence,
  graph,
  targetingModel = HH_S5F,
  field = NULL,
  exclude = NULL,
  junctionWeight = NULL,
  start = 1,
  end = nchar(sequence)
)
}
\arguments{
\item{sequence}{string defining the MRCA sequence to seed mutations from.}

\item{graph}{\code{igraph} object defining the seed lineage tree, with 
vertex annotations, whose edges are to be recreated.}

\item{targetingModel}{5-mer \link{TargetingModel} object to be used for computing 
probabilities of mutations at each position. Defaults to
\link{HH_S5F}.}

\item{field}{annotation to use for both unweighted path length exclusion 
and consideration as the MRCA node. If \code{NULL} do not 
exclude any nodes.}

\item{exclude}{vector of annotation values in \code{field} to exclude from 
potential MRCA set. If \code{NULL} do not exclude any nodes.
Has no effect if \code{field=NULL}.}

\item{junctionWeight}{fraction of the nucleotide sequence that is within the 
junction region. When specified this adds a proportional 
number of mutations to the immediate offspring nodes of the 
MRCA. Requires a value between 0 and 1. If \code{NULL} then 
edge weights are unmodified from the input \code{graph}.}

\item{start}{Initial position in \code{sequence} where mutations can 
be introduced. Default: 1}

\item{end}{Last position in \code{sequence} where mutations can 
be introduced. Default: last position (sequence length).}
}
\value{
A \code{data.frame} of simulated sequences with columns:
          \itemize{
            \item \code{name}:      name of the corresponding node in the input 
                                    \code{graph}.  
            \item \code{sequence}:  mutated sequence.
            \item \code{distance}:  Hamming distance of the mutated sequence from 
                                    the seed \code{sequence}.
          }
}
\description{
\code{shmulateTree} returns a set of simulated sequences generated from an input 
sequence and a lineage tree. The input sequence is used to replace the most recent 
common ancestor (MRCA) node of the \code{igraph} object defining the lineage tree. 
Sequences are then simulated with mutations corresponding to edge weights in the tree. 
Sequences will not be generated for groups of nodes that are specified to be excluded.
}
\examples{
# Load example lineage and define example MRCA
data(ExampleTrees, package="alakazam")
graph <- ExampleTrees[[17]]
sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"

# Simulate using the default human 5-mer targeting model
shmulateTree(sequence, graph)

# Simulate using the mouse 5-mer targeting model
# Exclude nodes without a sample identifier
# Add 20\% mutation rate to the immediate offsprings of the MRCA
shmulateTree(sequence, graph, targetingModel=MK_RS5NF,
             field="sample_id", exclude=NA, junctionWeight=0.2)
 
}
\seealso{
See \link{shmulateSeq} for imposing mutations on a single sequence. 
          See \link{HH_S5F} and \link{MK_RS5NF} for predefined 
          \link{TargetingModel} objects.
}
