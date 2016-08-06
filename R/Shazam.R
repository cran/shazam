# shazam package documentation and import directives

#' The shazam package
#'
#' Dramatic improvements in high-throughput sequencing technologies now enable 
#' large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
#' antigen-receptor proteins located on the surface of T and B lymphocytes. The \code{shazam}
#' package provides tools for advanced analysis of somatic hypermutation (SHM) in
#' immunoglobulin (Ig) sequences. The key functions in \code{shazam}, broken down topic, are 
#' described below.
#' 
#' @section  Mutational profiling:
#' \code{shazam} provides tools to quantify the extent and nature of SHM within
#' full length V(D)J sequences as well as sub-regions (eg, FWR and CDR).
#' Quantification of expected mutational loaded, under specific SHM targeting 
#' models, can also be performed along with model driven simulations of SHM.
#' 
#' \itemize{
#'   \item  \link{collapseClones}:           Build clonal consensus sequence.
#'   \item  \link{observedMutations}:        Compute observed mutation counts and frequencies.
#'   \item  \link{expectedMutations}:        Compute expected mutation frequencies.
#'   \item  \link{shmulateSeq}:              Simulate mutations in a single sequence.
#'   \item  \link{shmulateTree}:             Simulate mutations over a lineage tree.
#' }
#' 
#' @section  SHM targeting models:
#' Computational models and analyses of SHM have separated the process 
#' into two independent components: 
#' \enumerate{
#'   \item  A mutability model that defines where mutations occur.
#'   \item  A nucleotide substitution model that defines the resulting mutation.
#' }
#' Collectively these are what form the targeting model of SHM. \code{shazam} 
#' provides empirically derived targeting models for both humans and mice,
#' along with tools to build these mutability and substitution models from data.
#' 
#' \itemize{
#'   \item  \link{createTargetingModel}:     Build a 5-mer targeting model.
#'   \item  \link{plotMutability}:           Plot 5-mer mutability rates.
#'   \item  \link{HS5FModel}:                Human 5-mer SHM targeting model.
#'   \item  \link{MRS5NFModel}:              Mouse 5-mer SHM targeting model.
#' }
#'
#' @section  Quantification of selection pressure:
#' Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
#' novel method for quantifying antigen-driven selection in high-throughput
#' Ig sequence data. Targeting models created using \code{shazam} can be used 
#' to estimate the null distribution of expected mutation frequencies used
#' by BASELINe, providing measures of selection pressure informed by known 
#' AID targeting biases.
#'              
#' \itemize{
#'   \item  \link{calcBaseline}:             Calculate the BASELINe probability
#'                                           density functions (PDFs).
#'   \item  \link{groupBaseline}:            Combine PDFs from sequences grouped
#'                                           by biological or experimental relevance.
#'   \item  \link{summarizeBaseline}:        Compute summary statistics from BASELINe PDFs.
#'   \item  \link{testBaseline}:             Perform significance testing for the difference
#'                                           between BASELINe PDFs.
#'   \item  \link{plotBaselineDensity}:      Plot the probability density functions
#'                                           resulting from selection analysis.
#'   \item  \link{plotBaselineSummary}:      Plot summary stastistics resulting from 
#'                                           selection analysis.
#' }
#'
#' @section  Mutational distance calculation:
#' \code{shazam} provides tools to compute evolutionary distances between 
#' sequences or groups of sequences, which can leverage SHM targeting 
#' models. This information is particularly useful in understanding and 
#' defining clonal relationships.
#'              
#' \itemize{
#'   \item  \link{distToNearest}:            Tune clonal assignment thresholds by calculating 
#'                                           distances to nearest-neighbors.
#'   \item  \link{calcTargetingDistance}:    Construct a nucleotide distance matrix from a 
#'                                           5-mer targeting model.
#' }
#'
#' @name     shazam
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'   \item  Cui A, et al. A model of somatic hypermutation targeting in mice based on 
#'            high-throughput immunoglobulin sequencing data. Under review.
#'  }
#' 
#' @import   ggplot2
#' @import   graphics
#' @import   methods
#' @import   utils
#' @importFrom  alakazam    getAllele getGene getFamily getSegment
#'                          getAAMatrix getDNAMatrix
#'                          pairwiseDist pairwiseEqual seqDist seqEqual
#'                          isValidAASeq translateStrings gridPlot
#'                          getMRCA getPathLengths tableEdges
#' @importFrom  ape         mst
#' @importFrom  data.table  data.table setkey setkeyv
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       do n desc funs %>%
#'                          as_data_frame data_frame data_frame_
#'                          bind_cols bind_rows combine
#'                          filter filter_ select select_ arrange arrange_
#'                          group_by group_by_ ungroup
#'                          mutate mutate_ summarize summarize_
#'                          mutate_each mutate_each_ summarize_each summarize_each_
#'                          rename rename_ transmute transmute_
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  igraph      V E as_adjacency_matrix graph_from_data_frame
#'                          vertex_attr set_vertex_attr
#' @importFrom  lazyeval    interp
#' @importFrom  scales      log2_trans log10_trans trans_breaks trans_format
#'                          math_format percent scientific
#' @importFrom  tidyr       gather gather_ spread spread_
#' @importFrom  iterators   icount
#' @importFrom  SDMTools    wt.sd
#' @importFrom  seqinr      c2s s2c words translate
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          approx convolve weighted.mean p.adjust
#'                          dbeta pbeta qbeta rbeta
#' @importFrom  stringi     stri_dup stri_flatten stri_join stri_length
#'                          stri_count_boundaries stri_count_regex 
#'                          stri_extract_all_regex stri_extract_first_regex  
#'                          stri_replace_all_regex stri_replace_first_regex
NULL


#### Sysdata ####

# Ordered nucleotide character set
# NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

# IMGT V segment length
# VLENGTH <- 312

# 5x312 logical matrix of CDR positions
# CDR_Nuc_Mat

# 5x312 logical matrix of FWR positions
# FWR_Nuc_Mat

# 12x216 matrix of replacement and silent mutation permutations
# CODON_TABLE

# 1x24 vector of amino acid charge classes
# AMINO_ACIDS_CHARGE

# 1x24 vector of amino acid hydropathy classes
# AMINO_ACIDS_HYDROPATHY

# 1x24 vector of amino acid polarity classes
# AMINO_ACIDS_POLARITY

# TODO: What is this?
# CONST_I

# TODO: And what is this?
# BAYESIAN_FITTED

# Add built-in variables to global variables environment
utils::globalVariables(c("M1NDistance", "HS1FDistance", 
                         "U5NModel", "HS5FModel", "MRS5NFModel"), package="shazam")