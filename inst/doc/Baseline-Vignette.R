## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Load example data
library(shazam)
data(ExampleDb, package="alakazam")

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Collapse clonal groups into single sequences
clones <- collapseClones(ExampleDb, regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Count observed mutations and append OBSERVED columns to the output
observed <- observedMutations(clones, 
                              sequenceColumn="CLONAL_SEQUENCE",
                              regionDefinition=IMGT_V_NO_CDR3, nproc=1)
# Count observed mutations and append EXPECTED columns to the output
expected <- expectedMutations(observed, 
                              sequenceColumn="CLONAL_SEQUENCE",
                              targetingModel=HS5FModel,
                              regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Calculate selection scores using the output from expectedMutations
baseline <- calcBaseline(expected, testStatistic="focused", 
                         regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Calculate selection scores from scratch on subset
baseline <- calcBaseline(ExampleDb, testStatistic="focused", 
                         regionDefinition=IMGT_V_NO_CDR3, nproc=1)

# Subset the original data to switched isotypes
db_sub <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG"))
# Calculate selection scores from scratch on subset
baseline_sub <- calcBaseline(db_sub, testStatistic="focused", 
                             regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=FALSE, warning=FALSE, results="hide"--------------------------
#  # Calculate selection on charge class with the mouse 5-mer model
#  baseline <- calcBaseline(ExampleDb, testStatistic="focused",
#                           regionDefinition=IMGT_V_NO_CDR3,
#                           targetingModel=MRS5NFModel,
#                           targetingModel=CHARGE_MUTATIONS,
#                           nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Combine selection scores by time-point
grouped_1 <- groupBaseline(baseline, groupBy=c("SAMPLE"))

# Combine selection scores by time-point and isotype
grouped_2 <- groupBaseline(baseline_sub, groupBy=c("SAMPLE", "ISOTYPE"))

## ---- eval=FALSE, warning=FALSE, results="hide"--------------------------
#  # First group by subject and status
#  subject_grouped <- groupBaseline(baseline, groupBy=c("STATUS", "SUBJECT"))
#  
#  # Then group the output by status
#  status_grouped <- groupBaseline(subject_grouped, groupBy="STATUS")

## ---- eval=TRUE----------------------------------------------------------
testBaseline(grouped_1, groupBy="SAMPLE")

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Set sample and isotype colors
sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
isotype_colors <- c("IgM"="darkorchid", "IgD"="firebrick", 
                    "IgG"="seagreen", "IgA"="steelblue")

# Plot mean and confidence interval by time-point
plotBaselineSummary(grouped_1, "SAMPLE")

# Plot selection scores by time-point and isotype for only CDR
plotBaselineSummary(grouped_2, "SAMPLE", "ISOTYPE", groupColors=isotype_colors,
                    subsetRegions="CDR")

# Group by CDR/FWR and facet by isotype
plotBaselineSummary(grouped_2, "SAMPLE", "ISOTYPE", facetBy="group")

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Plot selection PDFs for a subset of the data
plotBaselineDensity(grouped_2, "ISOTYPE", groupColumn="SAMPLE", colorElement="group", 
                    colorValues=sample_colors, sigmaLimits=c(-1, 1))

