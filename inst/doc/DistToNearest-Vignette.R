## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Subset example data to one sample
library(shazam)
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, SAMPLE == "-1h")

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Use nucleotide Hamming distance and normalize by junction length
dist_ham <- distToNearest(db, model="ham", first=FALSE, normalize="length", 
                          nproc=1)

# Use genotyped V assignments and 5-mer model
dist_hs5f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", model="hs5f", 
                           first=FALSE, normalize="none", nproc=1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate Hamming distance histogram
library(ggplot2)
p1 <- ggplot(subset(dist_ham, !is.na(DIST_NEAREST)),
             aes(x=DIST_NEAREST)) + 
    theme_bw() + xlab("Hamming distance") + ylab("Count") +
    scale_x_continuous(breaks=seq(0, 1, 0.1)) +
    geom_histogram(fill="steelblue", color="white", binwidth=0.02) +
    geom_vline(xintercept=0.15, color="firebrick", linetype=3)
plot(p1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate hs5f distance histogram
p2 <- ggplot(subset(dist_hs5f, !is.na(DIST_NEAREST)),
             aes(x=DIST_NEAREST)) + 
    theme_bw() + xlab("HS5F distance") + ylab("Count") +
    scale_x_continuous(breaks=seq(0, 50, 5)) +
    geom_histogram(fill="steelblue", color="white", binwidth=1) +
    geom_vline(xintercept=7, color="firebrick", linetype=3)
plot(p2)

## ----fields, eval=TRUE, warning=FALSE------------------------------------
dist_fields <- distToNearest(ExampleDb, model="ham", first=FALSE, 
                             normalize="length", fields="SAMPLE", 
                             nproc=1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate grouped histograms
p3 <- ggplot(subset(dist_fields, !is.na(DIST_NEAREST)), 
             aes(x=DIST_NEAREST)) + 
    theme_bw() + xlab("Grouped Hamming distance") + ylab("Count") +
    geom_histogram(fill="steelblue", color="white", binwidth=0.02) +
    geom_vline(xintercept=0.15, color="firebrick", linetype=3) +
    facet_grid(SAMPLE ~ ., scales="free_y")
plot(p3)

## ----cross, eval=TRUE, warning=FALSE-------------------------------------
dist_cross <- distToNearest(ExampleDb, model="ham", first=FALSE, 
                            normalize="length", cross="SAMPLE", nproc=1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate cross sample histograms
p4 <- ggplot(subset(dist_cross, !is.na(CROSS_DIST_NEAREST)), 
             aes(x=CROSS_DIST_NEAREST)) + 
    theme_bw() + xlab("Cross-sample Hamming distance") + ylab("Count") +
    geom_histogram(fill="steelblue", color="white", binwidth=0.02) +
    geom_vline(xintercept=0.15, color="firebrick", linetype=3) +
    facet_grid(SAMPLE ~ ., scales="free_y")
plot(p4)

