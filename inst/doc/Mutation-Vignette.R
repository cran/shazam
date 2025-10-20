## ----eval=TRUE, warning=FALSE, message=FALSE----------------------------------
# Import required packages
library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)

# Load and subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call %in% c("IGHA", "IGHG") & sample_id == "+7d")

## ----eval=TRUE----------------------------------------------------------------
# Calculate R and S mutation counts
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
# Show new mutation count columns
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## ----eval=TRUE----------------------------------------------------------------
# Calculate combined R and S mutation frequencies
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq")) %>%
  head(n=4)

## ----eval=TRUE, warning=FALSE-------------------------------------------------
g1 <- ggplot(db_obs, aes(x=c_call, y=mu_freq, fill=c_call)) +
        geom_boxplot() + 
        labs(title="Total mutations", 
             x="Isotype", y="Mutation frequency") +
        scale_fill_manual(name="Isotype", values=IG_COLORS, limits=force) +
        theme_bw() 
plot(g1)

## ----eval=TRUE----------------------------------------------------------------
# Calculate R and S mutation counts for individual CDRs and FWRs
db_obs_v <- observedMutations(db, sequenceColumn="sequence_alignment",
                              germlineColumn="germline_alignment_d_mask",
                              regionDefinition=IMGT_V_BY_REGIONS,
                              frequency=FALSE, 
                              nproc=1)
# Show new FWR mutation columns
db_obs_v %>% 
  select(sequence_id, starts_with("mu_count_fwr")) %>%
  head(n=4)

# Calculate aggregate CDR and FWR V-segment R and S mutation frequencies
db_obs_v <- observedMutations(db_obs_v, sequenceColumn="sequence_alignment",
                              germlineColumn="germline_alignment_d_mask",
                              regionDefinition=IMGT_V,
                              frequency=TRUE, 
                              nproc=1)
# Show new CDR and FWR mutation frequency columns
db_obs_v %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## ----eval=TRUE, warning=FALSE-------------------------------------------------
g2 <- ggplot(db_obs_v, aes(x=c_call, y=mu_freq_cdr_s, fill=c_call)) +
        geom_boxplot() + 
        labs(title="CDR silent mutations", 
             x="Isotype", y="Mutation frequency") +
        scale_fill_manual(name="Isotype", values=IG_COLORS, limits=force) +
        theme_bw()
g3 <- ggplot(db_obs_v, aes(x=c_call, y=mu_freq_cdr_r, fill=c_call)) +
        geom_boxplot() + 
        labs(title="CDR replacement mutations",
             x="Isotype", y="Mutation frequency") +
        scale_fill_manual(name="Isotype", values=IG_COLORS, limits=force) +
        theme_bw()

alakazam::gridPlot(g2, g3, ncol=2)

## ----eval=TRUE----------------------------------------------------------------
# Calculate charge mutation frequency for the full sequence
db_obs_ch <- observedMutations(db, sequenceColumn="sequence_alignment",
                               germlineColumn="germline_alignment_d_mask",
                               regionDefinition=NULL,
                               mutationDefinition=CHARGE_MUTATIONS,
                               frequency=TRUE, 
                               nproc=1)
# Show new charge mutation frequency columns
db_obs_ch %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## ----eval=TRUE, warning=FALSE-------------------------------------------------
g4 <- ggplot(db_obs_ch, aes(x=c_call, y=mu_freq_seq_r, fill=c_call)) +
        geom_boxplot() + 
        labs(title="Charge replacement mutations", 
             x="Isotype", y="Mutation frequency") +
        scale_fill_manual(name="Isotype", values=IG_COLORS, limits=force) + 
        theme_bw()

plot(g4)

