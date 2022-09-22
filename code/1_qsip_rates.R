#------------
#   Stone et al. 
#   2022
#   Nutrients strengthen density dependence of per-capita growth and mortality rates in the soil bacterial community
#   Oecologia
# -----------
# Calculation of bacterial growth and mortality rates using quantitative stable isotope probing


# load packages
library(data.table)
library(ggplot2)

# load data
asv <- fread('../data/feature_table.csv')
tax <- fread('../data/taxonomy_table.csv')
sam <- fwrite('../data/experimental_data.csv')

# combine abundances with sample-level data
dat <- merge(asv, sam[,c('SampleID', 'isotope', 'treatment',
                         'ecosystem', 'sampleID', 'Density.g.ml',
                         'avg_16S_g_soil', 'timepoint')], 
             by = 'SampleID',
             all.y = T)

# add in taxononmic information
dat <- merge(dat, tax, by = 'taxonID', all.x = T)

# convert all isotope treatments to either 18O or 16O, remove sample W1_PP_13, it shouldn't be here
dat <- dat[, isotope := fifelse(isotope == '18O', 'label', 'light')
           ][timepoint == 0, `:=` (isotope = NA, treatment = NA)
             ][sampleID != 'W1_PP_13']

# identify replicates for abundance matching in the growth/death rate calculations
reps <- unique(dat[, .(sampleID, isotope, treatment, ecosystem, timepoint)
                   ])[order(sampleID)
                      ][, rep := 1:.N, by = .(isotope, treatment, ecosystem, timepoint)]

dat <- merge(dat, reps[, .(sampleID, rep)], all.x = T, by = 'sampleID')

# plot density curves
ggplot(unique(dat[!is.na(Density.g.ml), c(1, 3:10)]),
       aes(Density.g.ml, avg_16S_g_soil, color = isotope)) +
  geom_line(aes(group = sampleID), alpha = 0.35) +
  facet_grid(treatment ~ ., scales = 'free_y')
# plot interactively
# plotly::ggplotly()

# Adjust for spin artifact in the C treatment
# adjustment calculated by taking the difference in median of "good" curves and median of "bad" curves
to_adjust <- paste0('W1_', c('PP_4', 'PJ_5', 'MC_6', 'GL_4', 'GL_5', 'PJ_4',
                             'PP_5', 'MC_5', 'MC_23', 'PJ_23', 'GL_24', 'PP_22',
                             'PJ_22', 'MC_22', 'PP_9'))
density_correct <- 0.0123115812168361760115
dat[sampleID %in% to_adjust, Density.g.ml := Density.g.ml - density_correct]
#
rm(density_correct, to_adjust)

# show corrected curves
ggplot(unique(dat[!is.na(Density.g.ml), c(1, 3:10)]),
       aes(Density.g.ml, avg_16S_g_soil, color = isotope)) +
  geom_line(aes(group = sampleID), alpha = 0.35) +
  facet_grid(treatment ~ ., scales = 'free_y')
# plot interactively
# plotly::ggplotly()


# Variables to control qSIP calculations
########################################
min_frac_threshold <- 3
min_rep_threshold <- 2
prop_measures_for_correction <- 0.1
mu <- 0.6
########################################

# how many sequences total?
dat[, sum(seq_abund)]
# 37,547,072

# remove non-bacterial lineages - need to also remove mitochondria and chloroplasts
dat <- dat[Kingdom == 'Bacteria'
           ][, tax_string := paste(Phylum, Class, Order, Family, Genus, Species, sep = ';')
             ][!grepl('mitochond|chloroplast', tax_string, ignore.case = T)
               ][, tax_string := NULL]

# how many taxa total and remaining sequence reads?
format(dat[, .(num_tax = uniqueN(taxonID), num_seq = sum(seq_abund))], big.mark = ',')
# 99,465 ASVs and 34,886,320 sequence reads

# normalize relative abundances since the qPCR data exclude chloroplasts/mitochondria and Archaea
dat[, rel_abund := seq_abund / sum(seq_abund), by = SampleID]

# how many 16S sequences (per g soil) total across all samples?
format(dat[, sum(avg_16S_g_soil * rel_abund, na.rm =T)], big.mark = ',')
# 33,000,171

# include list of optional grouping variable(s) or taxonomic variables you want to keep
vars_to_keep <- c('timepoint', 'rep', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

# include variable(s) th at you want to group your data by (such as site or nutrient addition)
# mostly important for taxa filtering
grouping_vars <- c('ecosystem', 'treatment')

# rename the tube-level and fraction-level columns to better distinguish them
setnames(dat, 
         c('SampleID', 'sampleID', 'taxonID'), 
         c('sample_fraction', 'sample_id', 'taxon_id'))

# save time 0 abundances for last steps
dat_t0 <- dat[timepoint == 0
              ][, abund_16s_t0 := rel_abund * avg_16S_g_soil
                ][, .(ecosystem, rep, taxon_id, abund_16s_t0)]


# filter taxa by fraction frequency---------------------------------------------

dat <- dat[, frac_freq := .N, by = .(taxon_id, sample_id)
           ][frac_freq >= min_frac_threshold
             ][, frac_freq := NULL]

dat[, uniqueN(taxon_id)]
# 13,299

# WADs--------------------------------------------------------------------------

# calculate WAD
dat <- dat[, abund_16s := rel_abund * avg_16S_g_soil, by = .(taxon_id, sample_fraction)
           ][, tot_abund := sum(abund_16s), by = .(taxon_id, sample_id) # group by replicate here, not sample-fraction
             ][, weight := abund_16s / tot_abund, by = .(taxon_id, sample_fraction)
               ][, .(wad = sum(weight * Density.g.ml, na.rm = T),
                     wvd = sum(weight * (Density.g.ml - sum(weight * Density.g.ml, na.rm = T))^2, na.rm = T),   # weighted variance of density
                     abund_16s = sum(abund_16s, na.rm = T)), 
                 by = c('taxon_id', 'sample_id', 'isotope', vars_to_keep, grouping_vars)
                 ][wad > 0]


# filter taxa by replicate frequency--------------------------------------------

# for 18O - must occur in minimum number of reps per group
dat <- dat[, rep_freq := uniqueN(sample_id), by = c('taxon_id', grouping_vars)
           ][rep_freq >= min_rep_threshold
             ][, rep_freq := NULL]

dat[, uniqueN(taxon_id)]
# 6,302


# Calculate population growth/turnover------------------------------------------

# Correct unlabeled tube WADs using Ben's method
# shift in WAD of all 36 unlabeled replicates from their global mean using taxa common to all replicates of all treatments
light_shift <- dat[isotope == 'light'
                   ][, num_reps := uniqueN(sample_id), by = taxon_id
                     ][num_reps == 36
                       ][, .(sample_id, taxon_id, wad)]
# 13 globally shared ASVs

# for each taxon, calculate the average or "true" WAD
# note: use median if only 1 or 2 tubes seem abnormal, use mean if there are distinct groups
# then, calculate how much each tube differs from the true WADs on average
light_shift <- light_shift[, true_wad := median(wad), by = taxon_id
                           ][, diff_from_true := wad - true_wad
                             ][, .(tube_shift = mean(diff_from_true)), by = sample_id]

# correct for density shifts in unlabeled WAD values
dat <- light_shift[dat, on = 'sample_id', 
                   ][is.na(tube_shift), tube_shift := 0
                     ][, wad := wad - tube_shift]

# convert to wide format
all_cols <- setdiff(names(dat), c('isotope', 'wad'))
wide_formula <- paste0(paste(all_cols, collapse = ' + '), ' ~ isotope')
#
dat <- dcast(dat, as.formula(wide_formula), value.var = 'wad', fill = NA)

# average light WADs by taxon
dat[, light := mean(light, na.rm = T), by = taxon_id]

dat[is.nan(light), light := NA
    ][is.nan(label), label := NA]

# remove outlier taxa with very heavy light WADs (plotting against labeled WADs shows they are uninformative)
light_outlier <- max(boxplot(dat$light, plot = F)$stats)
dat <- dat[light < light_outlier]

# # correct for labeled tubes using quantile regression, assuming taxa near the bottom of the wad_label ~ wad_light line are unenriched
# # we expect that for unenriched taxa, the slope should be close to 1 - but may differ based on the preparation and handling of the labeld tube
# # 3 ways (applying to each sample):
# #   1. Look for the 5th percentile line, forcing slope to be 1, and adjust this line to the 1:1 line
# #       - no way to easily do this using rq but equivalent to Ember's median correction based on WAD differences
# #   2. Look for 5th percentile line, forcing intercept to 0 and adjust this line to the 1:1 line
# #       - rq(label ~ light + 0, tau = 0.05)
# #       - dat[, label := label * (1 / slope), by = sample_id]
# #   3. Look for 5th percentile line, allowing slope and intercept to deviate
# #       - rg(label ~ light, tau = 0.05)
# #       - dat[, label := (label - intercept) * (1/ slope)]

# # get slope and intercept coefficients to perform correction
# get_rq_slope <- function(l, h) as.list(coef(quantreg::rq(h ~ l, tau = 0.05)))
# #
# label_shift <- dat[!is.na(label) & light < light_outlier
#                    ][, get_rq_slope(light, label), by = sample_id]
# setnames(label_shift, old = c('(Intercept)', 'l'), new = c('intercept', 'slope'))
# 
# dat <- merge(dat, label_shift, by = 'sample_id', all.x = T)
# dat[, label_corrected := (label - intercept) * (1 / slope)][, `:=` (intercept = NULL, slope = NULL)]

# Correct for labeled data using Ember's method (this might be similar to Ben's)
# for each tube based on the proportion of negative differences in WADs
label_shift <- dat[, .(wad_diff = label - light), by = sample_id
                   ][!is.na(wad_diff)
                     ][order(wad_diff)
                       ][, .(shift = median(wad_diff[1:floor(prop_measures_for_correction * .N)])), by = sample_id]

dat <- merge(dat, label_shift, by = 'sample_id', all.x = T)
dat[, label_corrected := label - shift][, shift := NULL]

# calculate molecular weights
dat[, gc_prop := (1 / 0.083506) * (light - 1.646057)
    ][, mw_light := (0.496 * gc_prop) + 307.691
      ][, `:=` (mw_label = (((label - light) / light) + 1) * mw_light,
                mw_label_corrected = (((label_corrected - light) / light) + 1) * mw_light,
                mw_max = mw_light + 12.07747 * mu)]


# calculate population rates
# merge time 0 abundances back in
dat <- merge(dat, dat_t0, by = c('ecosystem', 'taxon_id', 'rep'), all.x = T)

dat[, `:=` (abund_16s_light = abund_16s * ((mw_max - mw_label) / (mw_max - mw_light)),
            abund_16s_light_corrected = abund_16s * ((mw_max - mw_label_corrected) / (mw_max - mw_light)))
    ][, `:=` (growth = log(abund_16s / abund_16s_light) / timepoint,
              mortality = log(abund_16s_light / abund_16s_t0) / timepoint,
              growth_corrected = log(abund_16s / abund_16s_light_corrected) / timepoint,
              mortality_corrected = log(abund_16s_light_corrected / abund_16s_t0) / timepoint)]


# Clean output------------------------------------------------------------------

# remove NA rate values - this will also remove all the unlabeled samples
dat <- dat[!is.na(growth) & !is.na(mortality), 
           !c('label', 'light', 'gc_prop', 'mw_light', 'mw_label', 'mw_max', 'abund_16s_light',
              'label_corrected', 'mw_label_corrected', 'abund_16s_light_corrected')]

# how many taxa are left after filtering and that occurred in both unlabeled and labeled samples?
dat[, uniqueN(taxon_id)]
# 2,277 unique ASVs for which growth and mortality rates could be determined (and which passed frequency thresholds in the time 7 samples)

# # save output in your desired format!
# fwrite(dat, file = '../data/qsip_growth_mortality.csv')
