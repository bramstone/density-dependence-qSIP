#------------
#   Stone et al. 
#   2022
#   Nutrients strengthen density dependence of per-capita growth and mortality rates in the soil bacterial community
#   Oecologia
# -----------
# relationship between density dependence and alpha diversity

library(data.table)
library(lme4)
library(lmerTest)
library(ggplot2)

dat <- fread('../data/qsip_growth_mortality.csv')
gcn <- fread('../data/qsip_taxa_16s_gcn_nstd.csv')
#
asv <- fread('../data/feature_table.csv')
tax <- fread('../data/taxonomy_table.csv')
sam <- fread('../data/experimental_data.csv')


dat <- merge(dat, gcn[, .(taxonID, num_16s, nstd)], all.x = T, by.x = 'taxon_id', by.y = 'taxonID')

dat[nstd > 0.15, num_16s := 2][, nstd := NULL]

dat[, `:=` (treatment = factor(treatment, levels = c('control', 'C', 'C_N')),
            growth = growth_corrected,
            mortality = mortality * -1,
            wvd_sqrt = sqrt(wvd),
            pop_log10 = log10(abund_16s / num_16s),
            pop_log10_t0 = log10(abund_16s_t0 / num_16s))]

# how many taxa initially?
dat[, uniqueN(taxon_id)]
dat[, uniqueN(taxon_id), by = ecosystem]
# 2300 with ~ 650-950 per ecosystem

# limit analyses to ASVs that occur in every treatment at least three times
dat <- dat[, n_trt := uniqueN(treatment), by = .(taxon_id, ecosystem)
           ][n_trt == 3
             ][, n_trt := NULL
               ][, n_rep := .N, b = .(taxon_id, treatment)
                 ][n_rep > 2
                   ][, n_rep := NULL]

# how many taxa now?
dat[, uniqueN(taxon_id)]
# 524


# combine abundances with sample-level data
asv <- merge(asv, sam[,c('SampleID', 'isotope', 'treatment',
                         'ecosystem', 'sampleID', 'Density.g.ml',
                         'avg_16S_g_soil', 'timepoint')], 
             by = 'SampleID',
             all.y = T)

# remove non-bacterial sequences from full community feature table
bact <- tax[Kingdom == 'Bacteria'
            ][, tax_string := paste(Phylum, Class, Order, Family, Genus, Species, sep = ';')
              ][!grepl('mitochond|chloroplast', tax_string, ignore.case = T)
                ][, taxonID]

asv <- asv[taxonID %in% bact]

# convert all isotope treatments to either 18O or 16O, remove sample W1_PP_13, it shouldn't be here
asv <- asv[, isotope := fifelse(isotope == '18O', 'label', 'light')
           ][timepoint == 0, `:=` (isotope = NA, treatment = NA)
             ][sampleID != 'W1_PP_13']

# identify replicates for abundance matching in the growth/death rate calculations
reps <- unique(asv[, .(sampleID, isotope, treatment, ecosystem, timepoint)
                   ])[order(sampleID)
                      ][, rep := 1:.N, by = .(isotope, treatment, ecosystem, timepoint)]

asv <- merge(asv, reps[, .(sampleID, rep)], all.x = T, by = 'sampleID')

# normalize relative abundances since the qPCR data exclude chloroplasts/mitochondria and Archaea
asv[, rel_abund := seq_abund / sum(seq_abund), by = SampleID]

asv <- asv[, .(rel_abund = sum(rel_abund * avg_16S_g_soil) / sum(avg_16S_g_soil)),
  by = .(taxonID, isotope, treatment, ecosystem, sampleID, timepoint)]


# per-sample density dependence calculation
set.seed(32122)

# density dependence
sam_g_lme <- lmer(growth ~ 1 +
                    (pop_log10_t0|sample_id),
                  data = dat, REML = T)
# doesn't converge - but that's OK because we're not using statistical output

sam_m_lme <- lmer(mortality ~ 1 +
                    (pop_log10_t0|sample_id),
                  data = dat, REML = T)

# generate net density dependence effect
dd <- data.table(sampleID = rownames(coef(sam_g_lme)$sample_id),
                 growth_dd = coef(sam_g_lme)$sample_id$pop_log10_t0,
                 mortality_dd = coef(sam_m_lme)$sample_id$pop_log10_t0
                 )[, net_dd := growth_dd - mortality_dd]


dat[, in_qsip := T]

asv <- merge(asv, dat[, .(taxon_id, sample_id, in_qsip)], all.x = T,
                  by.x = c('taxonID', 'sampleID'), by.y = c('taxon_id', 'sample_id'))

asv[is.na(in_qsip), in_qsip := F]


# function which calculates Pielou's evenness (from 0 - 1)
pielou <- function(x) {
  y <- x / sum(x)
  h_hat <- - sum(y * log(y))
  h_max <- log(length(y))
  h_hat / h_max
}

# calculate richness and diversity of total and qSIP community
alpha <- asv[, .(total_richness = uniqueN(taxonID),
                      qsip_richness = uniqueN(taxonID[in_qsip]),
                      total_invsimpson = vegan::diversity(rel_abund, 'invsimpson'),
                      qsip_invsimpson = vegan::diversity(rel_abund[in_qsip], 'invsimpson'),
                      total_pielou = pielou(rel_abund),
                      qsip_pielou = pielou(rel_abund[in_qsip])),
                  by = .(sampleID, isotope, treatment, ecosystem, timepoint)
             ][sampleID %in% dat$sample_id
               ][timepoint == 7]

alpha <- merge(alpha, dd, all.x = T, by = 'sampleID')

# run linear model to see if diversity is different between treatments
# remove outlier in growth density dependence relationship
rich_cndd_g <- lm(qsip_richness ~ growth_dd, alpha[growth_dd > -0.02])
even_cndd_g <- lm(qsip_invsimpson ~ growth_dd, alpha[growth_dd > -0.02])

rich_cndd_m <- lm(qsip_richness ~ mortality_dd, alpha)  
even_cndd_m <- lm(qsip_invsimpson ~ mortality_dd, alpha)

ggplot(alpha[growth_dd > -0.02], 
       aes(growth_dd, qsip_richness, color = treatment)) +
  geom_smooth(formula = y ~ x, method = lm, se = F, color = 'black') +
  geom_point() +
  xlab('Growth density dependence') +
  ylab('ASV richness') +
  scale_color_manual(values=hsv(c(0, .55, .125), c(0, 1, 1), c(.5, .6, .8)))

ggplot(alpha[growth_dd > -0.02], 
       aes(growth_dd, qsip_invsimpson, color = treatment)) +
  geom_smooth(formula = y ~ x, method = lm, se = F, color = 'black') +
  geom_point() +
  xlab('Growth density dependence') +
  ylab('Simpson evenness') +
  scale_color_manual(values=hsv(c(0, .55, .125), c(0, 1, 1), c(.5, .6, .8)))

ggplot(alpha, 
       aes(mortality_dd, qsip_richness, color = treatment)) +
  geom_smooth(formula = y ~ x, method = lm, se = F, color = 'black') +
  geom_point() +
  xlab('Mortality density dependence') +
  ylab('ASV richness') +
  scale_color_manual(values=hsv(c(0, .55, .125), c(0, 1, 1), c(.5, .6, .8)))

ggplot(alpha, 
       aes(mortality_dd, qsip_invsimpson, color = treatment)) +
  geom_smooth(formula = y ~ x, method = lm, se = F, color = 'black') +
  geom_point() +
  xlab('Mortality density dependence') +
  ylab('Simpson evenness') +
  scale_color_manual(values=hsv(c(0, .55, .125), c(0, 1, 1), c(.5, .6, .8)))

