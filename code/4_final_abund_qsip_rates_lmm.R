#------------
#   Stone et al. 
#   2022
#   Nutrients strengthen density dependence of per-capita growth and mortality rates in the soil bacterial community
#   Oecologia
# -----------
# relationship between final abundances and growth/mortality

library(data.table)
library(lme4)
library(lmerTest)
library(ggplot2)

dat <- fread('../data/qsip_growth_mortality.csv')
gcn <- fread('../data/qsip_taxa_16s_gcn_nstd.csv')

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


# test relationship between growth and final abundances
set.seed(82821)
g_lme_7 <- lmer(pop_log10 ~ growth*treatment + 
                  (1|Phylum/taxon_id/ecosystem),
                data = dat, REML = T)

# test relationship between mortality and final abundances
m_lme_7 <- lmer(pop_log10 ~ mortality_corrected*treatment + 
                  (1|Phylum/taxon_id/ecosystem),
                data = dat, REML = T)

dat$pred_growth <- predict(g_lme_7, newdata = dat, re.form = NA)
dat$pred_mortality <- predict(m_lme_7, newdata = dat, re.form = NA)

anova(g_lme_7)
anoga(m_lme_7)


# plot relationship between final abundances and growth
ggplot(dat, aes(growth, pop_log10)) +
  geom_vline(xintercept = 0, color = gray(.5)) +
  geom_point(size=rel(.75), color=hsv(.55, .6, .35)) +
  geom_line(aes(y = pred_growth), color=hsv(0,.8,1), size=.8) +
  ylab(expression(atop('Final population density', '(genomes g'^-1*' soil, log'[10]*'-scale)'))) +
  xlab('Per-capita growth rate') +
  facet_grid(. ~ treatment)

ggplot(dat, aes(mortality, pop_log10)) +
  geom_vline(xintercept = 0, color = gray(.5)) +
  geom_point(size=rel(.75), color=hsv(.55, .6, .35)) +
  geom_line(aes(y = pred_mortality), color=hsv(0,.8,1), size=.8) +
  ylab(expression(atop('Final population density', '(genomes g'^-1*' soil, log'[10]*'-scale)'))) +
  xlab('Per-capita mortality rate') +
  facet_grid(. ~ treatment)
