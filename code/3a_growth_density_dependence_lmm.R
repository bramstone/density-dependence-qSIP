#------------
#   Stone et al. 
#   2022
#   Nutrients strengthen density dependence of per-capita growth and mortality rates in the soil bacterial community
#   Oecologia
# -----------
# calculation of growth rates using linear mixed effects models
# Following the advice of Ieno et al. "Mixed effects models and extensions in ecology with R"

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
               ][, n_rep := .N, by = .(taxon_id, treatment)
                 ][n_rep > 2
                   ][, n_rep := NULL]

# how many taxa now?
dat[, uniqueN(taxon_id)]
# 524


# Estimate best random structure------------------------------------------------

set.seed(82621)

# No random variation structure
g_init <- lm(growth ~ pop_log10_t0*treatment,
             data = dat)

# just ecosystem
g_ie <- lmer(growth ~ pop_log10_t0*treatment + 
               (1|ecosystem),
             data = dat, REML = F)

g_pe <- lmer(growth ~ pop_log10_t0*treatment + 
               (pop_log10_t0|ecosystem),
             data = dat, REML = F)
# singular fit

# just taxon ID
g_it <- lmer(growth ~ pop_log10_t0*treatment + 
               (1|taxon_id),
             data = dat, REML = F)

g_pt <- lmer(growth ~ pop_log10_t0*treatment + 
               (pop_log10_t0|taxon_id),
             data = dat, REML = F)
# failed to converge when adjusted for 16S GCN

# just Phylum
g_ip <- lmer(growth ~ pop_log10_t0*treatment + 
               (1|Phylum),
             data = dat, REML = F)

g_pp <- lmer(growth ~ pop_log10_t0*treatment + 
               (pop_log10_t0|Phylum),
             data = dat, REML = F)

# taxon ID and ecosystem
g_ite <- lmer(growth ~ pop_log10_t0*treatment +
                (1|taxon_id/ecosystem),
              data = dat, REML = F)

g_pte <- lmer(growth ~ pop_log10_t0*treatment +
                (pop_log10_t0|taxon_id/ecosystem),
              data = dat, REML = F)
# failed to converge

# Phylum, taxon ID, and ecosystem
g_ipte <- lmer(growth ~ pop_log10_t0*treatment + 
                 (1|Phylum/taxon_id/ecosystem),
               data = dat, REML = F)

g_ppte <- lmer(growth ~ pop_log10_t0*treatment + 
                 (pop_log10_t0|Phylum/taxon_id/ecosystem),
               data = dat, REML = F)
# singular

# Phylum, taxon ID, and ecosystem, but separate slopes/intercepts
g_ipt_e <- lmer(growth ~ pop_log10_t0*treatment + 
                 (1|Phylum/taxon_id) +
                 (1|ecosystem),
               data = dat, REML = F)

g_ppt_e <- lmer(growth ~ pop_log10_t0*treatment + 
                 (pop_log10_t0|Phylum/taxon_id) +
                 (pop_log10_t0|ecosystem),
               data = dat, REML = F)
# failed to converge

# compare, keeping lm object at the end, excluding models that didn't converge or had singular fit
mod_compare <- anova(g_ie, 
                     g_pe,  # this is the model specified to answer H3
                     g_it,
                     # g_pt,
                     g_ip,
                     g_pp,
                     g_ite, 
                     # g_pte, 
                     g_ipte, 
                     # g_ppte,  
                     g_ipt_e, 
                     # g_ppt_e,  # leaving this out for parity with mortality model
                     g_init)

mod_compare <- data.table(model = rownames(mod_compare), mod_compare)
mod_compare[, delt_AIC := AIC - min(AIC)][order(delt_AIC)][]


# compare fixed effects
pull_terms <- function(x, aic) {
  y <- anova(get(x))
  z <- c(y$`F value`, y$`Pr(>F)`, aic)
  term <- c(paste0(rownames(y), '_F'), paste0(rownames(y), '_p'), 'delt_AIC')
  list(val = z, term = term)
}

fe_compare <- mod_compare[, delt_AIC := AIC - min(AIC)
                          ][order(delt_AIC), .(model, delt_AIC)
                            ][, pull_terms(model, delt_AIC), by = model]
fe_compare <- dcast(fe_compare[!is.na(val)], model ~ term, value.var = 'val')

fe_compare[order(delt_AIC), .(model, 
                              pop_log10_t0_F, pop_log10_t0_p,
                              treatment_F, treatment_p,
                              `pop_log10_t0:treatment_F`, `pop_log10_t0:treatment_p`)]


# Check residuals---------------------------------------------------------------

# Phylum, taxon ID, and ecosystem, but separate slopes/intercepts
g_optim <- lmer(growth ~ pop_log10_t0*treatment + 
                  (1|Phylum/Class/taxon_id) +
                  (1|ecosystem),
                data = dat, REML = F)

# get residuals
resid_data <- cbind(g_optim@frame, 
                    resid=resid(g_optim, scaled = T),
                    fitted=fitted(g_optim))

# residuals by fitted values, grouped by treatment
ggplot(resid_data, aes(fitted, resid)) + 
  geom_point(pch=1) +
  geom_hline(yintercept=0, color='red') +
  facet_wrap(~ treatment, ncol=2) +
  xlab('Fitted values') +
  ylab('Normalized residuals')

# residuals by pop
ggplot(resid_data, aes(pop_log10_t0, resid)) + 
  geom_point(pch=1) +
  geom_hline(yintercept=0, color='red') +
  facet_wrap(~ treatment, ncol=2) +
  xlab('Log10 population density') +
  ylab('Normalized residuals')

# overall OK, but...
# there is still a bit of trend in the fitted ~ residual plot


# Estimate best fixed structure-------------------------------------------------

# fit full model - include weighted variance of density as a covariate
g_full <- lmer(growth ~ pop_log10_t0*treatment + 
                 (1|Phylum/Class/taxon_id) +
                 (1|ecosystem),         
               data = dat, REML = F)

# remove 2-way interactions and compare with Log-likelihood test
g_no_2_int <- lmer(update(g_full, . ~ . - pop_log10_t0:treatment), data = dat, REML = F)
anova(g_full, g_no_2_int)
# full model is best

# fit final model with REML = T
g_optim <- lmer(growth ~ pop_log10_t0*treatment + 
                  (1|Phylum/taxon_id) +
                  (1|ecosystem),         
                data = dat, REML = T)

# look at significant terms
anova(g_optim)
# All terms significant

# plot
dat$pred_growth <- predict(g_optim, newdata = dat, re.form = NA)

log_breaks <- outer(1:10, 10^(1:5), '*')
log_labels <- log_breaks
log_labels[] <- scales::label_number_si()(log_labels)
log_labels[2:nrow(log_labels),] <- ''

ggplot(dat, aes(pop_log10_t0, growth)) +
  geom_hline(yintercept = 0, color = gray(.5)) +
  geom_point(size=rel(.75), color=hsv(.55, .6, .35)) +
  geom_line(aes(y = pred_growth), color=hsv(0,.8,1), size=.8) +
  xlab(expression(atop('Initial population density', '(genomes g'^-1*' soil, log'[10]*'-scale)'))) +
  ylab('Per-capita growth rate') +
  scale_x_continuous(breaks = c(log10(log_breaks)),
                     labels = c(log_labels)) +
  facet_grid(. ~ treatment)
