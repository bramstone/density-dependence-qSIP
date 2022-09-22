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


# Estimate best random structure------------------------------------------------

set.seed(82521)

# No random variation structure
m_init <- lm(mortality ~ pop_log10_t0*treatment,
             data = dat)

# just ecosystem
m_ie <- lmer(mortality ~ pop_log10_t0*treatment + 
               (1|ecosystem),
             data = dat, REML = F)

m_pe <- lmer(mortality ~ pop_log10_t0*treatment + 
               (pop_log10_t0|ecosystem),
             data = dat, REML = F)
# singular when adjusting for 16S GCN

# just taxon ID
m_it <- lmer(mortality ~ pop_log10_t0*treatment + 
               (1|taxon_id),
             data = dat, REML = F)

m_pt <- lmer(mortality ~ pop_log10_t0*treatment + 
               (pop_log10_t0|taxon_id),
             data = dat, REML = F)

# just Phylum
m_ip <- lmer(mortality ~ pop_log10_t0*treatment + 
               (1|Phylum),
             data = dat, REML = F)

m_pp <- lmer(mortality ~ pop_log10_t0*treatment + 
               (pop_log10_t0|Phylum),
             data = dat, REML = F)
# failed to converge

# taxon ID and ecosystem
m_ite <- lmer(mortality ~ pop_log10_t0*treatment +
                (1|taxon_id/ecosystem),
              data = dat, REML = F)

m_pte <- lmer(mortality ~ pop_log10_t0*treatment +
                (pop_log10_t0|taxon_id/ecosystem),
              data = dat, REML = F)
# singular fit - exclude - failed to converge when adjusting for 16S GCN

# Phylum, taxon ID, and ecosystem
m_ipte <- lmer(mortality ~ pop_log10_t0*treatment + 
                 (1|Phylum/taxon_id/ecosystem),
               data = dat, REML = F)

m_ppte <- lmer(mortality ~ pop_log10_t0*treatment + 
                 (pop_log10_t0|Phylum/taxon_id/ecosystem),
               data = dat, REML = F)
# singular fit

# Phylum, taxon ID, and ecosystem, but separate slopes/intercepts
m_ipt_e <- lmer(mortality ~ pop_log10_t0*treatment + 
                 (1|Phylum/taxon_id) +
                 (1|ecosystem),
               data = dat, REML = F)

m_ppt_e <- lmer(mortality ~ pop_log10_t0*treatment + 
                 (pop_log10_t0|Phylum/taxon_id) +
                 (pop_log10_t0|ecosystem),
               data = dat, REML = F)
# failed to converge

# compare, keeping lm object at the end, excluding models that didn't converge or had singular fit
mod_compare <- anova(m_ie, 
                     m_pe,  # this is the model specified to address H3
                     m_it,
                     m_pt,
                     m_ip,
                     # m_pp,
                     m_ite, 
                     # m_pte, 
                     m_ipte,
                     # m_ppte, 
                     m_ipt_e,
                     # m_ppt_e, 
                     m_init)

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

# get residuals
resid_data <- cbind(m_ipt_e@frame, 
                    resid=resid(m_ipt_e, scaled = T),
                    fitted=fitted(m_ipte))

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

# residuals look pretty good - but there is slightly lower variance in the C treatment


# Estimate best fixed structure-------------------------------------------------

# fit full model
m_full <- lmer(mortality ~ pop_log10_t0*treatment + 
                 (1|Phylum/taxon_id) +
                 (1|ecosystem),
               data = dat, REML = F)

# remove 2-way interaction and compare with Log-likelihood test
m_no_2_int <- lmer(update(m_full, . ~ . - pop_log10_t0:treatment), data = dat, REML = F)
anova(m_full, m_no_2_int)
# better with the 2-way interaction 

# fit final model with REML = T
m_optim <- lmer(mortality ~ pop_log10_t0*treatment + 
                  (1|Phylum/taxon_id) +
                  (1|ecosystem),
                data = dat, REML = T)

# look at significant terms
anova(m_optim)
# All terms significant


# plot
dat$pred_mortality <- predict(m_optim, newdata = dat, re.form = NA)

log_breaks <- outer(1:10, 10^(1:5), '*')
log_labels <- log_breaks
log_labels[] <- scales::label_number_si()(log_labels)
log_labels[2:nrow(log_labels),] <- ''

ggplot(dat, aes(pop_log10_t0, mortality)) +
  geom_hline(yintercept = 0, color = gray(.5)) +
  geom_point(size=rel(.75), color=hsv(.55, .6, .35)) +
  geom_line(aes(y = pred_mortality), color=hsv(0,.8,1), size=.8) +
  xlab(expression(atop('Initial population density', '(genomes g'^-1*' soil, log'[10]*'-scale)'))) +
  ylab('Per-capita mortality rate') +
  scale_x_continuous(breaks = c(log10(log_breaks)),
                     labels = c(log_labels)) +
  facet_grid(. ~ treatment)
