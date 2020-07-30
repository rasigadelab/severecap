#'#######################################################################
#'ANALYSIS CODE FOR
#'Risk factors of severity in community-acquired staphylococcal pneumonia
#'Gillet et al.
#'
#'(c) 2020, Jean-Philippe Rasigade, University of Lyon
#'jean-philippe.rasigade<at>univ-lyon1.fr

library(data.table, quietly = T)
library(lubridate, quietly = T)
library(stringr, quietly = T)
library(survival)

load("imputed.Rdata")

##############################################################
# Tabulate raw data (not imputed data)
# Toddlers, PVL, MRSA, death

table(d$death)
table(d[pvl == FALSE]$death)

tb <- table(d[, .(toddler, pvl)])
tb
fisher.test(tb)

tb <- table(d[, .(toddler, mrsa)])
tb
fisher.test(tb)

tb <- table(d[, .(toddler, death)])
tb
fisher.test(tb)

table(d[toddler&pvl == TRUE]$death)
table(d[(!toddler)&pvl == TRUE]$death)

##############################################################
# NON-TODDLER ANALYSES

d3 <- d[toddler == FALSE]

# Treatment
tb <- table(d3[, .(ther_adapt, mrsa)])
tb
fisher.test(tb)

tb <- table(d3[, .(ther_adapt, death)])
tb
fisher.test(tb)

tb <- table(d3[, .(ther_antitox, death)])
tb
fisher.test(tb)

summary(glm(death ~ pvl + mrsa, d3, family = "binomial"))
summary(glm(death ~ pvl + mrsa + ther_antitox + ther_adapt, d3, family = "binomial"))

##############################################################
# MEDIAN SURVIVAL TIMES

surv <- Surv(d3$length_of_stay, event = d3$death, type = "right")

medsurv <- d3[, .(length_of_stay, death, pvl)][
  , death_delay := length_of_stay
  ][death == FALSE, death_delay := Inf]

# QUANTILES AND WILCOXON

medsurv[death == TRUE, .(
  q1 = quantile(death_delay, 0.25),
  q2 = quantile(death_delay, 0.50),
  q3 = quantile(death_delay, 0.75)
), by = pvl]

wilcox.test(death_delay ~ pvl, medsurv[death == TRUE])

medsurv[death == FALSE, .(
  q1 = quantile(length_of_stay, 0.25),
  q2 = quantile(length_of_stay, 0.50),
  q3 = quantile(length_of_stay, 0.75)
), by = pvl]

wilcox.test(length_of_stay ~ pvl, medsurv[death == FALSE])

#############################################################################???
#' ## Survival curves, all patients
library(survminer, quietly = T)

surv <- Surv(d$length_of_stay, event = d$death, type = "right")

#' ### Global curves
ggsurvplot(survfit(surv ~ mrsa, d), risk.table = TRUE, xlim = c(0,90), title = "MRSA influence, all patients")
ggsurvplot(survfit(surv ~ pvl, d), risk.table = TRUE, xlim = c(0,90), title = "PVL influence, all patients")


#' ### Stratified curves, per age group
#' 
# Combine strata into new factor "S" for readability
marks <- c("-", "+")
d[ , S := factor(paste("PVL",  marks[pvl + 1], " ", "MRSA", marks[mrsa + 1], sep = ""))]
ggsurvplot(survfit(surv ~ S, d, subset = age <= 5), risk.table = F, title = "Age 0-5y")
ggsurvplot(survfit(surv ~ S, d, subset = age > 5 & age <= 20), risk.table = F, title = "Age 6-20y")
ggsurvplot(survfit(surv ~ S, d, subset = age > 20 & age <= 60), risk.table = F, xlim = c(0,90), title = "Age 21-60y")
ggsurvplot(survfit(surv ~ S, d, subset = age > 60), risk.table = F, xlim = c(0,90), title = "Age 61-84y")

##############################################################
# SURVIVAL MODELS

# Non-toddlers only, use imputed data to avoid row-wise deletion
dimp3 <- dimp[toddler == FALSE]

surv <- Surv(dimp3$length_of_stay, event = dimp3$death, type = "right")

# PVL and MRSA
mod <- coxph(surv ~ pvl + mrsa, dimp3)
mod
exp(confint(mod))

coxph(surv ~ pvl, dimp3)
coxph(surv ~ mrsa, dimp3)

# Effect of baseline covariates
coxph(surv ~ charlson + sex, dimp3)
coxph(surv ~ pvl + mrsa, dimp3)
coxph(surv ~ pvl + mrsa + charlson + sex, dimp3)
anova(coxph(surv ~ pvl + mrsa, dimp3), coxph(surv ~ pvl + mrsa + charlson + sex, dimp3))

# Effect of treatment
coxph(surv ~ ther_adapt, dimp3)
coxph(surv ~ ther_adapt, dimp3, subset = mrsa == TRUE)
coxph(surv ~ ther_adapt, dimp3, subset = mrsa == FALSE)

coxph(surv ~ ther_antitox, dimp3)
coxph(surv ~ ther_antitox, dimp3, subset = pvl == TRUE)
coxph(surv ~ ther_antitox, dimp3, subset = pvl == FALSE)

#################################################################
# MULTIVARIATE SURVIVAL MODEL - TABLE 3

# Log-transformed predictors
dimp3[, pct_0_l2 := log2(pct_0)]
dimp3[, lactates_0_l2 := log2(lactates_0)]

# List of predictors
preds <- c("charlson", "sex", "pvl", "mrsa", "sofa_0", "flu", "hemoptysis_0", "rash_0", "leukopenia",
           "pct_0_l2", "lactates_0_l2", "ther_adapt", "ther_antitox")

# BIVARIATE MODELS
modlist <- data.table(name = preds, t(sapply(preds, function(pred) {
  form <- as.formula(sprintf("surv ~ %s", pred))
  mod <- coxph(form, dimp3)
  exp(cbind(coef(mod), confint(mod)))
})))

# setnames(modlist, names(modlist), c("name", "coef", "expcoef", "se", "z", "p"))
setnames(modlist, names(modlist), c("name", "coef", "lo", "hi"))

modlist

modlist[ , HR_univ := sprintf("%.2f (%.2f to %.2f)", coef, lo, hi)]

# FULL MULTIVARIATE MODEL
fullmod <- coxph(surv ~ ., dimp3[, preds, with = F])
anova(fullmod)
extractAIC(fullmod)

fullmodlist <- data.table(name = preds, exp(cbind(coef = coef(fullmod), confint(fullmod))) )
setnames(fullmodlist, names(fullmodlist), c("name", "coef", "lo", "hi"))

fullmodlist[ , HR_full := sprintf("%.2f (%.2f to %.2f)", coef, lo, hi)]

# STEPWISE MULTIVARIATE MODEL
stepmod <- step(fullmod)
extractAIC(stepmod)

stepmodlist <- data.table(name = names(coef(stepmod)), exp(cbind(coef = coef(stepmod), confint(stepmod))) )
setnames(stepmodlist, names(stepmodlist), c("name", "coef", "lo", "hi"))

stepmodlist[ , HR_step := sprintf("%.2f (%.2f to %.2f)", coef, lo, hi)]

# COMBINE RESULTS AND BUILD TABLE 3
table3 <- modlist[, .(name, HR_univ)]
table3 <- merge(table3, fullmodlist[, .(name, HR_full)], all.x = TRUE, sort = FALSE)
table3 <- merge(table3, stepmodlist[, .(name, HR_step)], all.x = TRUE, sort = FALSE)

table3

####################################################################################
# PREDICTOR INTERACTION MATRIX
# Use Z-scores for common scale

modlist <- list()

for(pred in preds) {
  form <- as.formula(sprintf("~ . - %s", pred))
  mod <- update(fullmod, form)
  modlist[[pred]] <- mod
}

z <- data.table(name = preds, full = coef(summary(fullmod))[, "z"])

for(pred in preds) {
  cf <- coef(summary(modlist[[pred]]))[, "z"]
  zpred <- data.table(name = names(cf), zscore = cf)
  setnames(zpred, "zscore", pred)
  z <- merge(z, zpred, by = "name", all.x = TRUE, sort = FALSE)
}

print(z)
# xlclipboard(z)

# Take relative differences
zdiff <- data.table(z)
zdiff[ , (preds) := lapply(.SD, function(x) x - full), .SDcol = preds]
zdiff

# xlclipboard(zdiff)

##########################################################################
# INTERACTION MODEL - TABLE 4

# Add selected interactions suggested by the pairwise interaction analysis
fullmod <- coxph(surv ~ 
              charlson + sex + pvl + mrsa + sofa_0 + flu + hemoptysis_0 + rash_0 + leukopenia + 
              pct_0_l2 + lactates_0_l2 + ther_antitox + ther_adapt +
              ther_adapt:mrsa +
              sofa_0:lactates_0_l2 +
              flu:leukopenia +
              ther_antitox:pvl
            , dimp3)

extractAIC(fullmod)
fullmodlist <- data.table(name = names(coef(fullmod)), exp(cbind(coef = coef(fullmod), confint(fullmod))) )
setnames(fullmodlist, names(fullmodlist), c("name", "coef", "lo", "hi"))
fullmodlist[ , HR_full := sprintf("%.2f (%.2f to %.2f)", coef, lo, hi)]

stepmod <- step(fullmod)
extractAIC(stepmod)

stepmodlist <- data.table(name = names(coef(stepmod)), exp(cbind(coef = coef(stepmod), confint(stepmod))) )
setnames(stepmodlist, names(stepmodlist), c("name", "coef", "lo", "hi"))
stepmodlist[ , HR_step := sprintf("%.2f (%.2f to %.2f)", coef, lo, hi)]

table4 <- fullmodlist[, .(name, HR_full)]
table4 <- merge(table4, stepmodlist[, .(name, HR_step)], all.x = TRUE, sort = FALSE)

table4

# xlclipboard(table4)


