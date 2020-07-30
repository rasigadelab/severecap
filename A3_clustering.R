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
library(readxl, quietly = T)

load("imputed.Rdata")

##########################################################################
#' ## Model-based clustering
#' 
#' The technique objectivates the presence of clusters of patients grouped by age
#' 
#' ### All patients
library(mclust)
mc <- Mclust(d[]$age_days)
summary(mc)
plot(mc, what = "classification")
plot(mc, what = "density")

#' There's a clearly separated cluster of very young patients (<=3y of age), then 2 clusters, below and above 15000d (40y)
#'
#' ### PVL-positive
#' 
mc <- Mclust(d[pvl == T]$age_days)
summary(mc)
plot(mc, what = "classification")
plot(mc, what = "density")

#' In PVL-positive patients then is no separation in adult patients, the 2 clusters separate around 3-5y

#' ### PVL-negative
#' 
mc <- Mclust(d[pvl == F]$age_days)
summary(mc)
plot(mc, what = "classification")
plot(mc, what = "density")


##########################################################################
# Visualize mortality as a function of age (smoothing procedure)

# Use KDE smoothing
toddler <- FALSE

bw <- ifelse(toddler, 1/4, 8)
dsub <- if(toddler) d[toddler == TRUE][] else d[]

# Subselect PVL+/- cases and rerun analysis if needed
# dsub <- d[pvl == T]

age_all  <- dsub$age_days / 365
age_dead <- dsub[death == T]$age_days / 365

age_left  <- 0
age_right <- max(age_all)

kde_points <- 100
kde_dead <- density(age_dead, n = kde_points, bw = bw, from = age_left, to = age_right)
kde_all  <- density(age_all, n = kde_points, bw = bw, from = age_left, to = age_right)

mortality_x <- kde_all$x
dens_y_all  <- kde_all$y * nrow(dsub)
dens_y_dead <- kde_dead$y * sum(dsub$death)
mortality_y <- dens_y_dead / dens_y_all

# Mortality plots

plot(mortality_x, dens_y_all, type = "l", ylim = c(0, max(dens_y_all)))
lines(mortality_x, dens_y_dead, lty = 2)

plot(mortality_x, mortality_y, type = "l", ylim = c(0, max(mortality_y)))

# NB: confidence bands, see
# https://stats.stackexchange.com/questions/207129/compute-confidence-interval-for-univariate-kernel-density-estimator


# Bootstrap
B <- 200 # FOR TESTING
# B <- 10000 # FOR PUBLICATION
estimates_all   <- matrix(NA, nrow = kde_points, ncol = B)
estimates_dead  <- matrix(NA, nrow = kde_points, ncol = B)
estimates_ratio <- matrix(NA, nrow = kde_points, ncol = B)
b <- 1
while(b <= B) {
  cat(".")
  boot_indices <- sample(1:nrow(dsub), replace = T)
  
  if(sum(dsub[boot_indices]$death)) {
    estimates_all[,b]  <- density(dsub[boot_indices]$age_days / 365, n = kde_points, bw = bw, from = age_left, to = age_right)$y * nrow(dsub)
    estimates_dead[,b] <- density(dsub[boot_indices][death == T]$age_days / 365, n = kde_points, bw = bw, from = age_left, to = age_right)$y * sum(dsub[boot_indices]$death)
    estimates_ratio[,b] <- pmin(1, estimates_dead[,b] / estimates_all[,b])
    b <- b + 1
  }
}

# Mortality plot
ConfidenceBandsAll <- apply(estimates_all, 1, quantile, probs = c(0.025, 0.975))
ConfidenceBandsDead <- apply(estimates_dead, 1, quantile, probs = c(0.025, 0.975))
ConfidenceBandsRatio <- apply(estimates_ratio, 1, quantile, probs = c(0.025, 0.975))

addBand <- function(x, y1, y2, color, alpha = 0.4) {
  xshade <- c(x, rev(x))
  yshade <- c(y2, rev(y1))
  polygon(xshade,yshade, border = NA,col=adjustcolor(color, alpha))
}

# Density plots
{
  par(family = "sans")
  par(mfrow = c(1,2))

  # Patient density plot
  plot(mortality_x, dens_y_all, type = "l", ylim = c(0, max(ConfidenceBandsAll)),
       xlab = "Age (y)", ylab = "No. of patients (density)")
  addBand(mortality_x, ConfidenceBandsAll[1,], ConfidenceBandsAll[2,], "darkgreen")
  
  lines(mortality_x, dens_y_dead)
  addBand(mortality_x, ConfidenceBandsDead[1,], ConfidenceBandsDead[2,], "red")
  
  legend("topleft", bty = "n", fill = c("darkgreen", "red"), legend = c("All", "Dead"))

  # Mortality plot
  plot(mortality_x, mortality_y * 100, type = "l", ylim = c(0, max(ConfidenceBandsRatio)*100),
       xlab = "Age (y)", ylab = "Mortality rate (%)")
  addBand(mortality_x, ConfidenceBandsRatio[1,]*100, ConfidenceBandsRatio[2,]*100, "red")
  par(mfrow = c(1,1))
}


# Density plots (months)
{
  par(family = "sans")
  par(mfrow = c(1,2))
  
  # Patient density plot
  plot(mortality_x * 12, dens_y_all / 12, type = "l", ylim = c(0, max(ConfidenceBandsAll / 12)),
       xlab = "Age (months)", ylab = "No. of patients (density)")
  addBand(mortality_x*12, ConfidenceBandsAll[1,] / 12, ConfidenceBandsAll[2,] / 12, "darkgreen")
  
  lines(mortality_x * 12, dens_y_dead / 12)
  addBand(mortality_x*12, ConfidenceBandsDead[1,] / 12, ConfidenceBandsDead[2,] / 12, "red")
  
  legend("topright", bty = "n", fill = c("darkgreen", "red"), legend = c("All", "Dead"))
  
  # Mortality plot
  plot(mortality_x * 12, mortality_y * 100, type = "l", ylim = c(0, max(ConfidenceBandsRatio)*100),
       xlab = "Age (months)", ylab = "Mortality rate (%)")
  addBand(mortality_x * 12, ConfidenceBandsRatio[1,]*100, ConfidenceBandsRatio[2,]*100, "red")
  par(mfrow = c(1,1))
}




