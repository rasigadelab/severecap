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

# Load clinical data
print(load(file = "data.Rdata"))

# New criterions:
d[, toddler := age <= 3]
d[ , leukopenia := wbc_0 < 3]
d[ , bmi := weight /  (height/100)^2]
d[ , obese := bmi > 30]

# Tabulate
table(d[, .(pvl)])
table(d[, .(pvl)]) / nrow(d)

table(d[, .(toddler)])
table(d[, .(toddler)]) / nrow(d)

table(d[, .(toddler, pvl)])
table(d[, .(toddler, pvl)]) / nrow(d)

########################
# IMPUTE MISSING DATA

library(missRanger)

dimp <- d[, names(d)[-(1:4)], with = F]
for(j in 1:ncol(dimp)) {
  if(class(dimp[[j]]) == "logical") dimp[[j]] <- as.numeric(dimp[[j]])
}
dimp <- missRanger(dimp, num.trees = 100, seed = floor(pi * 1e6))

save(d, dimp, file = "imputed.Rdata")
