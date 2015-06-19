## R clean diatom data script
# This script will:
# * Remove other information after epithet (var, etc.) and use naming format:
#    genus.epithet
# * Combine variants specified by sp, spp, sp1, etc., rewrite epithet to sp.,
#    and sum abund_ml across variants
# * Combine diatom variants to epithet level (abund_ml is summed)
# * Do biovol merge * abund_ml on `NLAdiatoms_biovolumes_CD.csv` where only
#    species with both a genus and a non-sp epithet receive a biovolume

# Load required libraries
library(plyr)

# Load data
full_rough_data <- read.csv("../Clean data/full_rough.csv", header=T, sep=",")
diatom_biovols <- read.csv("../Clean data/NLAdiatoms_biovolumes_CD.csv", header=T, sep=",")

# Create diatom dataset
diatom_rough_data <- full_rough_data[full_rough_data$TAXATYPE=='Diatoms',]
taxanames <- diatom_rough_data$TAXANAME

## Standardize all taxanames
# Parsing functions
parse_tax_cf <- function(taxaname) {
  # Removes "cf." from taxaname
  taxaname <- gsub(" cf. ", " ", taxaname)
  taxaname <- gsub(" cff. ", " ", taxaname)
  taxaname
}
parse_tax_sp <- function(taxaname) {
  # Combines sp, spp, sp1, etc. variants
  taxaname <- gsub(" sp. [0-9]", " sp.", taxaname)
  taxaname <- gsub(" spp.", " sp.", taxaname)
  taxaname <- gsub(" sp.[0-9]", " sp.", taxaname)
  taxaname <- gsub("[?]", "", taxaname)
  taxaname
}
# Parse taxaname function
parse_tax <- function(taxaname) {
  # Combines parsing sub-functions
  taxaname <- parse_tax_cf(taxaname)
  taxaname <- parse_tax_sp(taxaname)
  taxaname
}
# Parse all taxanames
taxanames_new <- unlist(sapply(taxanames, FUN = parse_tax))

# Split and remake taxanames
split_join_tax <- function(taxaname) {
  split_string <- strsplit(as.character(taxaname), " ")[[1]]
  taxaname <- paste(split_string[1], split_string[2], sep=".")
  taxaname
}
# Split and join taxanames
taxanames_new <- unlist(sapply(taxanames_new, FUN = split_join_tax))

# Replace diatom taxanames with new taxanames
diatom_rough_data$TAXANAME <- taxanames_new

## Combine diatom row entries by uniqueness
# Subset data.frame and sum abund_ml
diatom_unique <- ddply(diatom_rough_data, .(SITE_ID,VISIT_NO,SAMPLE_CATEGORY,GENUS,TAXANAME,MESH_SIZE,T_GROUP,TAXATYPE,NTL,PTL,LAKE_ORIGIN,LON_DD,LAT_DD,BIOMASS,BIOVOLUME), colwise(sum, c("abund_ml")))
# Match biovolumes from NLAdiatoms_biovolumes_CD.csv
diatom_biovols_index <- match(diatom_unique$TAXANAME, diatom_biovols$species)
unit_biovols <- unlist(sapply(diatom_biovols_index, FUN=function(x) { diatom_biovols$biovolume[x] } ))
# Multiply unit-level biovolumes by abundance
BIOVOLUME <- unit_biovols * diatom_unique$abund_ml
# Add BIOVOLUME to unique diatoms
diatom_unique$BIOVOLUME <- BIOVOLUME

# Reorder diatom data columns
diatom_unique <- diatom_unique[colnames(full_rough_data)]

# Combine data into new full dataset
full_rough_data_new <- full_rough_data[full_rough_data$TAXATYPE!='Diatoms',]
full_rough_data_new <- rbind(full_rough_data_new, diatom_unique)

# Write out dataset to CSV
write.csv(x=full_rough_data_new, file="../Clean data/full_rough_diatoms.csv")