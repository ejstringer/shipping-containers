
# load --------------------------------------------------------------------
source('./code/libraries.R')

ship <- read.csv('./data/shipping_meta/Additional Diagnostics Extract.csv')

dna <- read.csv('./output/ASV_all_tax_count_arthropoda.csv')

rna <- read.table('./data/metabarcoding/ASV_cdna_tax_count.tsv', 
                  header = T, sep = '\t')

spp_names <- list.files('./data/species_specific/', '.csv')
spp <- lapply(paste0('./data/species_specific/', spp_names), read.csv)
names(spp) <- sub('\\.csv', '', sub('_2023', '', gsub(' ', '_', spp_names)))
spp %>% names


# simplify ship -----------------------------------------------------------


