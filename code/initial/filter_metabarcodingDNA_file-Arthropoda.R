source('./code/libraries.R')


# load data ---------------------------------------------------------------
system.time(dna <- read.delim('./data/metabarcoding/ASV_all_tax_count.tsv',
                              sep = '\t')) # takes roughly 17 minutes to load in


# filter data -------------------------------------------------------------

# only keep arthropods
dnaArthropoda <- dna[dna$Phylum == 'Arthropoda',]

# percent  ----------------------------------------------------------------

all <- nrow(dna)               # total number of asvs
arthrods <- nrow(dnaArthropoda)# arthropoda asvs

round(arthrods/all*100,1) # percent of asvs that are arthropoda


# arthropod data ----------------------------------------------------------
ncol(dnaArthropoda)
names(dnaArthropoda)[1:20]

head(dnaArthropoda)[1:16]


# save --------------------------------------------------------------------

write.csv(dnaArthropoda,
          './output/ASV_all_tax_count_arthropoda.csv',
          row.names = FALSE) # about 3 minutes

system.time(small <- read.csv('./output/ASV_all_tax_count_arthropoda.csv'))
            # 44 seconds to load in
