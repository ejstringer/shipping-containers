

# description -------------------------------------------------------------

## from this exploration of the data I show that Container.sample.ID from the
## shipping data links to the metabarcoding data, while the Container.number
## links to the species specific data.
##     I also show that there is a mismatch, where the shipping data contains a 
## sample id that has a .1 and .2 at the end, for sample 20210714_1354. My plan
## is to rename 20210714_1354.1  to 20210714_1354 as 20210714_1354.2 does not 
## have any associated eDNA data. 
##     My next step is to simplify the files by only including containers with
## genetic data associated with them, filtering further the arthropoda metabar-
## coding data and joining dna and rna data, and joining the species specific 
## data. This code will be in a file called simplifying_data_files.R

# load --------------------------------------------------------------------
source('./code/libraries.R')

ship <- read.csv('./data/shipping_meta/Additional Diagnostics Extract.csv')

dna <- read.csv('./output/ASV_all_tax_count_arthropoda.csv')

rna <- read.csv('./data/metabarcoding/ASV_cdna_tax_count.tsv')

spp <- read.csv('./data/species_specific/Khapra beetle_2023.csv')

table(spp$Method.Collection)

# match containers --------------------------------------------------------

containers <- ship %>% 
  select(Container.sample.ID, Container.number) %>% 
  separate(col = Container.sample.ID, sep = '_', remove = FALSE,
           into = c('con.no', 'v', 'c05', 'id1', 'id2')) %>%
  filter(complete.cases(id1)) %>% 
  mutate(id = paste(id1, id2, sep = '_')) %>% 
  unique()

table(unique(spp$Container.Number) %in% containers$Container.number)

containers$Container.sample.ID
dna_container_id <- colnames(dna)[grep('X20', colnames(dna))] %>% 
  sub('X', '', .)

xiny <- dna_container_id %in% containers$id
table(xiny)

which(dna_container_id == dna_container_id[!xiny])



yinx <- containers$id %in% dna_container_id
table(yinx)
containers$id[!yinx]
which(containers$id %in% containers$id[!yinx])


mismatch_position <- grep(dna_container_id[!xiny], containers$id)
mismatach_container <- containers$Container.sample.ID[mismatch_position]

missing <- ship %>% filter(Container.sample.ID %in% mismatach_container) %>% 
  mutate(id = 1:nrow(.)) %>% 
  select(Container.sample.ID, Arrival.date, 
         Loading.country, Destination.country, Origin.country, Vacuum, Is.Khapra.eDNA.positive, Is.Khapra.eRNA.positive
         ) %>% 
  pivot_wider(names_from = Container.sample.ID,
              values_from = c(Arrival.date, Loading.country, Destination.country, Origin.country,
                              Vacuum, Is.Khapra.eDNA.positive, Is.Khapra.eRNA.positive))
missing$
for (i in seq(1,length(missing), 2)) {
  x <- sort(missing[[i]][[1]]) %in% sort(missing[[i+1]][[1]]) 
  print(table(x))
}

ship$Container.sample.ID[grep('\\.1', ship$Container.sample.ID)]

table(ship$Sample.method, useNA = 'always') 
table(ship$Vacuum, useNA = 'always')

