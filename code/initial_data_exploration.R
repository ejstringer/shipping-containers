

# description -------------------------------------------------------------

## from this exploration of the data I show that Container.sample.ID from the
## shipping data links to the metabarcoding data, while the Container.number
## links to the species specific data.
##     I also show that there is a mismatch, where the shipping data contains a 
## sample id that has a .1 and .2 at the end, for sample 20210714_1354. My plan
## is to rename 20210714_1354.1  to 20210714_1354 as 20210714_1354.2 does not 
## have any associated eDNA data. Actually, 20210714_1354.2 does not seem to 
## have any data except for shipping data associated with it and a different 
## sample weight... I might just remove this, however it is in all the species
## specific data. 
##     My next step is to simplify the files by only including containers with
## genetic data associated with them, filtering further the arthropoda metabar-
## coding data and joining dna and rna data, and joining the species specific 
## data. This code will be in a file called curate_datasets.R

## Ok so Container.Sample.Id is just a paste(), the lab.code is the true simple
## sample id, which I will use as the index for my data frame links. 

## I have identify multiple data frames to separate the shipping container data
## into. 

# load --------------------------------------------------------------------
source('./code/libraries.R')

ship <- read.csv('./data/shipping_meta/Additional Diagnostics Extract.csv')

dna <- read.csv('./output/ASV_all_tax_count_arthropoda.csv')

rna <- read.table('./data/metabarcoding/ASV_cdna_tax_count.tsv', 
                  header = T, sep = '\t')

spp <- read.csv('./data/species_specific/Khapra beetle_2023.csv')
ant <- read.csv('./data/species_specific/Electric Ant_2023.csv')
bmsb <- read.csv('./data/species_specific/BMSB_2023.csv')

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
mismatach_container

missing <- ship %>% filter(Container.sample.ID %in% mismatach_container) %>% 
  mutate(id = 1:nrow(.)) %>% 
  select(Container.sample.ID, Arrival.date, 
         Loading.country, Destination.country, Origin.country, Vacuum, Is.Khapra.eDNA.positive, Is.Khapra.eRNA.positive
         ) %>% 
  pivot_wider(names_from = Container.sample.ID,
              values_from = c(Arrival.date, Loading.country, Destination.country, Origin.country,
                              Vacuum, Is.Khapra.eDNA.positive, Is.Khapra.eRNA.positive))

for (i in seq(1,length(missing), 2)) {
  x <- sort(missing[[i]][[1]]) %in% sort(missing[[i+1]][[1]]) 
  print(table(x))
}

ship$Container.sample.ID[grep('\\.1', ship$Container.sample.ID)]

table(ship$Sample.method, useNA = 'always') 
table(ship$Vacuum, useNA = 'always')

spp[grep('\\.', spp$Unique.Sample.Identifier),]

## find data difference in .1 and .2
the.samples <- ship %>% filter(grepl('\\.', Lab.code))

table(the.samples$Lab.code, the.samples$Identification.type)
table(the.samples$Lab.code, the.samples$Sample.method)
table(the.samples$Lab.code, the.samples$Sample.weight..g.)
table(the.samples$Identification.type, the.samples$Sample.weight..g.)


the.samples %>% filter(grepl('\\.2', Lab.code)) %>% View

# simplify ----------------------------------------------------------------
ship_simple <- ship %>% 
  filter(Container.number %in% ant$Container.Number | 
           Container.number %in% bmsb$Container.Number) %>% 
  mutate(sample_id = sub('C05246_', 'X', Lab.code),
         sample_id = ifelse(sample_id == '', 
                            paste0(Container.sample.ID, 'xx'), sample_id)) %>% 
  relocate(c(Identification.type, Original.container.sample.ID),
           .after = Container.number) %>% 
  relocate(sample_id)

# shipping ----------------------------------------------------------------

spp$Unique.Sample.Identifier %>% unique %>% length
spp$Container.Number %>% unique %>% length

shipx <- ship %>% mutate(Container.sample.ID = sub('\\.1', '', Container.sample.ID))

shipx$Container.sample.ID[grep('\\.1', ship$Container.sample.ID)] %>% length

table(shipx$Container.sample.ID == ship$Container.sample.ID)

sort(table(ship$Container.number))[1000]
ship %>% names
ship[,-c(grep('DNA', names(ship)),grep('RNA', names(ship)),
         which(names(ship) %in% c('Vacuum', 'Sweep', 'Vacuum.underfloor',
                                  'Sample.method.code', 'Project.code',
                                  'No.result', 'Has.undetermined.result',
                                  'Original.container.sample.ID')))] %>%
  filter(Container.number == 'FTAU1272657') %>% 
  relocate(Identification.type, .after = Container.number) %>% View

ship$Is.Khapra.eDNA.positive[which(ship$Container.sample.ID != ship$Original.container.sample.ID)]

table(ship$Container.number, ship$Arrival.sequence.number)

ship %>% 
  filter(Container.number)
arrange(container_id, Arrival.sequence.number) %>% 
  select(container_id, Arrival.sequence.number, 
         Arrival.date, Destination.country) %>% View
ship_simple$Identification.type %>% table
ship_simple %>% head(., 2)
ship_simple$Is.latest.arrival %>% table

test <- ship_simple %>% 
  select(Container.number, Arrival.date, Is.latest.arrival) %>% 
  group_by(Container.number) %>% 
  mutate(notlatest = isFALSE(Is.latest.arrival)) %>% 
  summarise(maxdate = max(Arrival.date),
            latest = sum(notlatest)) 
test$latest %>% sum

table(ship$Container.number %in% spp$Container.Number)
table(spp$Container.Number %in% ship$Container.number)

ship$Container.number[!ship$Container.number %in% spp$Container.Number]


## metabarcoding containers
dna_container_id <- colnames(dna)[grep('X20', colnames(dna))]
dna_container_id[grep('20210714_1354', dna_container_id)]

rna_container_id <- colnames(dna)[grep('X20', colnames(rna))]
rna_container_id[grep('20210714_1354', dna_container_id)]

colnames(dna)[1:25]  
colnames(dna)[1:25]  %in% sub('_cDNA', '', colnames(rna))
sub('_cDNA', '', colnames(rna)) %in% colnames(dna)
colnames(rna)[!sub('_cDNA', '', colnames(rna)) %in% colnames(dna)]


unique(ship_simple$container_id) %>% table %>% table
dna_container_id %>% length

spp$Container.Number %>% unique() %>% length
ship_simple$Container.number %>% unique() %>% table %>% table

ship_simple %>% names
ship_simple[, c(2,5, 12, 55, 62,63)] %>% 
  unique() %>% 
  arrange(Container.number, Arrival.sequence.number) %>% View
uniqship <- unique(ship_simple[, c(1:2,8, 9, 12,26, 64)])

duplicated(uniqship$Container.number) %>% table
ship %>% arrange(Identification.type) %>% head
dupsship <- uniqship %>%
  filter(Container.number %in% uniqship$Container.number[duplicated(uniqship$Container.number)]) %>% 
  mutate(#dup = duplicated(Container.number),
    inDna = container_id %in% dna_container_id) %>% 
  arrange(Container.number) #%>% 
#split(., .$Container.number)

table(dupsship$container_id %in% dna_container_id)

dupsship %>% split(., .$Container.number)

table(ship_simple$Identification.type, ship_simple$Sample.method.code, useNA = 'always')

ship$Lab.code %>% table

ship_simple %>% filter(Lab.code == '') %>% 
  select(sample_id, Lab.code, Container.sample.ID, Sample.method,Sample.weight..g., Lab.name,
         Identification.type) %>% unique %>% 
  split(., .$Container.sample.ID)

## metabarcoding -----

dna_container_id <- colnames(dna)[grep('X20', colnames(dna))]
dna_container_id[grep('20210714_1354', dna_container_id)]

rna_container_id <- colnames(dna)[grep('X20', colnames(rna))]
rna_container_id[grep('20210714_1354', dna_container_id)]

colnames(dna)[1:25]  
colnames(dna)[1:25]  %in% sub('_cDNA', '', colnames(rna))
sub('_cDNA', '', colnames(rna)) %in% colnames(dna)
colnames(rna)[!sub('_cDNA', '', colnames(rna)) %in% colnames(dna)]

rna_container_id %in% ship_simple$sample_id 
table(dna_container_id %in% ship_simple$sample_id)
dna_container_id[!dna_container_id %in% ship_simple$sample_id ]

table(unique(ship_simple$sample_id) %in% dna_container_id)


## ship extra containers ------
ship[!(ship$Container.number %in% ant$Container.Number | 
  ship$Container.number %in% bmsb$Container.Number),]

# goods -------------------------------------------------------------------


ship_simple %>% 
arrange(container_id, Arrival.sequence.number) %>% 
  select(container_id, Arrival.sequence.number, 
         Arrival.date, Goods.description, 
         Shipping.company.goods.description,
         Goods.khapra.classification,
         Goods.hitchhiker.risk.classification) %>% View
