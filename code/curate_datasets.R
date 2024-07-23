
# load --------------------------------------------------------------------
source('./code/libraries.R')

ship <- read.csv('./data/shipping_meta/Additional Diagnostics Extract.csv')

dna <- read.csv('./output/ASV_all_tax_count_arthropoda.csv')

rna <- read.csv('./data/metabarcoding/ASV_cdna_tax_count.tsv')

spp <- read.csv('./data/species_specific/Khapra beetle_2023.csv')


# shipping ----------------------------------------------------------------

ship_simple <- ship %>% 
  select(container_id, Container.number, Arrival.date,
         -Sample.method.code, Sample.method, Vacuum,
         Vacuum.underfloor, Is.Khapra.eDNA.positive, Is.Khapra.eRNA.positive,
         Is.BMSB.eDNA.positive, is.bm)
  separate(col = Container.sample.ID, sep = '_', remove = FALSE,
           into = c('a', 'v', 'b', 'id1', 'id2')) %>%
  mutate(container_id = paste0('X', id1,'_', id2)) %>% 
  select(-a, -v, -b, -id1, -id2, )
  
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

ship_simple %>% 
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
