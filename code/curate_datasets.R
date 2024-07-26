
# Steps:
## Ship
#   1. sample_id - rename X20210714_1354.1 and remove X20210714_1354.2 from ship
#   2. 
#   3. 

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


# ship --------------------------------------------------------------------
# simplify
ship_simple <- ship %>% 
  filter(Container.number %in% ant$Container.Number | 
           Container.number %in% bmsb$Container.Number) %>% 
  mutate(sample_id = sub('C05246_', 'X', Lab.code),
         sample_id = ifelse(sample_id == '', 
                            paste0(Container.sample.ID, 'xx'), sample_id),
         sample_id = sub('\\.1','', sample_id)) %>% 
  filter(!grepl('\\.2', sample_id)) %>% 
  relocate(c(Identification.type, Original.container.sample.ID),
           .after = Container.number) %>% 
  rename(container_id = Container.number) %>% 
  relocate(sample_id, container_id)  
  
ship_simple$Goods.hitchhiker.risk.classification[ship_simple$Shipping.company.goods.description != ''] %>% table

## collection --------------------------------------------------------------

collection <- ship_simple %>% 
  select(sample_id, container_id, Collection.date, 
         Sample.method, Sample.weight..g.) %>% 
  rename(collection_date = Collection.date,
         sample_method = Sample.method, sample_weight = Sample.weight..g.) %>% 
  unique()


## history -----------------------------------------------------------------

history <- ship_simple %>% 
  mutate(data_origin = ifelse(Goods.description != '' | Origin.country != '',
                             'ICS', 'Shipping company'),
         data_origin = ifelse(Goods.description != '' & 
                                Shipping.company.goods.description != '',
                             'Both', data_origin)) %>% 
  select(container_id, data_origin, Arrival.date, 
         Loading.country, Origin.country,
         Goods.description, Shipping.company.goods.description
         ) %>% 
  filter(Arrival.date != '') %>%
  setNames(tolower(gsub("\\.","_",names(.)))) %>%
  unique() %>% 
  group_by(container_id, arrival_date,data_origin,
           loading_country, origin_country, goods_description) %>% 
  summarise(company_goods_description = paste(shipping_company_goods_description,
                          collapse = ' ;; ')) %>%
  rename(ICS_goods_description = goods_description) %>% 
  ungroup() %>% 
  arrange(container_id, desc(arrival_date)) %>% 
  mutate(sampled = ifelse(duplicated(container_id), F, T)) %>%
  relocate(sampled, container_id, data_origin)

history[1:1000,] %>% View
history$container_id %>% unique %>% length
history$data_origin %>% table


## visual surveys ----------------------------------------------------------


