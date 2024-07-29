
# Steps:

# •	Renamed 20210714_1354.1 to 20210714_1354 and removed 20210714_1354.2
# •	Transshipments made destination country meaningless… but we know the 
#    container was in Australia for sampling
# •	Assume latest arrival date is associated with collection date 
# •	SWEEP method had sample ids that corresponded with different 
#    containers for other methods (VAC and UNDERFLOOR VAC), so added _s to all
#    SWEEP samples (these were not used for any genetic data) - REMOVED SWEEP 
#     data after discussing with Alejandro
# •	Used Furui for khapra beetle data to join with other species-specific data
#    (removed Olson columns). 
# •	Filtered metabarcoding data to include Percent Identity >= 90% and added 
#    a column called high_priority referring to Australia's insect species of 
#    high pest risk.


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

priority <- read.csv('./data/high_priority_insect_species.csv')

quality <- read.csv('./data/C05246_Master metadata_sampledata.csv')

# ship --------------------------------------------------------------------

# simplify
ship_simple <- ship %>% 
  left_join(quality[,c('Unique.Sample.Identifier',
                       'Conc..2..ng.ul.', 'X260.280', 'X260.230',
                       'Working.Dilution')],
            by = c('Lab.code' = 'Unique.Sample.Identifier')) %>% 
  filter(Container.number %in% spp$Electric_Ant$Container.Number | 
           Container.number %in% spp$BMSB$Container.Number) %>% 
  mutate(sample_id = sub('C05246_', 'X', Lab.code),
         sample_id = ifelse(sample_id == '', 
                            paste0(Container.sample.ID, 'xx'), sample_id),
         sample_id = sub('\\.1','', sample_id),
         sample_id = ifelse(Sample.method=='SWEEP',
                            paste0(sample_id, '_s'),
                            sample_id)) %>% 
  filter(!grepl('\\.2', sample_id),
         Sample.method != 'SWEEP') %>% 
  relocate(c(Identification.type, Original.container.sample.ID),
           .after = Container.number) %>% 
  rename(container_id = Container.number) %>% 
  relocate(sample_id, container_id)  
  
## collection --------------------------------------------------------------

collection <- ship_simple %>% 
  select(sample_id, container_id, Collection.date, 
         Sample.method, Sample.weight..g., 
         Conc..2..ng.ul., X260.280, X260.230, Working.Dilution,
         Identification.type, Lab.name) %>% 
  rename(collection_date = Collection.date,
         sample_method = Sample.method,
         sample_weight = Sample.weight..g.,
         visual_contents = Identification.type,
         visual_lab = Lab.name) %>% 
  group_by(sample_id) %>% 
  mutate(visual_contents = ifelse(duplicated(visual_contents),
                                  '',visual_contents)) %>% 
  group_by(sample_id, container_id, collection_date,
           sample_method, sample_weight,
           Conc..2..ng.ul., X260.280, X260.230, Working.Dilution,
           visual_lab) %>% 
  summarise(visual_contents = paste(visual_contents, collapse = '')) %>%
  ungroup %>% 
  mutate(visual_contents = gsub("LC","L;C",visual_contents),
         visual_contents = gsub("LS","L;S",visual_contents),
         visual_contents = gsub("NS","N;S",visual_contents)) %>%
  setNames(str_trim(tolower(gsub("[[:punct:]]+"," ",names(.))))) %>% 
  setNames(gsub(' ', '_', names(.))) %>% 
  unique() 
collection %>% names
table(collection$visual_contents)
table(duplicated(collection$sample_id))
table(collection$sample_method)

## history -----------------------------------------------------------------

history <- ship_simple %>% 
  mutate(data_origin = ifelse(Goods.description != '' | Origin.country != '',
                             'ICS', 'Shipping company'),
         data_origin = ifelse(Goods.description != '' & 
                                Shipping.company.goods.description != '',
                             'Both', data_origin)) %>% 
  select(container_id, data_origin, Arrival.date, 
         Loading.country, Origin.country,
         Goods.description, Shipping.company.goods.description,
         Goods.khapra.classification, Goods.hitchhiker.risk.classification
         ) %>% 
  mutate_at(.vars = vars(Goods.khapra.classification,
                         Goods.hitchhiker.risk.classification),
            .funs = list(~ ifelse(.=='Unknown', '', .))) %>% 
  filter(Arrival.date != '') %>%
  setNames(tolower(gsub("\\.","_",names(.)))) %>%
  unique() %>% 
  group_by(container_id, arrival_date,data_origin,
           loading_country, origin_country, goods_description,
           goods_khapra_classification, goods_hitchhiker_risk_classification) %>% 
  summarise(company_goods_description = paste(shipping_company_goods_description,
                          collapse = ' ;; ')) %>%
  rename(ICS_goods_description = goods_description,
         khapra_risk = goods_khapra_classification,
         hitchhiker_risk = goods_hitchhiker_risk_classification) %>% 
  ungroup() %>% 
  arrange(container_id, desc(arrival_date)) %>% 
  mutate(sampled = ifelse(duplicated(container_id), F, T)) %>%
  relocate(sampled, container_id, data_origin) %>% 
  relocate(khapra_risk, hitchhiker_risk, 
           .after = company_goods_description)

history %>% head
history$container_id %>% unique %>% length
table(history$data_origin)


## visual surveys ----------------------------------------------------------

ship_simple %>% names
# 1 4 7 31:53
table(ship_simple$Identification.type,
      ship_simple$Lab.name)
ship_simple %>% filter(Identification.type == '') %>% head

visual <- ship_simple[,c(1,4, 7, 31:53)] %>% 
  filter(Identification.type != '') %>% 
  setNames(str_trim(tolower(gsub("[[:punct:]]+"," ",names(.))))) %>% 
  setNames(gsub(' ', '_', names(.))) %>% 
  setNames(gsub('_animal', '', names(.))) %>% 
  setNames(gsub('_seed', '', names(.))) %>% 
  setNames(gsub('_contaminant', '', names(.))) %>%
  unique() %>% 
  split(., .$identification_type)

visual$SEED %>% names
colSelection <- list(ANIMAL = c(1,6:11,18:20, 5),
                     CONTAMINATION = c(1, 25:26, 5),
                     SEED = c(1, 6:11, 22:24,5))

visual_simple <- lapply(1:3, 
                        function(x) visual[[x]][colSelection[[x]]])
names(visual_simple) <- names(visual)

sapply(visual_simple, ncol)
sapply(visual_simple, nrow)
lapply(visual_simple, names)


# metabarcoding -----------------------------------------------------------


dna_simple <- dna %>% filter(PI. >= 90) %>%
  setNames(gsub('\\.', '_', names(.))) %>%
  setNames(gsub('QC_', 'QC', names(.))) %>%
  setNames(gsub('PI_', 'PI', names(.))) %>%
  rename_if(!grepl('X2', names(.)), tolower) %>% 
  rename(PI = pi, QC = qc, asv_id = id) %>% 
  mutate(high_priority = best_hit %in% str_trim(priority$species)) %>% 
  relocate(high_priority, .after = read_count)


rna_simple <- rna %>% filter(PI. >= 90) %>% 
  setNames(gsub('\\.', '_', names(.))) %>%
  setNames(gsub('_cDNA', '', names(.))) %>%
  setNames(gsub('QC_', 'QC', names(.))) %>%
  setNames(gsub('PI_', 'PI', names(.))) %>%
  rename_if(!grepl('X2', names(.)), tolower) %>% 
  rename(PI = pi, QC = qc, asv_id = id) %>% 
  mutate(high_priority = best_hit %in% str_trim(priority$species)) %>% 
  relocate(high_priority, .after = read_count)


# species specific --------------------------------------------------------
lapply(spp, names)

for (spp_name in 1:length(spp)) {
  spp[[spp_name]]$common_name <- names(spp)[spp_name]
}
spp_renamed <- spp
spp_renamed$Khapra_beetle <- spp$Khapra_beetle %>% 
  setNames(ifelse(grepl('NA_', names(.)),
                  paste0(names(.), '.Ct'), names(.))) %>% 
  setNames(gsub('Furui_', 'e', names(.))) %>% 
  setNames(gsub('_R', '.rep.', names(.))) %>% 
  select_if(!grepl('Olson', names(.)))

names(spp)

species_name <- data.frame(common_name = names(spp_renamed),
                           species =  c('Lymantria dispar asiatica',
                                             'Halyomorpha halys',
                                             'Wasmannia auropunctata',
                                             'Trogoderma granarium',
                                             'Lycorma delicatula'))

spp_specific <- do.call('bind_rows', spp_renamed) %>% 
  rename(container_id = Container.Number,
         lab_code = LAB.Code.for.Container,
         sample_id = Unique.Sample.Identifier,
         method_collection = Method.Collection) %>% 
  left_join(species_name) %>% 
  mutate(common_name = ifelse(common_name == 'BMSB',
                              'Brown marmorated stink bug',
                              common_name),
         common_name = gsub('_', ' ', tolower(common_name)),
         sample_id = sub('C05246_', '', sample_id),
         sample_id = sub('\\.1','', sample_id)) %>%
  filter(!grepl('\\.2', sample_id)) %>% 
  select(-method_collection, -container_id) %>% 
  relocate(sample_id, lab_code, common_name, species) %>% 
  unique()

spp_specific %>% names
spp_specific[,c(1,3)] %>% duplicated %>% table


# save --------------------------------------------------------------------

# project code: C05246

write.csv(collection, file ='./output/C05246_collection.csv', row.names = F)
write.csv(history, file ='./output/C05246_container_history.csv', row.names = F)

write.csv(visual_simple$ANIMAL,
          file ='./output/C05246_visual_animal.csv', row.names = F)
write.csv(visual_simple$CONTAMINATION,
          file ='./output/C05246_visual_contamination.csv', row.names = F)
write.csv(visual_simple$SEED,
          file ='./output/C05246_visual_seed.csv', row.names = F)

write.csv(spp_specific,
          file ='./output/C05246_genetic_species_specific.csv', row.names = F)
write.csv(dna_simple,
          file ='./output/C05246_genetic_metabarcoding_dna.csv', row.names = F)
write.csv(rna_simple,
          file ='./output/C05246_genetic_metabarcoding_rna.csv', row.names = F)
