
# congruence

# load --------------------------------------------------------------------

source('./code/libraries.R')

collection <- read.csv('./output/C05246_collection.csv')

spp <- read.csv('./output/C05246_genetic_species_specific.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR')

collection <- collection %>% filter(sample_method != 'VAC UNDERFLOOR')

animal <- read.csv('./output/C05246_visual_animal.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method)

rna <- read.csv('./output/C05246_genetic_metabarcoding_rna.csv')
dna <- read.csv('./output/C05246_genetic_metabarcoding_dna.csv')


# visual detections -------------------------------------------------------

vis_detected <- animal %>% 
  #filter(species %in% unique(spp$species)) %>% 
  mutate(species = ifelse(species %in% unique(spp$species), species,
                          'other'),
         specific = species %in% unique(spp$species)) %>% 
  group_by(container_id, species) %>% 
  summarise(specimens_visual = n(),
            detected_visual = sum(specific) > 0)
vis_detected %>% arrange(desc(species))


# metabarcoding detections ------------------------------------------------

rna_detected <- rna %>% filter(best_hit %in% unique(spp$species)) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = X20210506_0046:X20210811_2063, names_to = 'sample_id',
               values_to = 'reads_rna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_rna = reads_rna > 0) %>% 
  rename(species = best_hit)

rna_detected[,c('container_id', 'species')] %>% duplicated %>% table

names(dna) %>% tail
dna_detected <- dna %>% filter(best_hit %in% unique(spp$species)) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>%
  pivot_longer(cols = X20210504_0002:X20210813_2107, names_to = 'sample_id',
               values_to = 'reads_dna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_dna = reads_dna > 0) %>% 
  rename(species = best_hit) %>% 
  unique()

dna_detected[,c('container_id', 'species')] %>% duplicated %>% table

# overview ----------------------------------------------------------------
spp %>% names
spp %>% head

table(duplicated(collection$sample_id))
no.cont <- (table(collection$container_id))
table(no.cont)
con3 <- no.cont[no.cont>1]

collection[collection$container_id %in% names(con3),]
spp[spp$container_id %in% names(con3),]
spp %>% glimpse()

spp_detected <- spp %>% 
  select(-lab_code) %>% 
  group_by(container_id, common_name, species) %>% 
  summarise_if(is.double, mean) %>% 
  ungroup() %>% 
  mutate(detected = select(., eDNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         detected_specific = detected > 0) %>% 
  left_join(vis_detected) %>% 
  arrange(desc(detected_visual),
          desc(specimens_visual)) %>%
  left_join(rna_detected) %>%
  left_join(dna_detected) %>% 
  select_if(!grepl('rep', names(.))) %>% 
  glimpse
 spp_detected$container_id %>% table %>% table
 
tapply(spp_detected$detected_specific, 
       spp_detected$common_name, sum)

# detections
detection <- spp_detected %>% 
  group_by(common_name, species) %>% 
  summarise(specific = sum(detected_specific),
            dna = sum(detected_dna, na.rm = T),
            rna = sum(detected_rna, na.rm = T),
             visual = sum(detected_visual, na.rm = T),
            containers = n()) %>% 
  arrange(desc(specific)) %>% 
  pivot_longer(cols = dna:visual, names_to = 'data', values_to = 'detections')

# congruence
congruence <- spp_detected %>% 
  rowwise() %>% 
  mutate(congruence_dna = (detected_dna == T & detected_specific == T),
         congruence_rna = (detected_rna == T & detected_specific == T),
         congruence_vis = (detected_visual == T & detected_specific == T)) %>% 
  group_by(common_name, species) %>% 
  summarise(specific = sum(detected_specific),
            dna = sum(congruence_dna, na.rm = T),
            rna = sum(congruence_rna, na.rm = T),
            visual = sum(congruence_vis, na.rm = T),
            containers = n()
            ) %>% 
  arrange(desc(specific)) %>% 
  pivot_longer(cols = dna:visual, names_to = 'data', values_to = 'congruence')

detection %>% left_join(congruence) %>% 
  select(-containers, -species) %>% 
  mutate(prop_detection = detections/specific,
         prop_congruence = congruence/detections) %>% 
  arrange(data) %>% 
  relocate(common_name, data, prop_congruence)

spp_detected %>% filter(container_id %in% vis_detected$container_id) %>% arrange(common_name)
