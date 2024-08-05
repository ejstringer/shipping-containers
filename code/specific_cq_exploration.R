
# congruence

# load --------------------------------------------------------------------

source('./code/libraries.R')

arrival <- read.csv('./output/C05246_container_history.csv') %>% 
  filter(sampled) %>% select(container_id, arrival_date) 

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

spp <- read.csv('./output/C05246_collection.csv') %>% 
  left_join(arrival) %>% 
  left_join(specific) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  mutate(days_since = ymd(collection_date)- ymd(arrival_date)) %>% 
  relocate(container_id) %>% 
  relocate(days_since, .after = collection_date) %>% 
  pivot_longer(cols = eDNA_rep_1_Ct:eRNA_rep_3_Ct, names_to = 'molecule',
               values_to = 'cq') %>% 
  mutate(molecule = sub('_rep_', '_', molecule),
         molecule = sub('_Ct', '', molecule)) %>% 
  separate(molecule, into = c('molecule', 'rep'),sep = '_')
spp %>% names
spp %>% glimpse()

sppDNA <- filter(spp, #molecule == 'eDNA',
                 cq > 0, 
                 days_since > 0) %>% 
  mutate(seed = ifelse(grepl('SEED', visual_contents),
                       'Seed', 'No seed'))
sppDNA %>% names
sppDNA %>% nrow
sppDNA %>% 
  ggplot(aes(days_since, cq, colour = common_name))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  facet_grid(~molecule)
lmer(cq ~ days_since + molecule +(1|sample_id),data = sppDNA) %>% 
  summary

sppDNA %>% 
  ggplot(aes(common_name, cq, fill = molecule))+
  geom_boxplot()+
  theme_bw()

sppDNA %>% 
  ggplot(aes(seed, cq, fill = seed))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c('lightblue', 'lightgreen'))+
  facet_grid(molecule~common_name)

boxplot(sppDNA$cq)

