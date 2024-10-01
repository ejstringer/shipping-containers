
## This script it to create three datasets. First is the presence of DNA or RNA 
## and the congruence across two types of genetic data: species specific and 
## metabarcoding. An additional interest is whether containers contain multiple
## pest species or if they are dispersed across multitudes of containers, with 
## containers generally only having one or two pests detected.
##    The next dataset is for container variables of interest, pretty simple, 
## things like age of container, goods being transported, and risk country. 
##    The next dataset is species specific, filtered down to only containers 
## where DNA was detected. Here we are interested in how DNA cq relates to 
## detection of RNA but also how sample quality and yield may effect DNA cq.



# load --------------------------------------------------------------------

source('code/libraries.R')


collection <- read.csv('./output/C05246_collection.csv')

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

spp <- collection %>% 
  left_join(specific) %>% 
  relocate(container_id) %>% 
  select(-visual_lab) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  filter(!sample_id %in% c('X20210721_1510', 'X20210511_0088'))


  spp %>% 
    filter(is.na(eDNA_rep_1_Ct)) %>% 
  mutate(across(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
                ~ifelse(is.na(.x),34.56, .x))) %>%  # based on mean of 3/3 Furui
    data.frame # DNA NA is actually present

rna <- read.csv('./output/C05246_genetic_metabarcoding_rna.csv')%>% 
  filter(PI > 99)
dna <- read.csv('./output/C05246_genetic_metabarcoding_dna.csv') %>% 
  filter(PI > 99)

seed <- read.csv('./output/C05246_visual_seed.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method)


contam <- read.csv('./output/C05246_visual_contamination.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method)

container_history <- read.csv('./output/C05246_container_history.csv') 
container_trial <- read.csv('./data/shipping_meta/Approach rate trial container data summary.csv') 


spp5 <- data.frame(species = unique(spp$species)) %>% 
  separate(species,sep = ' ', into = c('genus'), remove = F)
priority <- read.csv('./data/high_priority_insect_species.csv') %>% 
  mutate(species = str_trim(species))

###########################################################################-
#           congruence                                                  ----
###########################################################################-



## metabarcoding -----------------------------------------------------------

rna_detected <- rna %>% #filter(high_priority) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = X20210506_0046:X20210811_2063, names_to = 'sample_id',
               values_to = 'reads_rna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(reads_rna = ifelse(reads_rna < 100, 0, reads_rna)) %>% 
  group_by(container_id, best_hit) %>% 
  summarise(reads_rna = sum(reads_rna)) %>% 
  mutate(detected_rna_spp = best_hit %in% spp5$species,
         reads_rna_priority =  ifelse(best_hit %in% priority$species,
                                      reads_rna, 0),
         detected_rna_priority =  best_hit %in% priority$species) %>% 
  rename(species = best_hit) %>% 
  group_by(container_id) %>% 
  summarise(total_reads_RNA = sum(reads_rna),
            pest_reads_RNA = sum(reads_rna_priority),
            richness_RNA = sum(reads_rna > 0),
            pests_metaRNA = sum(detected_rna_priority & reads_rna > 0),
            pests_meta5RNA = sum(detected_rna_spp & reads_rna > 0)) %>% 
  mutate(rel_abundance_RNA = pest_reads_RNA/total_reads_RNA)

rna_detected[,c('container_id')] %>% duplicated %>% table

names(dna) %>% tail
dna_detected <- dna %>% #filter(high_priority) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>%
  pivot_longer(cols = X20210504_0002:X20210813_2107, names_to = 'sample_id',
               values_to = 'reads_dna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(reads_rna = ifelse(reads_dna < 100, 0, reads_dna)) %>% 
  group_by(container_id, best_hit) %>% 
  summarise(reads_dna = sum(reads_dna)) %>% 
  mutate(detected_dna_spp = best_hit %in% spp5$species,
         reads_dna_priority =  ifelse(best_hit %in% priority$species,
                                      reads_dna, 0),
         detected_dna_priority =  best_hit %in% priority$species) %>% 
  rename(species = best_hit) %>% 
  group_by(container_id) %>% 
  summarise(total_reads_DNA = sum(reads_dna),
            simpsondna = simpson(reads_dna),
            shannondna = shannon(reads_dna),
            pest_reads_DNA = sum(reads_dna_priority),
            richness_DNA = sum(reads_dna > 0),
            pests_metaDNA = sum(detected_dna_priority & reads_dna > 0),
            pests_meta5DNA = sum(detected_dna_spp & reads_dna > 0)) %>% 
  mutate(rel_abundance_DNA = pest_reads_DNA/total_reads_DNA) %>% 
  left_join(rna_detected)

dna_detected[,c('container_id')] %>% duplicated %>% table
dna_detected$pests_meta5DNA %>% sum

dna_detected$rel_abundance_DNA[dna_detected$pests_metaDNA>0] %>% boxplot

colSums(dna_detected[,2:5], na.rm = T)
colMeans(dna_detected[,2:5], na.rm = T)


## species specific --------------------------------------------------------


specific_detection <- spp %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
 # filter(!sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  mutate(across(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
                ~ifelse(is.na(.x),34.56, .x))) %>%  # 3/3 mean Furui
  mutate(specific_DNA = select(., eDNA_rep_1_Ct:eDNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         specific_RNA = select(., eRNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         detected_specDNA = specific_DNA > 0,
         detected_specRNA = specific_RNA > 0) %>% 
  group_by(container_id) %>% 
  summarise(pests_sppDNA = sum(detected_specDNA),
            pests_sppRNA = sum(detected_specRNA),) %>% 
  arrange(container_id) 

  

## detections --------------------------------------------------------------
pest_detection <- specific_detection %>% 
    left_join(dna_detected)

### plots -------------------------------------------------------------------



  pest_detection %>% 
  ggplot(aes(pests_sppDNA, pests_metaDNA, 
                             colour = factor(pests_meta5DNA)))+
    geom_jitter(size = 2.5,alpha = 0.7, width = 0.25, height = 0.25)+
    theme_bw()


pest_detection %>% 
  filter(complete.cases(pests_metaRNA)) %>% 
  ggplot(aes(pests_sppDNA, pests_sppRNA, 
             colour = factor(pests_metaRNA)))+
  geom_jitter(size = 2.5,alpha = 0.6, width = 0.1, height = 0.1)+
  theme_bw()
  

pest_detection %>% 
  ggplot(aes(factor(pests_sppDNA), pests_meta5DNA, 
             colour = factor(pests_sppRNA)))+
  geom_jitter(aes(size = factor(pests_sppRNA)), alpha = 0.5)+
  theme_classic()+
  scale_y_continuous(limits = c(0,3))


t.test(pest_detection$richness_DNA,
       pest_detection$richness_RNA,
       paired = T)

pest_detection %>% 
  filter(complete.cases(richness_RNA)) %>% 
 # lm(richness_RNA ~ richness_DNA, data = .) %>% summary
  ggplot(aes(richness_DNA, richness_RNA))+
  geom_point()+
  theme_classic()+
  geom_smooth(method = 'lm')+
  geom_abline(slope = 1, intercept = 0, colour = 'red', lty = 2)
  
### save --------------------------------------------------------------------

write.csv(pest_detection, './output/pest_detection_methods.csv', row.names = F)
  
    
###########################################################################-
#           container variables                                         ----
###########################################################################-

## seeds -------------------------------------------------------------------

seedcat<- seed %>% 
  select(container_id, seed_category) %>%
  unique() %>% 
  group_by(container_id) %>% 
  summarise(category = paste(seed_category, collapse = '-')) %>% 
  mutate(category = sub('-Other', '', category),
         category = sub('Other-', '', category),
         category = ifelse(category == 'Other', 'other seed', 'pulse/cereal')) %>% 
  mutate_all(as.factor) %>% 
  rename(seed = category)
seedcat$container_id %>% duplicated %>% table
seedcat %>% summary

## contamination ----------
contam %>% head
contam %>% mutate_all(factor) %>% summary


contam_cat <- contam %>% 
  select(container_id, matrix_type) %>%
  unique() %>% 
  mutate(matrix_type= str_replace_all(matrix_type, "[^[:alnum:]]", " "),
         wood = ifelse(grepl('ood', matrix_type)|grepl('imber', matrix_type), 
                           1, 0),
         soil = ifelse(grepl('oils', matrix_type), 1, 0)) %>%
  select(-matrix_type) %>% 
  #unique %>% 
  group_by(container_id) %>% 
  summarise(wood = sum(wood),
            soil = sum(soil)) %>% 
  mutate(wood = ifelse(wood > 0, 1, 0),
         soil = ifelse(soil > 0, 1, 0))
contam_cat  %>%  summary 
## goods ICS -----------
container_history$container_id %>% unique %>% length
ics <-container_history %>% filter(data_origin != 'Shipping company')

icssampled <- ics %>% arrange(container_id, arrival_date) %>% 
  mutate(ics_sampled = ifelse(duplicated(container_id), F, T)) %>%
  filter( ymd(arrival_date) > ymd('2020-01-01')) %>% 
  relocate(ics_sampled) %>% 
  filter(ics_sampled) %>% 
  # filter(hitchhiker_risk != '') %>% 
  mutate(ICS_goods_description = tolower(ICS_goods_description),
         goods_risk = ifelse(grepl('food', hitchhiker_risk), 'food', 'other'),
         goods_risk = ifelse(grepl('animal', hitchhiker_risk), 'food', goods_risk),
         goods_risk = ifelse(grepl('kb', hitchhiker_risk), 'food', goods_risk),
         #  goods_risk = ifelse(grepl('soil', hitchhiker_risk), 'soil', goods_risk),
         goods_risk = ifelse(grepl('plant', hitchhiker_risk), 'plant', goods_risk),
         goods_risk = ifelse(grepl('unknown', hitchhiker_risk), 'unknown', goods_risk),
         goods_risk = ifelse(hitchhiker_risk == '', 'unknown', goods_risk)) %>% 
  mutate(goods_risk = ifelse(goods_risk == 'unknown' & grepl('food',ICS_goods_description),'food', goods_risk),
         goods_risk = ifelse(goods_risk == 'unknown' & grepl('beans',ICS_goods_description),'food', goods_risk),
         goods_risk = ifelse(goods_risk == 'unknown' & grepl('wood',ICS_goods_description) & ICS_goods_description != 'sherwood',
                             'plant', goods_risk),
         goods_risk = ifelse(goods_risk == 'unknown' & grepl('paper',ICS_goods_description),
                             'plant', goods_risk)) %>% 
  select(container_id, goods_risk)#glimpse
icssampled$goods_risk %>% table
icssampled$goods_risk %>% table %>% sum -382
icssampled$container_id %>% duplicated %>% table


## goods and risk country --------------------------------------------------

  goods_sum <- container_history %>%
    filter(sampled) %>% 
    mutate(timber = grepl('timber', tolower(company_goods_description)) |
             grepl('wood', tolower(company_goods_description)) |
             grepl('wood', tolower(ICS_goods_description)),
           food = grepl('seed', tolower(company_goods_description)) |
             grepl('wheat', tolower(company_goods_description))|
             grepl('seed', tolower(ICS_goods_description))|
             grepl('food', tolower(company_goods_description))|
             grepl('vegetables', tolower(company_goods_description))|
             grepl('tomato', tolower(company_goods_description))|
             grepl('cereal', tolower(company_goods_description))|
             grepl('animal', tolower(company_goods_description))|
             grepl('potato', tolower(company_goods_description))|
             grepl('nuts', tolower(company_goods_description))|
             grepl('palm', tolower(company_goods_description)),
           risk_country = ifelse(loading_country %in% c('Sri Lanka', 'Saudi Arabia', 
                                                        'Qatar', 'Pakistan', 'China', 'Israel',
                                                        'India', 'Greece','Egypt',
                                                        'Bangladesh', 'Bahrain'),
                                 'yes','no'),
           goods = case_when(
             food  ~ 'food',
             timber ~ 'timber',
             !(food | timber) ~ 'other_goods'
           )) %>% 
    select(container_id,risk_country, goods)
  goods_sum


## container specs ---------------------------------------------------------

  container <- container_trial %>% 
    select(Container.number,
           Container.grade,  Container.size,
           Age.at.earliest.sample.date) %>%
    rename(age = Age.at.earliest.sample.date, 
           container_id = Container.number) %>% 
    mutate(Container.grade = sub(' \\(steel floor\\)', '', Container.grade),
           Container.grade = sub('General purpose\\/', '', Container.grade),
           Container.grade = ifelse(Container.grade %in% c('Cotton quality', 
                                                           'Flexi tank',
                                                           'Unknown'),
                                    'General purpose', Container.grade),
           Container.size = factor(Container.size)) %>% 
    setNames(tolower(gsub('\\.', '_', names(.)))) %>% 
    left_join(goods_sum) %>% 
    left_join(icssampled) %>% 
    left_join(seedcat) %>% 
    left_join(contam_cat) %>% 
    mutate(goods_risk = ifelse(is.na(goods_risk), 'unknown', goods_risk),
           seed = ifelse(is.na(seed), 'no seed', as.character(seed)),
           soil = ifelse(is.na(soil), 0, soil),
           wood = ifelse(is.na(wood), 0, wood)) %>% 
    mutate_if(is.character, as.factor)
  container %>% summary
container$container_grade %>% table  

container$wood %>% table  

## save --------------------------------------------------------------------

write.csv(container, './output/pest_container_specs.csv', row.names = F)

###########################################################################-
#           DNA rep cq                                                  ----
###########################################################################-


## rna in container --------------------------------------------------------


pestrna <- spp %>% 
  mutate(pest = select(., eRNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         pest_rna = pest > 0) %>% 
  select(container_id, common_name, pest_rna)


## max rep cq --------------------------------------------------------------

maxrepcq <- spp %>% 
  mutate(across(eDNA_rep_1_Ct:eRNA_rep_3_Ct, ~ifelse(.x == 0, 99, .x))) %>%
  pivot_longer(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
               names_to = 'rep', values_to = 'cq') %>% 
  group_by(container_id, common_name) %>%
  filter(cq != 99) %>% 
  summarise(max_cq = min(cq),
            pos_dna = n()) %>%
  rename(min_cq = max_cq)

## dna reps ----------------------------------------------------------------
spp %>% names
spp[c(1,grep('eRNA', names(spp)))]
DNA_cq_sample <- spp %>%  
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  mutate(detected_DNA = select(., eDNA_rep_1_Ct:eDNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         detected_DNA = detected_DNA > 0) %>% 
  filter(detected_DNA) %>%
  mutate(x260_280 =  as.numeric(sub(',', '', x260_280)),
         x260_230 = as.numeric(x260_230)) %>% 
  select(container_id, species, common_name, 
         sample_weight, conc_2_ng_ul, x260_230, x260_280, 
         eDNA_rep_1_Ct:eDNA_rep_3_Ct) %>% 
  pivot_longer(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
               names_to = 'rep', values_to = 'cq') %>% 
  mutate(rep = sub('eDNA_rep_', '', rep),
         rep = sub('_Ct','', rep)) %>% 
  unique() %>% 
  left_join(pestrna) %>% 
  mutate(pest_rna = ifelse(pest_rna, 1, 0),
         cq_rank = case_when(
           cq== 0 ~ 'zero',
           cq < 35 ~ 'high',
           cq < 40 ~ 'med',
           cq >= 40 ~ 'low')) %>% 
  left_join(maxrepcq) # min cq not zero

DNA_cq_sample[16,]
DNA_cq_sample %>% filter(container_id == 'TGHU0965042')
## save --------------------------------------------------------------------

write.csv(DNA_cq_sample, './output/pest_dna_cq_rna.csv', row.names = F)

## plots -------------------------------------------------------------------

DNA_cq_sample %>% 
 # filter(cq > 0) %>% 
  ggplot(aes(cq, pest_rna))+
  geom_point()+
  theme_classic()+
  geom_smooth()

DNA_cq_sample %>% 
  filter(cq > 0) %>% 
  ggplot(aes(sample_weight, cq))+
  scale_x_log10()+
  geom_point()+
  theme_classic()+
  geom_smooth()


DNA_cq_sample %>% 
  filter(cq > 0) %>% 
  ggplot(aes(x260_230, cq))+
  scale_x_log10()+
  geom_point()+
  theme_classic()+
  geom_smooth()

DNA_cq_sample %>% 
  filter(cq > 0) %>% 
  ggplot(aes(x260_280, cq))+
  geom_point()+
  theme_classic()+
  geom_smooth()

names(DNA_cq_sample)

## ---------
