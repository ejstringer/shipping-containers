
## This script it to create three datasets. First is the presence of DNA or RNA 
## and the congruence across two types of genetic data: species specific and 
## metabarcoding. An additional interest is whether containers contain multiple
## pest species or if they are dispersed across multitudes of containers, with 
## containers generally only having one or two pests detected.
##    The next dataset is for container variables of interest, pretty simple, 
## things like age of container, goods being transported, and risk country. 
##    The next dataset is species specific, filtered down to only containers 
## where DNA was detected. Here we are interested in how DNA cq relates to 
## detection of RNA but also how sample quality and yeild may effect DNA cq.



# load --------------------------------------------------------------------

source('code/libraries.R')


collection <- read.csv('./output/C05246_collection.csv')

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

spp <- collection %>% 
  left_join(specific) %>% 
  relocate(container_id) %>% 
  select(-visual_lab) %>% 
 # filter(sample_method != 'VAC UNDERFLOOR') %>% 
 # filter(!sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  mutate(across(eRNA_rep_1_Ct,
                ~ifelse((common_name == 'brown marmorated stink bug' &
                            eRNA_rep_1_Ct > 0 & eRNA_rep_2_Ct == 0 &
                            eDNA_rep_1_Ct == 0), 0, .x))) %>% 
  mutate(across(eRNA_rep_1_Ct:eRNA_rep_3_Ct,
                ~ifelse((common_name == 'electric ant' &
                           eRNA_rep_1_Ct > 0 & 
                           eDNA_rep_1_Ct == 0), 0, .x))) 
  spp %>% 
    filter(is.na(eDNA_rep_1_Ct)) %>% 
  mutate(across(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
                ~ifelse(is.na(.x),34.56, .x))) # based on mean of 3/3 Furui

rna <- read.csv('./output/C05246_genetic_metabarcoding_rna.csv')%>% 
  filter(PI > 0.99)
dna <- read.csv('./output/C05246_genetic_metabarcoding_dna.csv') %>% 
  filter(PI > 0.99)

animal <- read.csv('./output/C05246_visual_animal.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method)

container_history <- read.csv('./output/C05246_container_history.csv') 
container_trial <- read.csv('./data/shipping_meta/Approach rate trial container data summary.csv') 


spp5 <- data.frame(species = unique(spp$species)) %>% 
  separate(species,sep = ' ', into = c('genus'), remove = F) 

###########################################################################-
#           congruence                                                  ----
###########################################################################-

## metabarcoding -----------------------------------------------------------

rna_detected <- rna %>% filter(high_priority) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = X20210506_0046:X20210811_2063, names_to = 'sample_id',
               values_to = 'reads_rna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_rna_spp = best_hit %in% spp5$species) %>% 
  rename(species = best_hit) %>% 
  group_by(container_id) %>% 
  summarise(pests_metaRNA = sum(reads_rna > 0),
            pests_meta5RNA = sum(detected_rna_spp &reads_rna > 0))

rna_detected[,c('container_id')] %>% duplicated %>% table

names(dna) %>% tail
dna_detected <- dna %>% filter(high_priority) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>%
  pivot_longer(cols = X20210504_0002:X20210813_2107, names_to = 'sample_id',
               values_to = 'reads_dna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_dna_spp = best_hit %in% unique(spp$species)) %>% 
  rename(species = best_hit) %>% 
  unique() %>% 
  group_by(container_id) %>% 
  summarise(pests_metaDNA = sum(reads_dna > 0),
            pests_meta5DNA = sum(detected_dna_spp &reads_dna > 0)) %>% 
  left_join(rna_detected)

dna_detected[,c('container_id')] %>% duplicated %>% table


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
                             colour = factor(pests_sppRNA)))+
    geom_jitter(aes(size = factor(pests_sppRNA)), alpha = 0.5)+
    theme_classic()+
    scale_x_continuous(breaks = 0:4)

pest_detection %>% 
  ggplot(aes(factor(pests_sppDNA), pests_meta5DNA, 
             colour = factor(pests_sppRNA)))+
  geom_jitter(aes(size = factor(pests_sppRNA)), alpha = 0.5)+
  theme_classic()+
  scale_y_continuous(limits = c(0,3))
  
### save --------------------------------------------------------------------

write.csv(pest_detection, './output/pest_detection_methods.csv', row.names = F)
  
    
###########################################################################-
#           container variables                                         ----
###########################################################################-


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
           Container.grade, 
           Age.at.earliest.sample.date) %>%
    rename(age = Age.at.earliest.sample.date, 
           container_id = Container.number) %>% 
    mutate(Container.grade = sub(' \\(steel floor\\)', '', Container.grade),
           Container.grade = sub('General purpose\\/', '', Container.grade),
           Container.grade = ifelse(Container.grade %in% c('Cotton quality', 
                                                           'Flexi tank',
                                                           'Unknown'),
                                    'General purpose', Container.grade)) %>% 
    setNames(tolower(gsub('\\.', '_', names(.)))) %>% 
    left_join(goods_sum) %>% 
    mutate_if(is.character, as.factor)
  container %>% summary
container$container_grade %>% table  

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
  summarise(max_cq = min(cq)) %>%
  mutate(max_cq = ifelse(max_cq == 99, 0, max_cq)) %>% 
  filter(max_cq > 0) %>% 
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
  #filter(sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  select(container_id, species, common_name, 
         sample_weight, conc_2_ng_ul, x260_230, x260_280, 
         eDNA_rep_1_Ct:eDNA_rep_3_Ct) %>% 
  pivot_longer(eDNA_rep_1_Ct:eDNA_rep_3_Ct,
               names_to = 'rep', values_to = 'cq') %>% 
  mutate(rep = sub('eDNA_rep_', '', rep),
         rep = sub('_Ct','', rep)) %>% 
  left_join(pestrna) %>% 
  mutate(pest_rna = ifelse(pest_rna, 1, 0),
         cq_rank = case_when(
           cq== 0 ~ 'zero',
           cq < 35 ~ 'high',
           cq < 40 ~ 'med',
           cq >= 40 ~ 'low')) %>% 
  left_join(maxrepcq) # min cq not zero



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


# models -------------------------------------------------------------------

## DNA cq -------
dna <- DNA_cq_sample %>% 
  select(container_id,pest_rna, min_cq, common_name,
         conc_2_ng_ul, x260_230, x260_280) %>% 
  unique() %>% 
  rename(cq = min_cq) %>% 
  left_join(container)
dna %>% head
lm(cq ~ x260_230, 
         data = dna) %>% plot

boxplot((dna$cq))
ggplot(dna, aes(x260_230, cq))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()
cor(dna$cq, dna$x260_230)
m <- glm(pest_rna ~ cq, 
      data = dna, family = binomial) 

summary(m)
x <- seq(20,50, 0.1)
ndata <- data.frame(cq = x)
y <- predict(m, newdata=data.frame(cq = x),
             type="response")
y <- predict(m, newdata=data.frame(cq = x),se.fit = T, type = 'response')
ndata 
ggplot(cbind(ndata, y), aes(cq, y))+
  geom_ribbon(aes(ymin = right_lwr, ymax = right_upr),
              fill = 'grey', alpha = 0.5)+
  geom_line(colour = 'pink',lwd = 2)+
  theme_classic()+
  xlab('rep(1-3) minimum cq')+
  ylab('RNA detected in container')


## DNA 01 -----

pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0))

pest %>% data.frame %>% head
m <- glm(dna ~ container_grade, data = pest, family = binomial) 
summary(m)


fam <- family(m)
fam
str(fam)

ilink <- fam$linkinv
ilink

## grad the inverse link function
ilink <- family(m)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(ndata, setNames(as_tibble(predict(m, ndata, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))
## show
ndata


x <- unique(pest$container_grade)#seq(20,50, 0.1)

ndata <- data.frame(container_grade = x)
y <- predict(m, newdata=data.frame(container_grade = x),
             type="response", se.fit =T)
#y <- predict(m, newdata=data.frame(cq = x),se.fit = T, type = 'response')

ggplot(ndata, aes(container_grade, fit_link))+
 # geom_point(data = dna,alpha = 0.5)+
  geom_errorbar(aes(ymin = fit_link - right_lwr, ymax =fit_link + right_upr),
              fill = 'grey', alpha = 0.5, width = 0, lwd = 1)+
  geom_point(colour = 'black',size = 4)+
  theme_classic()+
  ylab('probability of detecting DNA (se)')
