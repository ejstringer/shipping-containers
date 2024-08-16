
source('./code/libraries.R')


# load -----------------------------------------------------------------

container_history <- read.csv('./output/C05246_container_history.csv') 
container_trial <- read.csv('./data/shipping_meta/Approach rate trial container data summary.csv') 

collection <- read.csv('./output/C05246_collection.csv')

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

days_since_arrival <-  container_history %>% 
  filter(sampled) %>% select(container_id, arrival_date) %>% 
  left_join(collection) %>% 
  mutate(days_since = as.numeric(ymd(collection_date)- ymd(arrival_date)))


## adjust arrival times ----------------------------------------------------

daysminus <- days_since_arrival%>% 
  group_by(container_id) %>% 
  summarise(days_since = mean(days_since)) %>% 
  filter(days_since<0) 

adjustedDaySince  <- collection[,c('sample_id', 'container_id', 'collection_date')] %>%
  filter(!sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  select(-sample_id) %>% 
  unique() %>% 
  filter(container_id %in% daysminus$container_id) %>% 
  left_join(container_history) %>% 
  relocate(collection_date, .after = arrival_date) %>%
  filter(!sampled) %>% 
  arrange(container_id, desc(arrival_date)) %>% 
  mutate(sampled = ifelse(duplicated(container_id), F, T),
         days_since2 = as.numeric(ymd(collection_date)- ymd(arrival_date))) %>% 
  relocate(days_since2, .after = arrival_date) %>% 
  filter(sampled) %>% 
  select(container_id, arrival_date, loading_country) %>% 
  rename(new_arrival_date= arrival_date)

adjustedDaySince
adjust <- container_history$container_id %in% adjustedDaySince$container_id
container_history$sampled[adjust] <- FALSE

history <- container_history %>%
  left_join(adjustedDaySince) %>% 
  mutate(sampled = ifelse(!is.na(new_arrival_date) & arrival_date == new_arrival_date, 
                           TRUE, sampled)) %>%
  group_by(container_id) %>% 
  filter(sampled) %>% arrange

arrival <- history %>% 
  filter(sampled) %>% select(container_id, arrival_date) 

## spp ---------------

spp <- collection %>% 
  left_join(arrival) %>% 
  left_join(specific) %>% 
  filter(sample_method != 'VAC UNDERFLOOR',
      #   common_name == speciesNAME[4],
         !sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  mutate(days_since = as.numeric(ymd(collection_date)- ymd(arrival_date))) %>%
  mutate(across(eDNA_rep_1_Ct:eRNA_rep_3_Ct, ~ifelse(.x == 0, 99, .x))) %>% 
  relocate(container_id) %>% 
  relocate(days_since, .after = collection_date) %>% 
  rowwise() %>% 
  mutate(eDNA = paste(sort(c(eDNA_rep_1_Ct,
                             eDNA_rep_2_Ct,
                             eDNA_rep_2_Ct)), collapse = "-"),
         eRNA = paste(sort(c(eRNA_rep_1_Ct,
                             eRNA_rep_2_Ct,
                             eRNA_rep_2_Ct)), collapse = "-")) %>% 
  select(container_id, days_since, common_name, eDNA, eRNA) %>% 
  separate(eDNA, into = paste0('DNArep', 1:3),sep = '-') %>% 
  separate(eRNA, into = paste0('RNArep', 1:3),sep = '-') %>% 
  mutate(across(DNArep1:RNArep3, ~ifelse(.x == 99, 0, .x))) %>%
  mutate(DNArep1 = ifelse(DNArep1 == '', NA, 
                          DNArep1))
spp %>% names
spp$days_since %>% summary
spp$common_name %>% table
## reps -------
repcqDNA <- spp %>%
  select(container_id:DNArep3) %>% 
  pivot_longer(DNArep1:DNArep3, names_to = 'rep', values_to = 'eDNA') %>% 
  mutate(rep = sub('DNArep', '', rep))


repcqRNA <- spp %>%
  select(container_id, days_since, common_name, RNArep1:RNArep3) %>% 
  pivot_longer(RNArep1:RNArep3, names_to = 'rep', values_to = 'eRNA') %>% 
  mutate(rep = sub('RNArep', '', rep))


reps_DNA <- repcqDNA %>% left_join(repcqRNA) %>% 
  mutate(days_since = case_when(
    days_since == 0 ~ '0day',
    days_since < 8 ~ '7days',
    days_since >= 8 ~ '8days'
  ),
  eDNA_rank = case_when(
    eDNA == 0 ~ 'zero',
    eDNA < 35 ~ 'high',
    eDNA < 40 ~ 'med',
    # cqDNA < 40 ~ 'med-low',
    eDNA >= 40 ~ 'low'
  ),
  eRNA_rank = case_when(
    eRNA == 0 ~ 'zero',
    eRNA < 30 ~ 'high',
    eRNA < 40 ~ 'med',
    # cqRNA < 40 ~ 'med-low',
    eRNA >= 40 ~ 'low'
  )) 



## goods -------------
goods_sum <- history %>%
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

## khapra -----------



khapra <- spp %>% 
  mutate(pest = rowSums(apply(spp[4:9], 2, function(x) x!=0))>0,
         pest = ifelse(pest, 'present', 'absent'),
         pest_present = rowSums(apply(spp[7:9], 2, function(x) x!=0))>0,
         pest_present = ifelse(pest_present, 'present', 'absent')) %>% 
  select(container_id, common_name, pest, pest_present)

## container --------
container <- container_trial %>% 
  select(Container.number, Container.size, Container.owner,
         Container.grade, Date.of.manufacture,
         Age.at.earliest.sample.date) %>%
  rename(age = Age.at.earliest.sample.date) %>% 
  mutate(Container.age = dmy('08/08/2024') - dmy(Date.of.manufacture),
         Container.age = round(as.numeric(Container.age)/365.25, 2),
         Container.age = ifelse(Container.age < 10, '10yrs&under', '10yrs&older'),
         Container.grade = sub(' \\(steel floor\\)', '', Container.grade),
         Container.grade = sub('General purpose\\/', '', Container.grade),
         Container.grade = ifelse(Container.grade %in% c('Cotton quality', 
                                                         'Flexi tank',
                                                         'Unknown'),
                                  'uncategorized', Container.grade),
         Container.owner = ifelse(!Container.owner %in% c('TEXTAINER EQUIPMENT MANAGEMENT LTD',
                                                          'HAPAG LLOYD A.G',
                                                          #   'ORIENT OVERSEAS CONTAINER LINE LTD.', # removed because does not give probabilities for all good senarios
                                                          'TRITON CONTAINER INTERNATIONAL LTD',
                                                          'MAERSK A/S'),
                                  'z_other_owner', Container.owner),
         Container.size = ifelse(is.na(Container.size), Container.size,
                                 paste0(Container.size, 'ft')),
         Container.owner = str_split(Container.owner, 
                                     pattern = ' ', simplify = T)[,1]) %>% 
  select(-Date.of.manufacture) %>% 
  setNames(tolower(gsub('\\.', '_', names(.)))) %>% 
  rename(container_id = container_number) 
container


# model data ------------
model_data <- reps_DNA %>% 
  left_join(khapra) %>% 
  left_join(goods_sum) %>% 
  left_join(container) %>% 
  mutate_if(is.factor, as.character) %>% 
  as.data.frame %>% 
  filter(container_id != 'TCNU1629061') %>% 
  #mutate(eRNA = ifelse(eRNA == 'zero', 'absent', 'present')) %>%
  mutate(eDNA = as.numeric(eDNA),
         eRNA = as.numeric(eRNA),
         goods = factor(goods, levels = c("other_goods", "timber", "food")),
         eDNA_rank = factor(eDNA_rank,
                            levels = c("zero", "low",  "med",  "high")),
         eDNA_binary = ifelse(eDNA >0, 1, 0),
         present = ifelse(pest_present == 'absent', 0, 1),
         goods = ifelse(goods == 'other_goods', 
                        'other', 'food/timber'),
         # container_grade = ifelse(container_grade == 'General purpose',
         #                          'uncategorized', container_grade),
         eDNA_rank_no = as.numeric(eDNA_rank)-1) %>% 
  filter(complete.cases(goods), complete.cases(container_age))
names(model_data)  

write.csv(model_data, './output/curated_model_data.csv', row.names = F)

