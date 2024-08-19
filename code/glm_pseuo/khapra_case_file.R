source('./code/libraries.R')
history <- read.csv('./output/C05246_container_history.csv') %>% 
  filter(sampled) 

# arrival -----------------------------------------------------------------

container <- read.csv('./data/shipping_meta/Approach rate trial container data summary.csv') 
container_trial <- read.csv('./data/shipping_meta/Approach rate trial container data summary.csv') 


collection <- read.csv('./output/C05246_collection.csv')
arrival <- read.csv('./output/C05246_container_history.csv') %>% 
  filter(sampled) %>% select(container_id, arrival_date) 

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

spp <- read.csv('./output/C05246_collection.csv') %>% 
  left_join(arrival) %>% 
  left_join(specific) %>% 
  filter(sample_method != 'VAC UNDERFLOOR',
         common_name == 'khapra beetle',
         !sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>% 
  mutate(days_since = as.numeric(ymd(collection_date)- ymd(arrival_date))) %>% 
  relocate(container_id) %>% 
  relocate(days_since, .after = collection_date) %>% 
  pivot_longer(cols = eDNA_rep_1_Ct:eRNA_rep_3_Ct, names_to = 'molecule',
               values_to = 'cq') %>% 
  mutate(molecule = sub('_rep_', '_', molecule),
         molecule = sub('_Ct', '', molecule)) %>% 
  separate(molecule, into = c('molecule', 'rep'),sep = '_') %>% 
  pivot_wider(names_from = molecule, values_from = cq)
spp %>% names
spp %>% data.frame %>% head



# find duplicate containers -----------------------------------------------
spp$container_id %>% unique %>% length
spp$sample_id %>% unique %>% length

uniqueID <- unique(spp[,c('container_id', 'sample_id')])
dupContainers <-  uniqueID$container_id[duplicated(uniqueID$container_id)]

spp %>% filter(container_id %in% dupContainers) %>% View

# day since adjust --------------------------------------------------------
daysminus <- spp %>% 
  group_by(container_id) %>% 
  summarise(days_since = mean(days_since)) %>% 
  filter(days_since<0) 

xx  <- left_join(collection[,c('container_id', 'collection_date')], 
                 history) %>% 
  filter(container_id %in% daysminus$container_id) %>% 
  relocate(collection_date, .after = arrival_date)

xx %>% View

adjustedDaySince <- xx %>% filter(!sampled) %>% 
  arrange(container_id, desc(arrival_date)) %>% 
  mutate(sampled = ifelse(duplicated(container_id), F, T),
         days_since2 = as.numeric(ymd(collection_date)- ymd(arrival_date))) %>% 
  relocate(days_since2, .after = arrival_date) %>% 
  filter(sampled) %>% 
  select(container_id, days_since2)
  

# summarise ---------------------------------------------------------------

repcqDNA <-spp %>% filter(eDNA > 0) %>% 
  group_by(container_id) %>% 
  summarise(repDNA = n(),
            cqDNA = min(eDNA))

repcqRNA <-spp %>% filter(eRNA > 0) %>% 
  group_by(container_id) %>% 
  summarise(repRNA = n(),
            cqRNA = min(eRNA))

spp_DNA <- spp %>% left_join(adjustedDaySince) %>%
  relocate(days_since2, .after = days_since) %>% 
  mutate(days_since = ifelse(days_since < 0,
                             days_since2, days_since)) %>% 
  select(-days_since2) %>% 
  group_by(container_id) %>% 
  summarise(days_since = mean(days_since)) %>% 
  mutate(days_since = case_when(
    days_since == 0 ~ 'on the day',
    days_since < 8 ~ 'within the week',
    days_since >= 8 ~ 'longer'
  )) %>% 
  left_join(repcqDNA) %>% 
  left_join(repcqRNA) %>%
  mutate_if(is.numeric, funs(replace_na(., 0)))

spp_DNA %>% View
spp_DNA$days_since %>% table

# add history -------------------------------------------------------------
#https://www.agriculture.gov.au/sites/default/files/images/target_risk_country_for_host_of_khapra_beetle_0.png
view(history)
goods_sum <- history %>% 
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
                          'yes','no')) %>% 
  select(container_id, timber, food,risk_country)
goods_sum$country %>% table
goods_sum$food %>% table
goods_sum$timber %>% table
goods_sum %>% view

spp_DNA_container <- goods_sum %>% 
  left_join(spp_DNA) %>% 
  mutate(timber = ifelse(timber, 'yes', 'no'),
         food = ifelse(food, 'yes', 'no'),
         khapra = ifelse(repDNA >0 | repRNA > 0, 'present', 'absent'),
         #repDNA = ifelse(repDNA %in% c(1,2), '1-2', repDNA),
         #repRNA = ifelse(repRNA %in% c(1,2), '1-2', repRNA),
         repDNA = factor(repDNA),
         repRNA = factor(repRNA),
         cqDNA = case_when(
           cqDNA == 0 ~ 'zero',
           cqDNA < 35 ~ 'high',
           cqDNA < 40 ~ 'med',
          # cqDNA < 40 ~ 'med-low',
           cqDNA >= 40 ~ 'low'
         ),
         cqRNA = case_when(
           cqRNA == 0 ~ 'zero',
           cqRNA < 35 ~ 'high',
           cqRNA < 40 ~ 'med',
          # cqRNA < 40 ~ 'med-low',
           cqRNA >= 40 ~ 'low'
         )) 
spp_DNA_container %>% View
#levels(spp_DNA_container$repDNA) <- c('0rep', '1rep', '2rep', '3rep')
#levels(spp_DNA_container$repRNA) <- c('0rep', '1rep', '2rep', '3rep')
levels(spp_DNA_container$repDNA)
# container owner ---------------------------------------------------------



container %>% names
container_sum <- container %>% 
  select(Container.number, Container.size, Container.owner,
         Container.grade, Date.of.manufacture) %>%
  mutate(Container.age = dmy('08/08/2024') - dmy(Date.of.manufacture),
         Container.grade = sub(' \\(steel floor\\)', '', Container.grade),
         Container.grade = sub('General purpose\\/', '', Container.grade),
         Container.grade = ifelse(Container.grade %in% c('Cotton quality', 
                                                         'Flexi tank',
                                                         'Unknown'),
                                  'uncategorized', Container.grade),
         Container.owner = ifelse(!Container.owner %in% c('TEXTAINER EQUIPMENT MANAGEMENT LTD',
                                                          'HAPAG LLOYD A.G',
                                                          'ORIENT OVERSEAS CONTAINER LINE LTD.',
                                                          'TRITON CONTAINER INTERNATIONAL LTD',
                                                          'MAERSK A/S'),
                                  'other', Container.owner),
         Container.owner = str_split(Container.owner, 
                                     pattern = ' ', simplify = T)[,1]) %>% 
  select(-Date.of.manufacture) %>% 
  setNames(tolower(gsub('\\.', '_', names(.))))
str_split(container_sum$container_owner, pattern = ' ', simplify = T)[,1]
container_sum$container_number %>% length
container_sum$container_size %>% table
container_sum$container_grade %>% table
container_sum$container_owner %>% table %>% sort


# all ---------------------------------------------------------------------

spp_DNA_container_owner <- container %>% 
  select(Container.number, Container.size, Container.owner,
         Container.grade, Date.of.manufacture) %>%
  mutate(Container.age = dmy('08/08/2024') - dmy(Date.of.manufacture),
         Container.age = round(as.numeric(Container.age)/365.25, 2),
         Container.age = ifelse(Container.age < 10, 'less than 10yrs', 'more than 10yrs'),
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
                                  '_other_owner', Container.owner),
         Container.size = ifelse(is.na(Container.size), Container.size,
                                 paste0(Container.size, 'ft')),
         Container.owner = str_split(Container.owner, 
                                     pattern = ' ', simplify = T)[,1]) %>% 
  select(-Date.of.manufacture) %>% 
  setNames(tolower(gsub('\\.', '_', names(.)))) %>% 
  rename(container_id = container_number) %>% 
  left_join(spp_DNA_container) %>% 
  mutate(goods = case_when(
    food == 'yes' ~ 'food',
    timber == 'yes' ~ 'timber',
    !(food == 'yes'| timber == 'yes') ~ 'other_goods'
  )) %>% 
  select(-food, -timber) %>% 
  mutate_if(is.factor, as.character) %>%
  replace(is.na(.), '*') %>% 
  select(-container_id)
spp_DNA_container_owner %>% names
spp_DNA_container_owner %>% glimpse

write.table(spp_DNA_container_owner, './output/case_file_khapra.txt', 
            sep = '\t',  row.names = F)


spp_DNA_container_owner$repDNA %>% table


case_file <- spp_DNA_container_owner[,c(2:4,6:14)]

a <- lapply(case_file, unique )

length(a)
data.frame(a[[1]],
           a[[2]],
           a[[10]],
           a[[11]],
           a[9],
           a[[8]])

case_nodes <- do.call(c, a) %>% 
  data.frame(column = names(.), value = ., row.names = NULL) %>% 
  mutate(column = gsub('[[:digit:]]+', '', column),
         row = row_number()) %>% 
  pivot_wider(names_from = column, values_from = value) %>% 
  as.data.frame() %>% 
replace(is.na(.), values = '*')

case_nodes <- spp_DNA_container_owner[,c(2:4, 10:14)] %>%
  unique %>% glimpse()

write.table(case_nodes, './output/case_nodes_khapra.txt', 
            sep = '\t',  row.names = F)

spp_DNA_container_owner[,c('container_owner')]


# CPT ---------------------------------------------------------------------

test <- spp_DNA_container_owner

table(test$timber, test$food)/2005
table(test$khapra, test$food)/2005
table(test$goods, test$khapra)/2005

table(test$timber, test$risk_country)

table(test$goods, test$container_age) # remove Orient
table(test$risk_country, test$container_owner)
table(test$container_size, test$container_owner)
table(test$container_age, test$container_owner)

table(test$cqDNA, test$goods)
table(test$cqRNA, test$goods)

table(test$cqDNA, test$khapra)
table(test$cqRNA, test$khapra)

table(test$cqDNA, test$repDNA)
table(test$cqRNA, test$repRNA)


table(test$repRNA, test$repDNA)
table(test$cqRNA, test$repRNA)



# NEW data structure ------------------------------------------------------


daysminus <- spp %>% 
  group_by(container_id) %>% 
  summarise(days_since = mean(days_since)) %>% 
  filter(days_since<0) 

xx  <- left_join(collection[,c('container_id', 'collection_date')], 
                 history) %>% 
  filter(container_id %in% daysminus$container_id) %>% 
  relocate(collection_date, .after = arrival_date)

adjustedDaySince <- xx %>% filter(!sampled) %>% 
  arrange(container_id, desc(arrival_date)) %>% 
  mutate(sampled = ifelse(duplicated(container_id), F, T),
         days_since2 = as.numeric(ymd(collection_date)- ymd(arrival_date))) %>% 
  relocate(days_since2, .after = arrival_date) %>% 
  filter(sampled) %>% 
  select(container_id, days_since2)


## spp ---------------
spp <- read.csv('./output/C05246_collection.csv') %>% 
  left_join(arrival) %>% 
  left_join(specific) %>% 
  filter(sample_method != 'VAC UNDERFLOOR',
         common_name == 'khapra beetle',
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
  select(container_id, days_since, eDNA, eRNA) %>% 
  separate(eDNA, into = paste0('DNArep', 1:3),sep = '-') %>% 
  separate(eRNA, into = paste0('RNArep', 1:3),sep = '-') %>% 
  mutate(across(DNArep1:RNArep3, ~ifelse(.x == 99, 0, .x))) %>% 
  left_join(adjustedDaySince) %>%
  relocate(days_since2, .after = days_since) %>% 
  mutate(days_since = ifelse(days_since < 0,
                             days_since2, days_since)) %>% 
  select(-days_since2) %>% 
  mutate(DNArep1 = ifelse(DNArep1 == '', NA, 
                          DNArep1))
  

spp %>% names
#spp %>% view

## reps -------
repcqDNA <- spp %>%
  select(container_id:DNArep3) %>% 
  pivot_longer(DNArep1:DNArep3, names_to = 'rep', values_to = 'eDNA') %>% 
  mutate(rep = sub('DNArep', '', rep))


repcqRNA <- spp %>%
  select(container_id, days_since, RNArep1:RNArep3) %>% 
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
reps_DNA  

## goods -------------
goods_sum <- history %>% 
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


## khapra -----------



  khapra <- spp %>% 
  mutate(khapra = rowSums(apply(spp[3:8], 2, function(x) x!=0))>0,
         khapra = ifelse(khapra, 'present', 'absent'),
         stillkhapra = rowSums(apply(spp[6:8], 2, function(x) x!=0))>0,
         stillkhapra = ifelse(stillkhapra, 'present', 'absent')) %>% 
  select(container_id, khapra, stillkhapra)

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
## join ALL ----------

case_file <- reps_DNA %>% 
  left_join(khapra) %>% 
  left_join(goods_sum) %>% 
  left_join(container) %>% 
  mutate_if(is.factor, as.character) %>% 
  as.data.frame %>% 
  replace(is.na(.), values = '*') %>%
  filter(container_id != 'TCNU1629061') %>% 
  #mutate(eRNA = ifelse(eRNA == 'zero', 'absent', 'present')) %>%
  select(-container_id) 
names(case_file)  

write.table(case_file, './output/casefiles/casefile_reps_khapra.txt', 
            sep = '\t',  row.names = F)
# // ~->[CASE-1]->~
# figures ------------
case_file %>% 
  group_by(risk_country) %>% 
  summarise(kh = sum(khapra == 'present')) %>% 
  # mutate(risk_country = ifelse(risk_country == '*',
  #                              'unknown', risk_country)) %>% 
  ggplot(aes(x=risk_country, y=kh))+
  geom_bar(stat = 'identity') +
  theme_classic()
  
table(case_file$eRNA, case_file$eDNA)

  
df3<- repcqDNA %>% 
  left_join(repcqRNA) %>% 
  left_join(goods_sum) %>% 
  mutate(eDNA = as.numeric(eDNA),
         eRNA = as.numeric(eRNA)) %>% 
  filter(complete.cases(goods)) 

  ggplot(filter(df3, eDNA >0),
         aes(risk_country, eDNA, fill = paste(risk_country,
                                              goods)))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c(rep('khaki',3), 'grey', rep('khaki',2)))+
  facet_grid('Khapra Beetle'~goods)+
  theme(legend.position = 'none')+
  xlab('High Risk Country')
  
    ggsave(file = './figures/khapra_eDNA_country.png')

  ggplot(filter(df3,eRNA >0),
         aes(risk_country, eRNA, colour = paste(risk_country,
                                              goods)))+
  geom_boxplot(fill = 'khaki')+
  theme_bw()+
  scale_colour_manual(values = c(rep('black',2), 'grey50', rep('black',3)))+
  facet_grid('Khapra Beetle'~goods)+
  theme(legend.position = 'none')+
  xlab('High Risk Country') 
  
    ggsave('./figures/khapra_eRNA_country.png')


df2<-repcqDNA %>% 
  left_join(repcqRNA) %>% 
  left_join(container) %>%
  left_join(goods_sum) %>% 
  mutate(eDNA = as.numeric(eDNA)) %>% 
  filter(eDNA > 0, complete.cases(goods),
         complete.cases(container_size)) 

  ggplot(df2, aes(container_size, eDNA))+
  geom_boxplot(fill = 'lightblue')+
  theme_bw()+
  facet_grid('Khapra Beetle'~goods) 
    ggsave('./figures/khapra_eDNA_size.png')
  
    ggplot(df2, aes(container_age, eDNA))+
    geom_boxplot(fill = 'lightgreen')+
    theme_bw()+
    facet_grid('Khapra Beetle'~goods) 
      ggsave('./figures/khapra_eDNA_age.png')
  

case_file %>% 
  filter(#container_age=='10yrs&under',
         container_size== '40ft',
         risk_country == 'yes', 
         goods == 'food')

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
         present = ifelse(stillkhapra == 'absent', 0, 1),
         goods = ifelse(goods == 'other_goods', 
                        'other', 'food/timber'),
         # container_grade = ifelse(container_grade == 'General purpose',
         #                          'uncategorized', container_grade),
         eDNA_rank_no = as.numeric(eDNA_rank)-1) %>% 
  filter(complete.cases(goods), complete.cases(container_age))
names(model_data)  


ggplot(model_data %>% filter(complete.cases(goods)),
       aes(x = eDNA_rank_no,
           fill = goods))+
  geom_histogram(bins =  4, colour = 'black')+
  theme_classic()+
  facet_grid(goods~risk_country,
             scale = 'free')


m <- glm(eDNA_binary ~  goods+risk_country, data = model_data, family = binomial)
summary(m, corr = FALSE)

m <- glm(present ~  eDNA_rank_no, data = model_data, family = binomial)
summary(m, corr = FALSE)

m <- glm(present ~ eDNA_rank_no + risk_country * goods, data = model_data, family = binomial)
summary(m, corr = FALSE)

m <- glm(present ~ eDNA_rank_no + risk_country * goods, data = model_data, family = binomial)
summary(m, corr = FALSE)

m <- glm(present ~ age + container_grade, data = model_data, family = binomial)
summary(m, corr = FALSE)


ggplot(model_data,
       aes(goods, fill = goods))+
  geom_bar(aes(y = (..count..)/sum(..count..)),
           colour = 'black')+
  ylab('frequency')+
  theme_classic()
model_data %>% 
ggplot(aes(risk_country, fill = risk_country))+
  geom_bar(colour = 'black')+
  ylab('frequency')+
  theme_classic()+
  theme(legend.position = 'none')+
  facet_grid(~goods)

plogis(-12.42)

ggplot(filter(model_data, eDNA >0,
              complete.cases(goods)),
       aes(x = eDNA, 
           colour = risk_country))+
  geom_histogram()+
   theme_bw()+
  facet_grid(risk_country~goods)+
  theme(legend.position = 'right')


# barplots ----------------------------------------------------------------


model_data %>% 
  group_by(goods, risk_country) %>% 
  summarise(value = n()/nrow(.)) %>% 
  ggplot(aes(y = value,
             fill = risk_country, x = goods))+
  geom_col(position = "dodge") +
  geom_text(aes(label = round(value,3)), vjust = -0.2,
            position = position_dodge(width = 0.9),
            size = 3) +
  ylab("Frequency") +
  xlab("Transporting cargo") +
  scale_fill_manual(values = c("lightblue",'lightblue4'),
                    name = "High Risk Country") +
  theme_light()+
  ggtitle('Port of Brisbane')+
  theme(legend.position = c(0.74,0.815),
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(colour = 'grey'))-> p1;p1

ggsave('./figures/frequency_goods_country.png',p1,
       units = 'cm', width = 9,
       height = 9)


model_data %>% 
  group_by(eDNA_rank) %>% 
  summarise(value = n()/nrow(.)) %>% 
  ggplot(aes(y = value,
             x = eDNA_rank, fill = eDNA_rank))+
  geom_col(position = "dodge") +
  geom_text(aes(label = round(value,3)), vjust = -0.2,
            position = position_dodge(width = 0.9),
            size = 3) +
  ylab("Frequency") +
  xlab("DNA cq level") +
  scale_fill_manual(values = c('grey',
                               "palegreen3",
                               "navajowhite1",
                               'salmon1'),
                    name = "DNA cq") +
  theme_light()+
  theme(legend.position = 'none', #c(0.85,0.73),
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(colour = 'grey'))-> p2;p2

ggsave('./figures/frequency_DNAcq_levels.png',p2,
       units = 'cm', width = 9,
       height = 8)

model_data %>% 
  group_by(present) %>% 
  summarise(value = n()/nrow(.)) %>% 
  ggplot(aes(y = value,
             x = present, fill = factor(present)))+
  geom_col(position = "dodge") +
  geom_text(aes(label = round(value,3)), vjust = -0.2,
            position = position_dodge(width = 0.9),
            size = 3) +
  ylab("Frequency") +
  xlab("") +
  scale_fill_manual(values = c('grey',
                               "palegreen3",
                               "navajowhite1",
                               'salmon1'),
                    name = "DNA cq") +
  theme_light()+
  theme(legend.position = 'none', #c(0.85,0.73),
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(colour = 'grey'))


# models Hypothesised -----------------------------------------------------
plot(model_data$age, model_data$eDNA)

full <- glm(present ~ eDNA_rank_no + goods *risk_country + container_grade + age,
            data = model_data, family = 'binomial')
null <- glm(present ~ eDNA_rank_no,
          data = model_data, family = 'binomial')

m3 <- glm(present ~ eDNA_rank_no + goods *risk_country +  age,
          data = model_data, family = 'binomial')

m2 <- glm(present ~ eDNA_rank_no + goods *risk_country + container_grade,
          data = model_data, family = 'binomial')

m1 <- glm(present ~ eDNA_rank_no + goods *risk_country,
          data = model_data, family = 'binomial')

m5 <- glm(present ~ eDNA_rank_no +  container_grade,
          data = model_data, family = 'binomial')

m6 <- glm(present ~ eDNA_rank_no + container_grade,
          data = model_data, family = 'binomial')

m7 <- glm(present ~ eDNA_rank_no + container_grade + goods,
          data = model_data, family = 'binomial')

m8 <- glm(present ~ eDNA_rank_no,
          data = model_data, family = 'binomial')


# null <- glm(present ~ 1,
#           data = model_data, family = 'binomial')
AIC(m1, m2)
aictb <- AIC(m1,m2, m3, null, full)
form <- sapply(list(m1, m2,m3, null, full),
               function(x)
                  as.character(formula(x))[3])
aictb %>%
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         formula = form,
         model = rownames(aictb)) %>% 
  arrange(delta) %>% 
  relocate(model) 
  flextable() %>% 
  autofit() %>% 
  italic(j = 1) %>% 
  border_remove() %>% 
  hline(i = c(5), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header', 
        border = fp_border(color = "grey40", width = 2)) %>% 
  save_as_docx(., path = './figures/model_selection.docx')
full %>% summary
m3 %>% summary
glm(present ~ eDNA_rank_no + risk_country * container_grade,
          data = model_data, family = 'binomial') %>% 
  summary(., corr = FALSE)
# dragon figure -----------------------------------------------------------
m2 <- glm(present ~ eDNA_rank_no + risk_country * goods,
          data = model_data, family = 'binomial')
summary(m2, corr = FALSE)
plogis(m2$coefficients)
flextable(data.frame(coef = names(m2$coefficients),
                     propability = plogis(m2$coefficients)))
gedplot_model <- data.frame(eDNA_rank_no = seq(0,3.49,0.01),
                            risk_country = rep(c("yes", 'no'), each = 350*2),
                            goods = rep(c('other','food/timber'), each = 350)
) %>%   #gedplot %>% 
  # mutate(`Trapping grid` = sub('side', '', paste(site, grid)),
  #        `Survey location` = sub('side', '', paste(site, grid))) %>% 
  mutate(present = predict(m2, newdata=data.frame(eDNA_rank_no, 
                                                    risk_country, 
                                                    goods),
                             type="response"))
summary(model_data)

myCol <- palette.colors(palette = "Okabe-Ito")[c(6,7,8,5)]
#viridisLite::viridis(4) 
          # c('grey50',
          #  "palegreen3",
          #  "navajowhite1",
          #  'salmon1') #viridisLite::viridis(4)[c(3,1)]#rep(c('#009503', '#65019F'), each = 1)

  ggplot(model_data, 
         aes(x = eDNA_rank_no, y = present,
             colour = goods))+
  geom_jitter(height = 0,#aes(shape = risk_country), 
              colour = 'grey60', width = 0.3,
              size = 1, alpha = 0.2)+
  geom_line(data = gedplot_model, size = 1,
            aes(linetype = risk_country))+
  scale_colour_manual(values = myCol,
                      name = 'Transporting')+
  scale_linetype_manual(values = c(2,1,2,1),
                        name = 'High Risk Country')+
  theme_bw()+
    #geom_hline(yintercept = 0.5, colour = 'red',alpha = 0.5)+
  # guides(colour = guide_legend(override.aes = list(shape = NA,
  #                                                  size = 2)))+
  ylab('Detection probability')+
  xlab('DNA cq level')+
  theme(legend.position = c(0.20,0.65),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.3, units = 'cm'),
        legend.background = element_rect(colour = 'black',),
        legend.key.width = unit(1.2, units = 'cm'),
        legend.key.height = unit(0.4, units = 'cm')
  )-> pp2 ;pp2

  ggsave('./figures/khapra_probability.png',
         plot = pp2, units = 'cm', width = 14,
         height = 9)  
 
  m2 <- glm(present ~ eDNA_rank_no + risk_country*goods + container_grade + age,
            data = model_data, family = 'binomial')
  summary(m2, corr = FALSE)
  anova(m2)
  
  gedplot_model <- data.frame(age = seq(0,24.9,0.1),
                              risk_country = 'no',
                              goods = 'food/timber',
                              container_grade = rep(unique(model_data$container_grade),
                                                    each = 250*4),
                              eDNA_rank_no = rep(0:3, each = 250)) %>% 
    mutate(present = predict(m2, newdata=data.frame(eDNA_rank_no,
                                                    risk_country,
                                                    goods,
                                                    age, 
                                                    container_grade),
                             type="response"))
  dna.labs <- paste('eDNA lvl:', 0:3)
  names(dna.labs) <- as.character(0:3)
  ggplot(filter(model_data), 
         aes(x = age, y = present,
             colour = container_grade))+
    geom_jitter(height = 0, 
                colour = 'grey60', width = 0.3,
                size = 1, alpha = 0.2)+
    geom_line(data = gedplot_model, size = 1)+
    scale_colour_manual(values = myCol,
                        name = 'Container grade')+
    theme_bw()+
    facet_wrap(~eDNA_rank_no, ncol = 4,
               labeller= labeller(eDNA_rank_no = dna.labs))+
    #geom_hline(yintercept = 0.5, colour = 'red',alpha = 0.5)+
    # guides(colour = guide_legend(override.aes = list(shape = NA,
    #                                                  size = 2)))+
    ylab('Detection probability')+
    xlab('Age of container')+
    theme(legend.position = c(0.17,0.77),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.3, units = 'cm'),
          legend.background = element_rect(colour = 'black',),
          legend.key.width = unit(1.2, units = 'cm'),
          legend.key.height = unit(0.4, units = 'cm')
    )-> pp3 ;pp3
  
  ggsave('./figures/khapra_probability_agegrade.png',
         plot = pp3, units = 'cm', width = 16,
         height = 11)  
  