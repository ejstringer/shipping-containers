
# load --------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
pest_detection_spp <- read.csv('./output/pest_detection_methods_species.csv')
container <- read.csv('./output/pest_container_specs.csv')



# models SPECIFIC---------------------


## level 1: 0/1 DNA -------------------------------------------------------
### data ------

# with metabarcoding simpson
pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, 
                        levels = unique(container$container_grade)[c(2,3,1)])) %>% 
  filter(complete.cases(age), complete.cases(container_size)) %>% 
  mutate(food = ifelse(goods_risk == 'food', 1, 0),
         seed = ifelse(seed == 'no seed', 0, 1))

### glm model -----
mFull <- glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest + no_countries, data = pest,
             family = binomial) 
summary(mFull)
anova(mFull)

m <- glm(dna ~ grade + simpsondna_nopest, data = pest,
         family = binomial) 
summary(m)


# Next step: Simpson as response ------------------------------------------

### glm model -----
mFull_simpson <- lm(simpsondna_nopest ~ grade + age + food + seed + wood + soil, 
                     data = pest) 
summary(mFull_simpson)
anova(mFull_simpson)

m <- lm(simpsondna_nopest ~ grade + age, data = pest) 
summary(m)


# species specific separated by species model -----------------------------

pest_spp <- pest_detection_spp %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, 
                        levels = unique(container$container_grade)[c(2,3,1)])) %>% 
  filter(complete.cases(age), complete.cases(container_size)) %>% 
  mutate(food = ifelse(goods_risk == 'food', 1, 0),
         seed = ifelse(seed == 'no seed', 0, 1))

pest_spp %>% group_by(common_name) %>% summarise(dna = sum(dna))

glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest, 
             data = filter(pest_spp, common_name=="khapra beetle"),
             family = binomial) %>% summary()

glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest, 
    data = filter(pest_spp, common_name=="asian spongy moth"),
    family = binomial) %>% summary()

glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest, 
    data = filter(pest_spp, common_name=="brown marmorated stink bug"),
    family = binomial) %>% summary()

glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest, 
    data = filter(pest_spp, common_name=="electric ant"),
    family = binomial) %>% summary()

glm(dna ~ grade + age + food + seed + wood + soil + simpsondna_nopest, 
    data = filter(pest_spp, common_name=="spotted lantern fly"),
    family = binomial) %>% summary()



