

# load --------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
container <- read.csv('./output/pest_container_specs.csv')
# congruence --------------------------------------------------------------
pestcount <-table(pest_detection$pests_sppDNA, 
      pest_detection$pests_metaDNA,
      dnn = c('pests_sppDNA', 'pests_metaDNA')) %>% 
  as.data.frame 

pest_detectionDNA <- pest_detection %>% 
  filter(!is.na(pests_metaDNA)) %>% 
  mutate_all(as.factor) %>% 
  left_join(pestcount) 

pDNA <-  pest_detectionDNA %>%  
  mutate(xx = paste(pests_sppDNA, pests_metaDNA)) %>% 
                filter(!duplicated(xx), Freq >1)

  ggplot(pest_detectionDNA, aes(pests_sppDNA, pests_metaDNA, 
             colour = factor(pests_meta5DNA)))+
  geom_jitter(size = 2.5,alpha = 0.7, width = 0.2, height = 0.2)+
  theme_bw()+
  theme(legend.position = c(.85,.8),
        legend.background = element_rect(colour = 'grey'))+
  geom_text(aes(label = Freq), colour = 'black',size = 5, 
            nudge_y = 0.3, nudge_x = -0.15, 
            data = pDNA) +
  guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific DNA detections')+
  ylab('Metabarcoding DNA detections')

  pestcount2 <-table(pest_detection$pests_sppRNA, 
                    pest_detection$pests_metaRNA,
                    dnn = c('pests_sppRNA', 'pests_metaRNA')) %>% 
    as.data.frame 
  
  pest_detectionDNA2 <- pest_detection %>% 
    filter(!is.na(pests_metaDNA)) %>% 
    mutate_all(as.factor) %>% 
    left_join(pestcount2) 
  
  pDNA2 <-  pest_detectionDNA2 %>%  
    mutate(xx = paste(pests_sppDNA, pests_metaDNA)) %>% 
    filter(!duplicated(xx), Freq >1)
  

pest_detection %>% 
  #filter(!is.na(pests_metaRNA))
  ggplot(aes(pests_sppRNA, pests_metaRNA))+
  geom_jitter(size = 2.5,alpha = 0.6, width = 0.1, height = 0.1)+
  theme_bw()+
  coord_cartesian(xlim=c(0, 2))+
  coord_cartesian(ylim=c(0, 2))+
  scale_y_continuous(breaks = 0:3)+
  scale_x_continuous(breaks = 0:3  )
#guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific RNA detections')+
  ylab('Metabarcoding RNA detections')



pest_detection %>% 
  ggplot(aes(factor(pests_sppDNA), pests_meta5DNA, 
             colour = factor(pests_sppRNA)))+
  geom_jitter(aes(size = factor(pests_sppRNA)), alpha = 0.5)+
  theme_classic()+
  scale_y_continuous(limits = c(0,3))

# models ------------------------------------------------------------------


## level 1: 0/1 DNA -------------------------------------------------------
### data ------
pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, levels = unique(container$container_grade)[c(2,3,1)])) %>% 
  filter(complete.cases(age), complete.cases(container_size)) %>% 
  mutate(food = ifelse(goods_risk == 'food', 1, 0),
         seed = ifelse(seed == 'no seed', 0, 1))

### glm model -----
mFull <- glm(dna ~ grade + age + food + seed + wood + soil, data = pest,
         family = binomial) 
summary(mFull)

m <- glm(dna ~ grade, data = pest,
         family = binomial) 
summary(m)
car::Anova(m, test="LR", type="III") 

chisq.test(table(pest$grade, pest$dna))

### predict -------
x <- unique(pest$grade)

y <- predict(m, newdata=data.frame(grade = x),
             type="response", se.fit =T)
ndata <- data.frame(grade = x, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))

### plot ---------
ggplot(ndata, aes(grade, fit))+
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                colour = viridisLite::viridis(4)[2], alpha = 0.7, width = 0, lwd = 1)+
  geom_point(colour = viridisLite::viridis(4)[2],size = 3)+
  theme_classic()+
  ylim(0,0.5)+
  ylab('Probability of detecting DNA')+
  xlab('Container grade')


## level 2: /1 DNA cq -----------------------------------------------------

### data ----------
dna <- DNA_cq_sample %>% 
  select(container_id,pest_rna, min_cq, common_name,
         conc_2_ng_ul, x260_230, x260_280) %>% 
  unique() %>% 
   mutate(pure260_230 = ifelse(x260_230 > 1.8 & x260_230 < 2.21, 'pure', 'low quality'),
          pure260_280 = ifelse(x260_280 > 1.7 & x260_280 < 2.01, 'pure', 'low quality')) %>% 
  rename(cq = min_cq) %>% 
  left_join(container)
dna %>% head

### lm model ---------
mfull2 <- lm(cq ~ conc_2_ng_ul+x260_230+x260_280+age, 
   data = dna)
summary(mfull2)

m2 <- lm(cq ~ x260_230, data = dna)
summary(m2)
summary(m2)$r.squared
### predict -------

x <- seq(0.1,2.1, 0.01)
y <- predict(m2, newdata = data.frame(x260_230 = x), se.fit = T, type='response')


ndata<- data.frame(x260_230 = x, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))
ndata %>% head

ggplot(ndata, aes(x260_230, fit))+
  geom_point(data = dna, aes(y = cq), alpha = 0.5, size = 2)+
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = 'grey', alpha = 0.5)+
  geom_line(colour = 'pink',lwd = 2)+
  
  theme_classic()+
  ylab('DNA cq')+
  xlab('260 / 230 ration')

## level 3: /1 DNA 0/1 RNA ------------------------------------------------

### glm model ---------

m3 <- glm(pest_rna ~ cq, 
         data = dna, family = binomial) 

summary(m3)

### predict -----------

x <- seq(20,50, 0.1)
y <- predict(m3, newdata = data.frame(cq = x), se.fit = T, type='response')


ndata<- data.frame(cq = x, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))
ndata %>% head


### plot --------------

ggplot(ndata, aes(cq, fit))+
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = 'grey', alpha = 0.5)+
  geom_line(colour = 'pink',lwd = 2)+
  geom_point(data = dna, aes(y = pest_rna, x = cq))+
  theme_classic()+
  xlab('DNA cq')+
  ylab('RNA detection probability')




