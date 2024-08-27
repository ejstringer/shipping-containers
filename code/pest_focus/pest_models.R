

# load --------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
container <- read.csv('./output/pest_container_specs.csv')
# congruence --------------------------------------------------------------
## data ----------------

# DNA
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

# RNA
pestcount2 <-table(pest_detection$pests_sppRNA, 
                   pest_detection$pests_metaRNA,
                   dnn = c('pests_sppRNA', 'pests_metaRNA')) %>% 
  as.data.frame 

pest_detectionRNA <- pest_detection %>% 
  filter(!is.na(pests_metaRNA)) %>% 
  mutate_all(as.factor) %>% 
  left_join(pestcount2) 

pRNA <-  pest_detectionRNA %>%  
  mutate(xx = paste(pests_sppRNA, pests_metaRNA)) %>% 
  filter(!duplicated(xx), Freq >1)

## plot DNA ---------
  ggplot(pest_detectionDNA, aes(pests_sppDNA, pests_metaDNA, 
             colour = factor(pests_meta5DNA)))+
  geom_jitter(size = 2.5,alpha = 0.7, width = 0.2, height = 0.2)+
  theme_bw()+
  theme(legend.position = c(.85,.8),
        legend.background = element_rect(colour = 'grey'))+
  geom_text(aes(label = Freq), colour = 'black',size = 3, 
            nudge_y = 0.3, nudge_x = -0.15, 
            data = pDNA) +
  guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific DNA detections')+
  ylab('Metabarcoding DNA detections')

ggsave('./figures/pest_counts_comparison_DNA.png')

pest_detectionDNA %>%
  select(container_id, pests_sppDNA, pests_metaDNA) %>% 
  pivot_longer(pests_sppDNA:pests_metaDNA, names_to = 'method',
               values_to = 'pests') %>% 
  mutate(pests = as.numeric(as.character(pests)),
         Method = ifelse(grepl('meta', method), 'Metabarcoding','Spp specific'),
         method = sub('pests_', '', sub('DNA','', method))) %>% 
  ggplot(aes(x = pests, fill = Method))+
  geom_histogram(colour = 'black', position = position_dodge(), binwidth = 0.5)+
  theme_classic()+
  scale_fill_manual(values = c('grey50', 'grey90'))+
  xlab('Pests detected with eDNA')+
  ylab('Frequency')

ggsave('./figures/pest_counts_hist_DNA.png')

## data RNA -------

pest_detectionRNA %>% 
  #filter(!is.na(pests_metaRNA))
  mutate(pests_sppRNA = factor(pests_sppRNA, 
                               levels = c('0', '1', '2'))) %>% 
  ggplot(aes(pests_sppRNA, pests_metaRNA, colour = pests_meta5RNA))+
  geom_jitter(size = 2.5,alpha = 0.6, width = 0.1, height = 0.1)+
  theme_bw()+
  geom_text(aes(label = Freq), colour = 'black',size = 5, 
            nudge_y = 0.3, nudge_x = -0.1, 
            data = pRNA) +
  scale_x_discrete(drop = FALSE)+
 guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific RNA detections')+
  ylab('Metabarcoding RNA detections')

ggsave('./figures/pest_counts_comparison_RNA.png')

pest_detectionRNA %>%
  select(container_id, pests_sppRNA, pests_metaRNA) %>% 
  pivot_longer(pests_sppRNA:pests_metaRNA, names_to = 'method',
               values_to = 'pests') %>% 
  mutate(pests = as.numeric(as.character(pests)),
         Method = ifelse(grepl('meta', method), 'Metabarcoding','Spp specific'),
         method = sub('pests_', '', sub('RNA','', method))) %>% 
  ggplot(aes(x = pests, fill = Method))+
  geom_histogram(colour = 'black', position = position_dodge(), 
                 binwidth = 0.5)+
  scale_fill_manual(values = c('grey50', 'grey90'))+
  theme_classic()+
  xlab('Pests detected with eRNA')+
  ylab('Frequency')

ggsave('./figures/pest_counts_hist_RNA.png')
# models ------------------------------------------------------------------


## level 1: 0/1 DNA -------------------------------------------------------
### data ------
pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, 
                        levels = unique(container$container_grade)[c(2,3,1)])) %>% 
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

### flextable ----------
tb <- summary(mFull)$coefficients %>% 
  data.frame %>% 
  mutate(coefficient = rownames(.),
         coefficient = ifelse(coefficient == '(Intercept)', 
                              'Intercept', coefficient),
         coefficient = sub('grade', '', coefficient),
         r2_MckZav = DescTools::PseudoR2(mFull, which = 'McKelveyZavoina'),
         r2_MckZav = ifelse(duplicated(r2_MckZav), NA, r2_MckZav)) %>% 
  relocate(coefficient) %>% 
  data.frame(row.names = NULL) %>% 
  setNames(gsub('__','_', tolower(gsub('\\.', '_', names(.))))) %>% 
  rename(p_value = pr__z_) %>% 
  mutate_if(is.numeric, round, 3)


fxtb1 <- tb %>% 
  rename_all(~ gsub("_", " ", .))%>% 
  flextable() %>% 
  autofit() %>% 
  bold(part = 'header') %>% 
  bold(i= 2:3) %>% 
  border_remove() %>% 
  set_caption(caption = 'Level 1 analysis: DNA presence and shipping container characteristics') %>% 
  footnote(i = 1, j = 1, 
           value = as_paragraph(
             c("Response variable = DNA detected (0,1)")
           ),
           ref_symbols = c( " * "),
           part = "header",inline=T) %>% 
  fontsize(size = 9, part = "footer") %>% 
  align(align = 'right', part = 'footer') %>% 
  hline(part = 'header',
        border = fp_border(color = "grey30",
                           width = 3, 
                           style = 'solid')) %>% 
  hline(i = nrow(tb),border = fp_border(color = "grey30",
                                        width = 2, 
                                        style = 'solid')) %>% 
  compose(
    part = "header", j = 6,
    value = as_paragraph(('R'),as_sup('2'), 
                         as_sub(' McKelvey & Zavoina')));fxtb1 

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

ggsave('./figures/pest_level1_analysis_containers.png')

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

### flextable -------

tb <- summary(mfull2)$coefficients %>% 
  data.frame %>% 
  mutate(coefficient = rownames(.),
         coefficient = ifelse(coefficient == '(Intercept)', 
                              'Intercept', coefficient),
         coefficient = gsub('_', ' ', coefficient),
         r2 = summary(mfull2)$adj.r.squared,
         r2= ifelse(duplicated(r2), NA, r2)) %>% 
  relocate(coefficient) %>% 
  data.frame(row.names = NULL) %>% 
  setNames(gsub('__','_', tolower(gsub('\\.', '_', names(.))))) %>% 
  rename(p_value = pr__t_) %>% 
  mutate_if(is.numeric, round, 3)

fxtb2 <- tb %>% 
  rename_all(~ gsub("_", " ", .))%>% 
  flextable() %>% 
  autofit() %>% 
  bold(part = 'header') %>% 
  bold(i = 3) %>% 
  border_remove() %>% 
  set_caption(caption = 'Level 2 analysis: DNA cq and sample characteristics') %>% 
  footnote(i = 1, j = 1, 
           value = as_paragraph(
             c("Response variable = DNA cq")
           ),
           ref_symbols = c( " * "),
           part = "header",inline=T) %>% 
  fontsize(size = 9, part = "footer") %>% 
  align(align = 'right', part = 'footer') %>% 
  hline(part = 'header',
        border = fp_border(color = "grey30",
                           width = 3, 
                           style = 'solid')) %>% 
  hline(i = nrow(tb),border = fp_border(color = "grey30",
                                        width = 2, 
                                        style = 'solid')) %>% 
  compose(
    part = "header", j = 6,
    value = as_paragraph(('R'),as_sup('2'))); fxtb2


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

ggsave('./figures/pest_level2_analysis_samples.png')
## level 3: /1 DNA 0/1 RNA ------------------------------------------------

### glm model ---------

m3 <- glm(pest_rna ~ cq, 
         data = dna, family = binomial) 

summary(m3)

### flextable -------

tb <- summary(m3)$coefficients %>% 
  data.frame %>% 
  mutate(coefficient = rownames(.),
         coefficient = ifelse(coefficient == '(Intercept)', 
                              'Intercept', coefficient),
         r2_MckZav = DescTools::PseudoR2(m3, which = 'McKelveyZavoina'),
         r2_MckZav = ifelse(duplicated(r2_MckZav), NA, r2_MckZav)) %>% 
  relocate(coefficient) %>% 
  data.frame(row.names = NULL) %>% 
  setNames(gsub('__','_', tolower(gsub('\\.', '_', names(.))))) %>% 
  rename(p_value = pr__z_) %>% 
  mutate_if(is.numeric, round, 3)

fxtb3<-tb %>% 
  rename_all(~ gsub("_", " ", .))%>% 
  flextable() %>% 
  autofit() %>% 
 # set_caption('Response variable = RNA detected (0,1)', ) %>% 
  bold(part = 'header') %>% 
  bold(i = 2) %>% 
  border_remove() %>% 
  set_caption(caption = 'Level 3 analysis: RNA presence and DNA cq') %>% 
  footnote(i = 1, j = 1, 
                  value = as_paragraph(
                    c("Response variable = RNA detected (0,1)")
                  ),
                  ref_symbols = c( " * "),
                  part = "header",inline=T) %>% 
  fontsize(size = 9, part = "footer") %>% 
  align(align = 'right', part = 'footer') %>% 
  hline(part = 'header',
        border = fp_border(color = "grey30",
                           width = 3, 
                           style = 'solid')) %>% 
  hline(i = nrow(tb),border = fp_border(color = "grey30",
                           width = 2, 
                           style = 'solid')) %>% 
  compose(
    part = "header", j = 6,
    value = as_paragraph(('R'),as_sup('2'), as_sub(' McKelvey & Zavoina')));fxtb3
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

ggsave('./figures/pest_level3_analysis_RNA.png')


## save flextables ---------

save_as_docx(fxtb1, fxtb2, fxtb3, path = './figures/pest_model_tables.docx')
