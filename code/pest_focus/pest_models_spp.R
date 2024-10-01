# load --------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
container <- read.csv('./output/pest_container_specs.csv')

# models SPECIFIC----------------------------------------------------------

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
  theme_bw()+
  ylim(0,0.3)+
 # geom_hline(yintercept = 0.19, lty =2)+
  ylab('Probability of detecting DNA')+
  xlab('Container grade')

ggsave('./figures/pest_focus/spp_level1_analysis_containers.png',
       width= 5, height = 4)

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

ggsave('./figures/pest_focus/spp_level2_analysis_samples.png',
       width = 5, height = 4)
## level 3: /1 DNA 0/1 RNA ------------------------------------------------

### glm model ---------


m3full <- glm(pest_rna ~ cq + common_name, 
          data = dna, family = binomial) 
summary(m3full)

m3 <- glm(pest_rna ~ cq, 
         data = dna, family = binomial) 

summary(m3)

### flextable -------

tb <- summary(m3full)$coefficients %>% 
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
  mutate_if(is.numeric, round, 3) %>% 
  mutate(coefficient= sub('common_name', '', coefficient))

fxtb3<-tb %>% 
  rename_all(~ gsub("_", " ", .)) %>% 
  flextable() %>% 
  autofit() %>% 
 # set_caption('Response variable = RNA detected (0,1)', ) %>% 
  bold(part = 'header') %>% 
  bold(i = which(tb$p_value<0.05)[-1]) %>% 
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

ndata <- data.frame(cq = x,
                    common_name = rep(unique(dna$common_name),
                                      each = length(x)))
y <- predict(m3full, 
             newdata = ndata,
             se.fit = T, type='response')


ndata<- cbind(ndata, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))
ndata %>% head


### plot --------------

ggplot(ndata, aes(cq, fit, colour = common_name, fill = common_name))+
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = 'grey', colour = 'grey80',
              alpha = 0.5)+
  geom_line(lwd = 1.5)+
  facet_wrap(~common_name, ncol = 1)+
  geom_point(data = dna, aes(y = pest_rna, x = cq))+
  theme_bw()+
  xlab('DNA cq')+
  ylab('RNA detection probability')+
  theme(legend.position = 'none',
        strip.background = element_blank())

ggsave('./figures/pest_focus/spp_level3_analysis_RNA_species.png', units = 'cm',
       height = 21, width = 14)


## save flextables ---------

save_as_docx(fxtb1, fxtb2, fxtb3, path = './figures/pest_focus/spp_model_tables.docx')
