# load --------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
container <- read.csv('./output/pest_container_specs.csv')

# models META ------------------------------------------------------------------

## level 1a: diversity -------------------------------------------------------
### data ------
pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, 
                        levels = unique(container$container_grade)[c(2,3,1)])) %>% 
  filter(complete.cases(age), complete.cases(container_size)) %>% 
  mutate(food = ifelse(goods_risk == 'food', 1, 0),
         seed = ifelse(seed == 'no seed', 0, 1),
         diversity = richness_DNA,
         metaDNA01 = ifelse(pests_metaDNA>0, 1,0),
         prop = pests_metaDNA/diversity,
         log_pestDNA_r.abu = log(rel_abundance_DNA)) 

### glm model -----
mFull<- lm(metaDNA01 ~ grade + age + food + seed + wood + soil,
       data = pest) 
mFull %>% summary


### flextable ----------
tb <- summary(mFull)$coefficients %>% 
  data.frame %>% 
  mutate(coefficient = rownames(.),
         coefficient = ifelse(coefficient == '(Intercept)', 
                              'Intercept', coefficient),
         coefficient = sub('grade', '', coefficient),
         r2 = summary(mFull)$adj.r.squared,
         r2= ifelse(duplicated(r2), NA, r2)) %>% 
  relocate(coefficient) %>% 
  data.frame(row.names = NULL) %>% 
  setNames(gsub('__','_', tolower(gsub('\\.', '_', names(.))))) %>% 
  rename(p_value = pr__t_) %>% 
  mutate_if(is.numeric, round, 3)


fxtb1a <- tb %>% 
  rename_all(~ gsub("_", " ", .))%>% 
  flextable() %>% 
  autofit() %>% 
  bold(part = 'header') %>% 
  bold(i= which(tb$p_value < 0.05)[-1]) %>% 
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
    value = as_paragraph(('R'),as_sup('2')));fxtb1a


## level 2: pests DNA -------------------------------------------------------

### model--------------
dna <- filter(pest, rel_abundance_DNA >0)
m <- lm(log(rel_abundance_DNA) ~ simpsondna, data = dna) 
m %>% summary
plot(m, which =4)


### flextable ----------
tb <- summary(m)$coefficients %>% 
  data.frame %>% 
  mutate(coefficient = rownames(.),
         coefficient = ifelse(coefficient == '(Intercept)', 
                              'Intercept', coefficient),
         coefficient = sub('grade', '', coefficient),
         r2 = summary(mFull)$adj.r.squared,
         r2= ifelse(duplicated(r2), NA, r2)) %>% 
  relocate(coefficient) %>% 
  data.frame(row.names = NULL) %>% 
  setNames(gsub('__','_', tolower(gsub('\\.', '_', names(.))))) %>% 
  rename(p_value = pr__t_) %>% 
  mutate_if(is.numeric, round, 3)


fxtb2a <- tb %>% 
  rename_all(~ gsub("_", " ", .))%>% 
  flextable() %>% 
  autofit() %>% 
  bold(part = 'header') %>% 
  bold(i= which(tb$p_value < 0.05)[-1]) %>% 
  border_remove() %>% 
  set_caption(caption = 'Level 2 analysis: DNA pest reads relative abundance and diversity') %>% 
  footnote(i = 1, j = 1, 
           value = as_paragraph(
             c("Response variable = relative abundance of pest DNA reads")
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
    value = as_paragraph(('R'),as_sup('2'))) %>% 
  compose(i = 2,j = 1,
    value = as_paragraph(('Simpson index')));fxtb2a

### predict -------
summary(dna)
x <- seq(0,1, 0.1)
z <- rep(seq(0,65,10), each = 11)
y <- predict(m, newdata = data.frame(simpsondna = x), se.fit = T, type='response')


ndata<- data.frame(simpsondna = x,
                   #richness_DNA = z, 
                   fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))
ndata %>% head
### plot -------------
ggplot(ndata, aes(simpsondna, fit))+
  geom_point(data = dna, aes(y = log(rel_abundance_DNA)), alpha = 0.5, size = 2)+
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = 'grey', alpha = 0.5)+
  geom_line(colour = 'pink',lwd = 2)+
  #facet_wrap(~factor(richness_DNA))+
  theme_classic()+
  ylab(expression('Relative abundance of pest reads   '[`log-scale`]))+
  xlab('Simpson index')

ggsave('./figures/pest_focus/meta_level2_analysis_samples.png',
       width = 5, height = 4)


## level 3: RNA -------------------------------------------------------

mx <- min(pest$rel_abundance_DNA[complete.cases(pest$pests_metaRNA) & pest$rel_abundance_DNA>0])
pestrna <- dna %>% filter(complete.cases(pests_metaRNA)) %>% 
  mutate(rna = ifelse(pest_reads_RNA>0, 1,0),
         dna = log(rel_abundance_DNA))

m3 <- glm(rna~dna, pestrna, family = binomial) 
m3%>% 
  summary
plot(pestrna$rna~pestrna$dna)


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
  mutate_if(is.numeric, round, 3) %>% 
  mutate(coefficient= sub('common_name', '', coefficient))

fxtb3a<-tb %>% 
  rename_all(~ gsub("_", " ", .)) %>% 
  flextable() %>% 
  autofit() %>% 
  # set_caption('Response variable = RNA detected (0,1)', ) %>% 
  bold(part = 'header') %>% 
  bold(i = which(tb$p_value<0.05)[-1]) %>% 
  border_remove() %>% 
  set_caption(caption = 'Level 3 analysis: RNA presence and DNA relative abundance') %>% 
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
  compose(i = 2,j = 1,
          value = as_paragraph('DNA reads ', 
                               as_sub('log-scale'))) %>% 
  compose(
    part = "header", j = 6,
    value = as_paragraph(('R'),as_sup('2'), as_sub(' McKelvey & Zavoina')));fxtb3a
### predict -----------

x <- seq(-12,0, 0.1)

ndata <- data.frame(dna = x)
y <- predict(m3, 
             newdata = ndata,
             se.fit = T, type='response')


ndata<- cbind(ndata, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))
ndata %>% head


### plot --------------

ggplot(ndata, aes(dna, fit))+
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = 'grey', colour = 'grey80',
              alpha = 0.5)+
  geom_line(colour = 'pink', lwd = 1.5)+
  geom_point(data = pestrna, aes(y = rna, x = dna))+
  theme_bw()+
  xlab(expression('Relative abundance of pest reads   '[`log-scale`]))+
  ylab('RNA detection probability')+
  theme(legend.position = 'none',
        strip.background = element_blank())

ggsave('./figures/pest_focus/meta_level3_analysis_RNA_species.png',
       height = 5, width = 7)


## save flextables ---------

save_as_docx(fxtb1a, fxtb2a, fxtb3a, path = './figures/pest_focus/meta_model_tables.docx')

