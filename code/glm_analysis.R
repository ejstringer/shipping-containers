
source('./code/libraries.R')
# load --------------------------------------------------------------------

sppdata <- read.csv('./output/curated_model_data.csv')

species <- sort(unique(sppdata$common_name))
species

model_data <- filter(sppdata, common_name == species[3])

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
  theme(legend.position = c(0.26,0.815),
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

summary(m2)
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
  relocate(model) %>% 
  flextable() %>% 
  autofit() %>% 
  italic(j = 1) %>% 
  border_remove() %>% 
  hline(i = c(5), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header', 
        border = fp_border(color = "grey40", width = 2))-> fxtb
fxtb  
save_as_docx(fxtb, path = './figures/model_selection.docx')


# null --------------------------------------------------------------------



# m1 -----------------------------------------------------------
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


# full --------------------------------------------------------------------



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


## casefile -------------
case_file <- model_data %>%
  replace(is.na(.), values = '*') %>%
  mutate(goods = ifelse(goods == 'other', 'other_goods', goods),
         grade = ifelse(container_grade == 'Food quality', container_grade,
                        'other_grade')) %>%
  mutate_if(is.character, as.factor) %>%
  select(pest_present, eDNA_rank, goods, grade, risk_country, container_age)
names(case_file)
summary(case_file)

write.table(case_file, './output/casefiles/casefile_khapra.txt', 
            sep = '\t',  row.names = F)

# // ~->[CASE-1]->~
