# models -------------------------------------------------------------------
source('./code/libraries.R')
DNA_cq_sample <- read.csv('./output/pest_dna_cq_rna.csv')
pest_detection <- read.csv('./output/pest_detection_methods.csv')
container <- read.csv('./output/pest_container_specs.csv')

## DNA cq -------
dna <- DNA_cq_sample %>% 
  select(container_id,pest_rna, min_cq, common_name,
         conc_2_ng_ul, x260_230, x260_280) %>% 
  unique() %>% 
  # mutate(#x260_230 = ifelse(x260_230 > 1.8 & x260_230 < 2.21, 'pure', 'low quality'),
  #        x260_280 = ifelse(x260_280 > 1.7 & x260_280 < 2.01, 'pure', 'low quality')) %>% 
  rename(cq = min_cq) %>% 
  left_join(container)
dna %>% head
lm(cq ~ conc_2_ng_ul+x260_230+x260_280+age, 
   data = dna) %>% summary
plot(cq ~ x260_230, data = dna, col = 'grey', pch = 16)
abline(lm(cq ~ x260_230, data = dna), col = 'black', lty = 2, lwd = 3)
boxplot(dna$x260_230)
boxplot(dna$x260_280)

m <- lm(cq ~ x260_230, data = dna)
x <- seq(0.1,2.1, 0.01)
y <- predict(m, newdata = data.frame(x260_230 = x), se.fit = T, type='response')


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


dna$cq %>% boxplot
boxplot((dna$cq))
ggplot(dna, aes(x260_230, cq))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()
cor(dna$cq, dna$x260_230)
## RNA cq -------------
m <- glm(pest_rna ~ cq, 
         data = dna, family = binomial) 

summary(m)
x <- seq(20,50, 0.1)
ndata <- data.frame(cq = x)
y <- predict(m, newdata=data.frame(cq = x),
             type="response")

ndata 

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


ggplot(cbind(ndata, y), aes(cq, y))+
  geom_ribbon(aes(ymin = right_lwr, ymax = right_upr),
              fill = 'grey', alpha = 0.5)+
  geom_line(colour = 'pink',lwd = 2)+
  geom_point(data = dna, aes(y = pest_rna, x = cq))+
  theme_classic()+
  xlab('rep(1-3) minimum cq')+
  ylab('RNA detected in container')


## DNA 01 -----

pest <- pest_detection %>% 
  left_join(container) %>% 
  mutate(dna = ifelse(pests_sppDNA > 0, 1, 0), 
         grade = factor(container_grade, levels = unique(container$container_grade)[c(2,1,3)])) %>% 
  filter(complete.cases(age), complete.cases(container_size))

pest %>% data.frame %>% head
m <- glm(dna ~ container_grade + age + goods_risk + seed, data = pest,
         family = binomial) 
summary(m)

with(summary(m), 1 - deviance/null.deviance)

mfull <- glm(dna ~ grade +age+goods_risk, data = pest,
             family = binomial) 
summary(mfull)
m <- glm(dna ~ seed, data = pest,
         family = binomial) 
summary(m)
car::Anova(m, test="LR", type="III") 

chisq.test(table(pest$container_grade, pest$dna))

m3 <- glm(dna ~ age , data = pest,
          family = binomial) 

m4 <- glm(dna ~ goods_risk, data = pest,
          family = binomial) 

m5 <- glm(dna ~ container_size, data = pest,
          family = binomial) 

summary(m)

AIC(m, m2, m3,m4,m5) %>%
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         # formula = form,
         model = rownames(.)) %>% 
  arrange(delta) %>% 
  relocate(model)

x <- unique(pest$container_grade)#seq(20,50, 0.1)

ndata <- data.frame(container_grade = x)

#y <- predict(m, newdata=data.frame(cq = x),se.fit = T, type = 'response')
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
                right_upr = ilink(fit_link + (1.96 * se_link)),
                right_lwr = ilink(fit_link - (1.96 * se_link))) %>% 
  mutate_if(is.numeric, round, 3)
## show
ndata

m %>% summary

y <- predict(m, newdata=data.frame(grade = x),
             type="response", se.fit =T)
ndata <- data.frame(grade = x, fit = y$fit, se = y$se.fit) %>% 
  mutate(lwr = fit - (se*1.96),
         upr = fit + (se*1.96))

ggplot(ndata, aes(grade, fit))+
  #geom_jitter(data = pest, aes(y = dna),alpha = 0.5,  width = 0.25, height = 0.075,   colour = 'lightblue')+
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                colour = viridisLite::viridis(4)[2], alpha = 0.7, width = 0, lwd = 1)+
  geom_point(colour = viridisLite::viridis(4)[2],size = 3)+
  theme_classic()+
  ylim(0,0.5)+
  ylab('Probability of detecting DNA')+
  xlab('Container grade')

