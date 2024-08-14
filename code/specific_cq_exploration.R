
# congruence

# load --------------------------------------------------------------------

source('./code/libraries.R')

arrival <- read.csv('./output/C05246_container_history.csv') %>% 
  filter(sampled) %>% select(container_id, arrival_date) 

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

spp <- read.csv('./output/C05246_collection.csv') %>% 
  left_join(arrival) %>% 
  left_join(specific) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  filter(!(common_name == 'brown marmorated stink bug' & 
           eRNA_rep_1_Ct > 0 & eRNA_rep_2_Ct == 0 &
           eDNA_rep_1_Ct == 0)) %>% 
  mutate(days_since = ymd(collection_date)- ymd(arrival_date)) %>% 
  relocate(container_id) %>% 
  relocate(days_since, .after = collection_date) %>% 
  pivot_longer(cols = eDNA_rep_1_Ct:eRNA_rep_3_Ct, names_to = 'molecule',
               values_to = 'cq') %>% 
  mutate(molecule = sub('_rep_', '_', molecule),
         molecule = sub('_Ct', '', molecule)) %>% 
  separate(molecule, into = c('molecule', 'rep'),sep = '_')
spp %>% names
spp %>% glimpse()

sppDNA <- filter(spp,
                 days_since > 0) %>% 
  mutate(seed = ifelse(grepl('SEED', visual_contents),
                       'Seed', 'No seed'))


lm_eqn <- function(x, y){
  m <- lm(y ~ x);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
  #                  list(a = format(unname(coef(m)[1]), digits = 2),
  #                       b = format(unname(coef(m)[2]), digits = 2),
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  eq <- substitute(italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


sppDNA %>% names
sppDNA %>% nrow

df1 <- sppDNA %>% 
  filter(cq > 0, days_since <40)

df1$x <- as.numeric(df1$days_since)
df1$y <- df1$cq

df1 %>% 
  group_by(common_name, molecule) %>% 
  mutate(text = lm_eqn(x,y)) %>% 
  ungroup() %>% 
  ggplot( aes(x, cq, colour = common_name))+
  geom_point(alpha = 0.7, colour = 'salmon')+
  theme_bw()+
  facet_grid(common_name~molecule) +
  geom_text(x = 10, y = 16.5, aes(label = text), parse = TRUE,
            size = 3,
            colour = 'black')+
  geom_smooth(method = 'lm', se = T, colour = 'black')+
  xlab('days since arrival')+
  ylab('Cq')+
  theme(legend.position = 'none',
        strip.background = element_rect(fill = 'grey90', 
                                        colour = 'grey90'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) -> p1

ggsave('./figures/days_since_arrival_style.png', p1, units = 'cm',
       height = 24, width = 14)


lmer(cq ~ days_since + molecule +(1|sample_id),
     data = filter(sppDNA, cq>0)) %>% 
  summary

sppDNA %>% 
  filter(cq > 0) %>% 
  ggplot(aes(common_name, cq, fill = molecule))+
  geom_boxplot()+
  theme_bw()

sppDNA %>% 
  pivot_wider(names_from = molecule, values_from = cq) %>% 
  mutate(visual = case_when(
    container_id == "TEMU3100647" ~ "TEMU3100647",
    container_id == "TGBU7753399" ~ "TGBU7753399",
    !container_id %in% c("TEMU3100647", "TGBU7753399") ~ 'other'
  )) %>% 
  filter(eDNA > 0 & eRNA > 0) %>% 
  ggplot(aes(eDNA, eRNA, colour = visual))+
  geom_smooth(method = 'lm', se = T, aes(group = species))+
  geom_point(size = 2, alpha = 0.7)+
  geom_abline(intercept = 0, slope = 1,
              lty = 1, size = 1, colour = 'green',
              alpha = 0.5)+
  theme_bw()+
  #scale_colour_manual(values = c('blue', 'green'))
  facet_grid(~common_name, scale = 'free_y')


