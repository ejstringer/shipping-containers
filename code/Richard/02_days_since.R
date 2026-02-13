collection <- read.csv('./output/C05246_collection.csv')

specific <- read.csv('./output/C05246_genetic_species_specific.csv')

arrival <- read.csv('./output/C05246_container_history.csv') %>% 
  filter(sampled) %>% select(container_id, arrival_date) 

spp <- collection %>% 
  left_join(specific) %>% 
  relocate(container_id) %>% 
  select(-visual_lab) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  filter(!sample_id %in% c('X20210721_1510', 'X20210511_0088')) %>%
  left_join(arrival) %>% 
  mutate(days_since = ymd(collection_date)-ymd(arrival_date)) %>% 
  select(container_id, common_name,days_since, eDNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
  pivot_longer(cols = eDNA_rep_1_Ct:eRNA_rep_3_Ct, names_to = 'molecule',
               values_to = 'cq') %>% 
  mutate(molecule = sub('_rep_', '_', molecule),
         molecule = sub('_Ct', '', molecule)) %>% 
  separate(molecule, into = c('molecule', 'rep'),sep = '_')

# days_since_save <- spp %>% select(container_id, days_since) %>% unique()
# write.csv(days_since_save, './output/days_since_arrival.csv',
#           row.names = FALSE)

spp %>% head
spp$days_since %>% boxplot
sppDNA <- filter(spp,
                 days_since > 0)

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

df1 <- sppDNA %>% 
  filter(cq > 0, days_since <40)

df1$x <- as.numeric(df1$days_since)
df1$y <- df1$cq

df1 %>% 
  group_by(common_name, molecule) %>% 
  mutate(text = lm_eqn(x,y)) %>% 
  ungroup() %>% 
  ggplot( aes(x, cq, colour = common_name))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  facet_grid(molecule~common_name) +
  geom_text(x = 26, y = 21, aes(label = text), parse = TRUE,
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
p1
ggsave('./figures/Richard/days_since_arrival_style.png', p1, units = 'cm',
       height = 14, width = 24)
