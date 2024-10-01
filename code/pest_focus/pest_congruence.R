

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
  mutate(pests_sppDNA_binary = ifelse(pests_sppDNA > 0, 1,0),
         pests_metaDNA_binary = ifelse(pests_metaDNA > 0, 1,0)) %>% 
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
  mutate(pests_sppRNA_binary = ifelse(pests_sppRNA > 0, 1,0),
         pests_metaRNA_binary = ifelse(pests_metaRNA > 0, 1,0)) %>% 
  mutate_all(as.factor) %>% 
  left_join(pestcount2) 

pRNA <-  pest_detectionRNA %>%  
  mutate(xx = paste(pests_sppRNA, pests_metaRNA)) %>% 
  filter(!duplicated(xx), Freq >1)

## plot DNA ---------
colourFriendly <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(pest_detectionDNA, 
       aes(pests_sppDNA, pests_metaDNA, 
           colour = factor(pests_meta5DNA)))+
  geom_jitter(size = 2,alpha = 0.5, width = 0.2, height = 0.2)+
  theme_bw()+
  theme(legend.position = c(.85,.8),
        legend.background = element_rect(colour = 'grey'))+
  geom_text(aes(label = Freq), colour = 'black',size = 3, 
            nudge_y = 0.3, nudge_x = -0.15, 
            data = pDNA) +
  scale_color_manual(values = colourFriendly[c(1,4,7)])+
  guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific DNA detections')+
  ylab('Metabarcoding DNA detections')


ggsave('./figures/pest_focus/pest_counts_comparison_DNA.png',
       width = 5.2, height = 3.5)

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
  ylab('Frequency')+
  theme(legend.position = c(0.8,0.8))

ggsave('./figures/pest_focus/pest_counts_hist_DNA.png',
       width = 4.2, height = 3.5)

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
  scale_color_manual(values = colourFriendly[c(1,4,7)])+
  guides(colour=guide_legend(title='Meta/Spp specific'))+
  xlab('Species specific RNA detections')+
  ylab('Metabarcoding RNA detections')

ggsave('./figures/pest_focus/pest_counts_comparison_RNA.png',
       width = 5, height = 3.5)

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

ggsave('./figures/pest_focus/pest_counts_hist_RNA.png',
       width = 4.2, height = 3.5)

pest_detection %>% head

irr::icc(
  pest_detection[c(2,4)], model = "twoway", 
  type = "agreement", unit = "single"
)

irr::icc(
  mutate_all(pest_detectionRNA[c(3,6)],as.numeric), model = "twoway", 
  type = "agreement", unit = "single"
)

