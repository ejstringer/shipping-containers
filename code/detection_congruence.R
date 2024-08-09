
# congruence

# load --------------------------------------------------------------------

source('./code/libraries.R')

collection <- read.csv('./output/C05246_collection.csv')

spp <- read.csv('./output/C05246_genetic_species_specific.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  filter(!(common_name == 'brown marmorated stink bug' & 
             eRNA_rep_1_Ct > 0 & eRNA_rep_2_Ct == 0 &
             eDNA_rep_1_Ct == 0))

collection <- collection %>% filter(sample_method != 'VAC UNDERFLOOR')

animal <- read.csv('./output/C05246_visual_animal.csv') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method)

rna <- read.csv('./output/C05246_genetic_metabarcoding_rna.csv')%>% 
  filter(PI > 0.99)
dna <- read.csv('./output/C05246_genetic_metabarcoding_dna.csv') %>% 
  filter(PI > 0.99)



# visual detections -------------------------------------------------------

vis_detected <- animal %>% 
  #filter(species %in% unique(spp$species)) %>% 
  mutate(species = ifelse(species %in% unique(spp$species), species,
                          'other'),
         specific = species %in% unique(spp$species)) %>% 
  group_by(container_id, species) %>% 
  summarise(specimens_visual = n(),
            detected_visual = sum(specific) > 0)
vis_detected %>% arrange(desc(species))


# metabarcoding detections ------------------------------------------------

rna_detected <- rna %>% filter(best_hit %in% unique(spp$species)) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = X20210506_0046:X20210811_2063, names_to = 'sample_id',
               values_to = 'reads_rna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_rna = reads_rna > 0) %>% 
  rename(species = best_hit)

rna_detected[,c('container_id', 'species')] %>% duplicated %>% table

names(dna) %>% tail
dna_detected <- dna %>% filter(best_hit %in% unique(spp$species)) %>%
  select_if(grepl('X20', names(.)) | grepl('best_hit', names(.))) %>% 
  group_by(best_hit) %>% 
  summarise_if(is.numeric, sum) %>%
  pivot_longer(cols = X20210504_0002:X20210813_2107, names_to = 'sample_id',
               values_to = 'reads_dna') %>% 
  left_join(collection[,c('sample_id', 'container_id', 'sample_method')]) %>% 
  relocate(container_id) %>% 
  filter(sample_method != 'VAC UNDERFLOOR') %>% 
  select(-sample_id, -sample_method) %>% 
  mutate(detected_dna = reads_dna > 0) %>% 
  rename(species = best_hit) %>% 
  unique()

dna_detected[,c('container_id', 'species')] %>% duplicated %>% table

# overview ----------------------------------------------------------------
spp %>% names
spp %>% head

table(duplicated(collection$sample_id))
no.cont <- (table(collection$container_id))
table(no.cont)
con3 <- no.cont[no.cont>1]

collection[collection$container_id %in% names(con3),]
spp[spp$container_id %in% names(con3),]
spp %>% glimpse()

spp_detected <- spp %>% 
  select(-lab_code) %>% 
  group_by(container_id, common_name, species) %>% 
  summarise_if(is.double, mean) %>% 
  ungroup() %>% 
  mutate(specific_all = select(., eDNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         specific_DNA = select(., eDNA_rep_1_Ct:eDNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         specific_RNA = select(., eRNA_rep_1_Ct:eRNA_rep_3_Ct) %>% 
           rowSums(na.rm = TRUE),
         detected_specDNA = specific_DNA > 0,
         detected_specRNA = specific_RNA > 0,
         detected_specific = specific_all > 0) %>% 
  left_join(vis_detected) %>% 
  arrange(desc(detected_visual),
          desc(specimens_visual)) %>%
  left_join(rna_detected) %>%
  left_join(dna_detected) %>% 
  select_if(!grepl('rep', names(.))) %>% 
  glimpse
 spp_detected$container_id %>% table %>% table
 
tapply(spp_detected$detected_specDNA, 
       spp_detected$common_name, sum)

# detections
detection <- spp_detected %>% 
  group_by(common_name, species) %>% 
  summarise(specd = sum(detected_specDNA),
            specr = sum(detected_specRNA),
            dna = sum(detected_dna, na.rm = T),
            rna = sum(detected_rna, na.rm = T),
             visual = sum(detected_visual, na.rm = T),
            containers = n()) %>% 
  arrange(desc(specd)) %>% 
  pivot_longer(cols = specd:visual, names_to = 'data', values_to = 'detections')

# congruence
congruence <- spp_detected %>% 
  rowwise() %>% 
  mutate(congruence_dna = (detected_dna == T & detected_specific == T),
         congruence_rna = (detected_rna == T & detected_specific == T),
         congruence_vis = (detected_visual == T & detected_specific == T)) %>% 
  group_by(common_name, species) %>% 
  summarise(specific = sum(detected_specific),
            dna = sum(congruence_dna, na.rm = T),
            rna = sum(congruence_rna, na.rm = T),
            visual = sum(congruence_vis, na.rm = T),
            containers = n()
            ) %>% 
  arrange(desc(specific)) %>% 
  pivot_longer(cols = dna:visual, names_to = 'data', values_to = 'congruence')

detection %>% left_join(congruence) %>% 
  select(-containers, species) %>% 
  mutate(prop_detection = detections/specific,
         prop_congruence = congruence/detections) %>% 
  arrange(species, prop_detection) %>% 
  relocate(common_name, species, data, prop_congruence) %>% 
  mutate_if(is.numeric, round, 2)



# venn diagram ------------------------------------------------------------

library(VennDiagram)


containersSPP <- spp_detected$container_id[spp_detected$detected_specific]
containersDNA <- spp_detected$container_id[spp_detected$detected_dna &
                                             complete.cases(spp_detected$detected_dna)] 

containersRNA <- spp_detected$container_id[spp_detected$detected_rna &
                                             complete.cases(spp_detected$detected_rna)]

containersVIS <- spp_detected$container_id[spp_detected$detected_visual &
                                             complete.cases(spp_detected$detected_visual)]

containersSPP_dna <- spp_detected$container_id[spp_detected$detected_specDNA]
containersSPP_rna <- spp_detected$container_id[spp_detected$detected_specRNA]

mycol <- c('lightblue', 'lightgreen', 'pink', 'gold')

venn.diagram(
  x = list(containersSPP_dna, containersSPP_rna,
           containersDNA, containersRNA),
  category.names = c("sppDNA" , "sppRNA", 'DNA', 'RNA'),
  filename = './figures/venn/spp_meta.png',
  output=TRUE,
  disable.logging = T,
  imagetype="png" ,
  # height = 480 , 
  # width = 480 , 
  # resolution = 300,
  col = mycol,
  fill = sapply(mycol, alpha, alpha = 0.3)
  
)



venn.diagram(
  x = list(containersSPP, unique(c(containersDNA, containersRNA)),
           containersVIS),
  category.names = c("Spp specific" , "Metabarcoding", 'Visual'),
  filename = './figures/venn/spp_meta_vis_venn_diagramm.png',
  output=TRUE,
  disable.logging = T,
  imagetype="png" ,
   height = 880 , 
   width = 880 , 
   resolution = 300,
  compression = "lzw",
  cat.fontface = 'bold',
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#FF6347'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#FF6347',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 120),
  cat.dist = c(0.065, 0.065, 0.095),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#FF6347'),
  rotation = 1
)


