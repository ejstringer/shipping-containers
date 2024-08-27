
cqcopy <- read.csv('./data/Khapra beetle_2023_FS.csv')

cqcopy %>% names
test <- cqcopy[c(3,grep('Furui_DNA', names(cqcopy)))]

test %>% 
  filter(complete.cases(Furui_DNA_R1)) %>% 
  unique() %>% 
  mutate(across(Furui_DNA_R1:Copy_Furui_DNA_R3, ~ as.numeric(.x))) %>% 
  pivot_longer(cols = Furui_DNA_R1:Copy_Furui_DNA_R3) %>%
  separate(col = name, sep = 'DNA', into = c('var', 'rep'))%>% 
  mutate_if(is.character, ~ gsub('_', '', .)) %>% 
  mutate_if(is.character, ~ gsub('R', '', .)) %>% 
  filter(value > 0) %>% 
  pivot_wider(values_from = value, names_from = var) -> furui
  ggplot(furui, aes(Furui, CopyFurui, colour = rep))+
  geom_smooth(alpha = 0.1)+
  geom_point(size = 2)+
  theme_light()+
 # facet_grid(~rep)+
  scale_y_log10()+
  xlab('cq')+
  ylab('copy number')
  
  p <- ggplot(furui, aes(CopyFurui, Furui))+
    geom_point(size = 5, alpha = 0.1)+
    geom_smooth(alpha = 0.1, colour = 'red')+
    theme_classic()+
    # facet_grid(~rep)+
    scale_x_log10()+
    ggtitle('Furui')+
    ylab('cq')+
    xlab('copy number');p

  lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 5),
                          b = format(unname(coef(m)[2]), digits = 4),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  
  df <- data.frame(x = log(furui$CopyFurui), y = furui$Furui)

  
  p1 <- p + geom_text(x = 1, y = 25, label = lm_eqn(df), parse = TRUE)
  p1 
  
  lm(CopyFurui ~ Furui, furui) %>% summary
  lm(log(CopyFurui) ~ Furui, furui) %>% summary
  m<- lm(Furui ~ log(CopyFurui), furui)
  m %>% summary
  m
  plot(log(furui$CopyFurui), furui$Furui)
    abline(m, col = 'red')
  
  test <- cqcopy[c(3,grep('Olson', names(cqcopy)))]
  
  
  test %>% 
    pivot_longer(cols = Olson_R1:Copy_Olson_R3) %>%
    separate(col = name, sep = 'son', into = c('var', 'rep'))%>% 
    mutate_if(is.character, ~ gsub('_', '', .)) %>% 
    mutate_if(is.character, ~ gsub('R', '', .)) %>% 
    filter(value > 0) %>% 
    pivot_wider(values_from = value, names_from = var) -> olson
  
  
    p<-ggplot(olson,aes(CopyOl, Ol))+
      geom_point(size = 5, alpha = 0.1)+
      geom_smooth(alpha = 0.1, colour = 'red')+
      theme_classic()+
      # facet_grid(~rep)+
      scale_x_log10()+
      ggtitle('Olson')+
      ylab('cq')+
      xlab('copy number');p

    
    
    log(m$coefficients)
    df <- data.frame(x = log(olson$CopyOl), y = olson$Ol)
    p2 <- p + geom_text(x = 2, y = 30, label = lm_eqn(df), parse = TRUE)
    
    pg <- grid.arrange(p1, p2, ncol = 2)
    
    ggsave('./figures/copy_number_formula_khapra.png', pg)
    