#analysis of different clones biofilm production measured using resazurin assay. All clones are from independent replicates. 3 technical replicates of assay. Assays conducted without phage.

library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)

biof <- read.csv('data/macro_biofilm_291124.csv',header=T)

ggplot(biof, aes(x = treat, y = lod, group = interaction(treat, res), col = res)) +
  geom_boxplot() +
  geom_point(position = position_dodge(0.75)) +
  facet_wrap(~run)

biof$lod <- log(biof$od)

biof2 <- group_by(biof, treat, res, rep) %>%
  summarise(odm = mean(od),
            lod = log(odm))

a <- ggplot(biof2, aes(x = treat, y = lod, group = interaction(treat, res))) +
  MicrobioUoE::geom_pretty_boxplot(aes(col =res, fill =res)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size = 2) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8, colour = 'black'), axis.title = element_text(size = 10, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  ylab('Biofilm production\nlog(fluroescence 560/590nm)') +
  scale_x_discrete(labels = c('Ancestor', 'No phage\nNo macrophage', 'No phage\nMacrophage', 'Phage\nMacrophage', 'Phage\nNo macrophage')) +
  xlab('Evolution Treatment')+
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 2, labels = c('Susceptible', 'Resistant')) +
  palettetown::scale_color_poke(pokemon = 'articuno', spread = 2, labels = c('Susceptible', 'Resistant')) +
  labs(col = 'Phage resistance', fill = 'Phage resistance') +
  labs(title = '(a)')

(a + b) / c + patchwork::plot_layout()
  
a + b + patchwork::plot_layout(nrow = 2)

m1 <- lm(lod ~ treat * res, data = biof2)
m2 <- lm(lod ~ treat + res, data = biof2)
anova(m1,m2,test='F')
m3 <- lm(lod ~ res, data = biof2)
anova(m2,m3,test='F')
m4 <- lm(lod ~ treat, data = biof2)
anova(m2,m4,test='F')
emmeans(m2, pairwise ~ res)
emmeans(m2, pairwise ~ treat)
#so looks like phage resistance reduces biofilm formation. And evolving with macrophages and phages reduces biofilm formation more than evolving without macrophages but with phages. Suggests some incraese in evo with phage and some decrease in evo without phage. Will see if this is linked to growth rates. Big variance in bacteria-macrophage data. 

biof3 <- biof2

biof3$treat[biof3$treat == 'A'] <- 'Ancestor'
biof3$treat[biof3$treat == 'B'] <- 'Bacteria-only'
biof3$treat[biof3$treat == 'BM'] <- 'Macrophage'
biof3$treat[biof3$treat == 'BP'] <- 'Bacteriophage'
biof3$treat[biof3$treat == 'BMP'] <- 'Bacteriophage+Macrophage'

m2 <- lm(lod ~ treat + res, data = biof3)

mens2 <- data.frame(emmeans::emmeans(m2, pairwise ~ treat)$contrasts)
mens2 <- mens2[,c(-3,-4)]
mens2 <- mutate(mens2,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

mens2$p.value[mens2$p.value == 0] <- "<0.001"

men_flex2 <- flextable(mens2) %>%
  set_header_labels(contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#do growth rate assays next to check for effect of different growth rates

#write.csv(biof2, 'est_biofilm.csv', row.names = F)

#now added growth rates

b_g <- read.csv('est_biofilm.csv', header = T)

ggplot(b_g, aes(x = grate, y = lod)) +
  geom_point(aes(fill = interaction(treat, res)), shape = 21, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 13, colour = 'black'), axis.title = element_text(size = 15, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 15, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  ylab('Biofilm production log(fluroescence 560/590nm)') +
  xlab('Growth rate hr-1') +
  scale_fill_manual(values = c("#a06aa0", '#db596b', 'pink', '#e79652', '#85cdc4', '#31646f', 'lightgreen'), labels = c('A', 'M', 'P-N', 'P-Y', 'B', 'MP-N', 'MP-Y'))  +
  scale_color_manual(values = c("#a06aa0", '#db596b', 'pink', '#e79652', '#85cdc4', '#31646f', 'lightgreen'), labels = c('A', 'M', 'P-N', 'P-Y', 'B', 'MP-N', 'MP-Y')) +
  labs(color = 'Clone', fill = 'Clone')


m1 <- lm(lod ~ grate, data = b_g)
m2 <- lm(lod ~ 1, data = b_g)
anova(m1,m2) #no correlation between growth rate and biofilm formation
