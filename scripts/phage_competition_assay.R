library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(flextable)

pcom <- read.csv('data/macro_phagecomp.csv',header=T)
#SM = spent media from 24 hr macrophage growth. AM = activated macrophage media. M = macrophages. B = just bacteria in fresh media with phage
#done Dec 2024 

pcom2 <- mutate(pcom, 
                m141 = log(den_14/160),
                mpnm = log(den_pnm/880),
                rfit = m141/mpnm)

ggplot(pcom2, aes(x = treat, y = rfit)) +
  MicrobioUoE::geom_pretty_boxplot(col ='black', fill = 'black') +
  geom_point(shape = 21, fill = 'white', size = 4) + 
  theme_bw() +
  theme(axis.text = element_text(size = 13, colour = 'black'), axis.title = element_text(size = 15, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 15, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  ylab('Relative growth rate of 14-1 to PNM') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Activated\nmacrophage media', 'Control\n(fresh media\nno macrophages)', 'Macrophages', 'Spent\nmacrophage media'))

m1 <- lm(rfit ~ treat, data = pcom2)
plot(m1)
m2 <- lm(rfit ~ 1, data = pcom2)
anova(m1,m2,test='F')

p141 <- pcom2[,c(1,9)]
p141$phage <- 'p141'
names(p141)[2] <- 'm'
p141$rep = rep(1:5, 4)
p141$rep2 = interaction(p141$treat, p141$rep)

pnm <- pcom2[,c(1,10)]
pnm$phage <- 'pnm'
names(pnm)[2] <- 'm'
pnm$rep = rep(1:5, 4)
pnm$rep2 = interaction(pnm$treat, pnm$rep)

phage_m <- merge(p141,pnm,all=T)

ggplot(phage_m, aes(x = treat, y = m)) +
  MicrobioUoE::geom_pretty_boxplot(col ='black', fill = 'black') +
  geom_point(shape = 21, aes(fill = phage), size = 3, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 13, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 13, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 13), legend.title = element_text(size = 13)) +
  ylab('Bacteriophage growth rate (m)') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Activated\nmacrophage media', 'Control\n(fresh media,\nno macrophages)', 'Macrophages', 'Spent\nmacrophage media')) +
  palettetown::scale_fill_poke(pokemon = 'wigglytuff', spread = 4, labels = c('14-1', 'PNM')) +
  labs(fill = 'Bacteriophage')

m1 <- lmer(m ~ phage * treat + (1|rep2), data = phage_m)
m2 <- lmer(m ~ phage + treat + (1|rep2), data = phage_m)
anova(m1,m2)
m3 <- lmer(m ~ phage + (1|rep2), data = phage_m)
anova(m2,m3) #sig effect of treat
m4 <- lmer(m ~ treat + (1|rep2), data = phage_m)
anova(m2,m4) #sig effect of phage
emmeans::emmeans(m2,pairwise ~ phage)
emmeans::emmeans(m2,pairwise ~ treat) 

phage_m2 <- phage_m
phage_m2$treat[phage_m2$treat == 'am'] <- 'Activated Media'
phage_m2$treat[phage_m2$treat == 'b'] <- ' No Macrophage Control'
phage_m2$treat[phage_m2$treat == 'm'] <- 'Macrophages'
phage_m2$treat[phage_m2$treat == 'sm'] <- 'Spent Macrophage Media'

m2 <- lmer(m ~ phage + treat + (1|rep2), data = phage_m2)

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
