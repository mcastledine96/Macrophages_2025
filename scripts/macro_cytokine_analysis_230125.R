library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)

#standard curve data
sc <- read.csv('data/tnf_alpha_SC_220125.csv', header = T)
#cytokine measurements
tnf <- read.csv('data/cytokine_tnfalpha_220125.csv', header = T)

#calculate od difference

sc <- mutate(sc,
             od = od_450 - od_620)

tnf <- mutate(tnf,
             od = od_450 - od_620)

#calculate average for duplicates and check values are within 20% of this value

sc_m <- group_by(sc, standard, cyto) %>%
  summarise(m_od = mean(od))

sc <- sc %>%
  left_join(sc_m, by = 'standard')

sc <- mutate(sc,
             st = m_od*0.2)

sc$InRange <- with(sc, od >= (m_od - st) & od <= (m_od + st)) 

tnf$rep2 <- interaction(tnf$treat, tnf$brep)

tnf_m <- group_by(tnf, brep, treat, phager, macro_e, rep2) %>%
  summarise(m_od = mean(od))

tnf <- tnf %>%
  left_join(tnf_m, by = 'rep2')

tnf <- mutate(tnf,
             st = m_od*0.2)

tnf$InRange <- with(tnf, od >= (m_od - st) & od <= (m_od + st))

#all duplicates are within the accepted range of the mean - specified in the ELISA booklet

#sc data

sc_m <- filter(sc_m, ! standard == 'blank')

ggplot(sc_m, aes(x = cyto, y = m_od)) +
  geom_point() 

m1 <- lm(cyto ~ m_od + I(m_od^2), data = sc_m)
m2 <- lm(cyto ~ m_od, data = sc_m)
anova(m1,m2) #sig quadratic effect

library(ciTools)

newdata <- with(sc_m, expand_grid(m_od = seq(min(sc_m$m_od),max(sc_m$m_od),length.out=100)))
new <- add_ci(newdata, m1, alpha = 0.05, type = "boot", includeRanef = FALSE, nSims = 2000)

ggplot() +
  geom_point(data = sc_m, aes(x = m_od, y = cyto)) +
  geom_line(data = new, aes(x = m_od, y = pred))

#predict cytokine levels 
n <- data.frame(m_od = tnf_m$m_od)
tnf_m$cyto <- predict(m1, newdata = n)

tnf_t <- filter(tnf_m, ! treat == 'macro')
m <- filter(tnf_m, treat == 'macro')

tnf_t <- mutate(tnf_t,
                cyto = cyto * 100)

tnf_m <- merge(tnf_t, m, all = T)

tnf_m2 <- group_by(tnf_m, treat, phager, macro_e) %>%
  summarise(m = mean(cyto),
            se = sd(cyto)/sqrt(3))

#first see if any difference between treatments

tnf_m$treat2 <- tnf_m$treat
tnf_m$treat2[tnf_m$treat2 == 'bmpr'] <- 'bmp'
tnf_m$treat2[tnf_m$treat2 == 'bmps'] <- 'bmp'
tnf_m$treat2[tnf_m$treat2 == 'bpr'] <- 'bp'
tnf_m$treat2[tnf_m$treat2 == 'bps'] <- 'bp'

tnf_m3 <- group_by(tnf_m, treat2) %>%
  summarise(m = mean(cyto),
            se = sd(cyto)/sqrt(3))

b <- ggplot() +
  geom_bar(data = tnf_m3, aes(x = treat2, y = m), stat = "identity", fill = 'white', col = 'black') +
  theme_bw() +
  geom_errorbar(data = tnf_m3,  aes(x = treat2, ymin = m-se, ymax = m+se), width = 0.3)  +
  geom_point(data = tnf_m, aes(x = treat2, y = cyto, fill = phager), shape = 21, size = 2, position = position_jitter(0.1)) + 
  ylab('TNF-alpha (pg/mL)') +
  xlab('Evolution Treatment') +
  scale_x_discrete(labels = c('Ancestor', 'No phage\nNo macrophage', 'No phage\nMacrophage', 'Phage\nMacrophage', 'Phage\nNo macrophage', 'Macrophage\ncontrol')) +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 4, labels = c('Susceptible', 'Resistant', 'Not relevant')) +
  labs(fill = 'Phage resistance') +
  theme(axis.text = element_text(size = 8, colour = 'black'), axis.title = element_text(size = 10, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  labs(title = '(b)')

m1 <- lm(cyto ~ treat2, data = tnf_m)
m2 <- lm(cyto ~ 1, data = tnf_m)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ treat2)

m1 <- lm(log(cyto) ~ treat2, data = tnf_m)
m2 <- lm(log(cyto) ~ 1, data = tnf_m)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ treat2) #only difference is to the macrophage treatment

tnf_m2 <- tnf_m
tnf_m2$treat2[tnf_m2$treat2 == 'anc'] <- "Ancestor"
tnf_m2$treat2[tnf_m2$treat2 == 'b'] <- "Bacteria-only"
tnf_m2$treat2[tnf_m2$treat2 == 'bm'] <- "Macrophage"
tnf_m2$treat2[tnf_m2$treat2 == 'bmp'] <- "Bacteriophage+Macrophage"
tnf_m2$treat2[tnf_m2$treat2 == 'bp'] <- "Bacteriophage"
tnf_m2$treat2[tnf_m2$treat2 == 'macro'] <- "Macrophage Control"

m1 <- lm(log(cyto) ~ treat2, data = tnf_m2)

mens2 <- data.frame(emmeans::emmeans(m1, pairwise ~ treat2)$contrasts)
mens2 <- mens2[,c(-3,-4)]
mens2 <- mutate(mens2,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

men_flex2 <- flextable(mens2) %>%
  set_header_labels(contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')


#then look at resistance

tnf_ta <- filter(tnf_t, ! treat == 'anc')
tnf_ta <- filter(tnf_ta, ! treat == 'b')
tnf_ta <- filter(tnf_ta, ! treat == 'bm')

m1 <- lm(log(cyto) ~ phager * macro_e, data = tnf_ta)
m2 <- lm(log(cyto) ~ phager + macro_e, data = tnf_ta)
anova(m1,m2,test='F')
m3 <- lm(log(cyto) ~ macro_e, data = tnf_ta) 
anova(m2,m3)
m4 <- lm(log(cyto) ~ phager, data = tnf_ta) #phager least sig so drop
anova(m2,m4)
m5 <- lm(log(cyto) ~ 1, data = tnf_ta)
anova(m3,m5)
