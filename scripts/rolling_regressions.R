### analysis of growth rates of clones taken from different treatments of macrophage project. All independent.

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)
library(zoo) #install.packages(zoo)
library(broom) #install.packages(broom)
library(growthcurver) # install.packages(growthcurver)
library(nls.multstart)

#plate 3 - BMPS wells C8, E8 and F8. D8 is duplicate of C8. 

#C3, D3, E3 = Bs
#C4, D4, E4 = BMs
#C5, D5, E5 = BPRs
#C6, D6, E6 = BPSs
#C7, D7, E7 = BMPRs
#C8, D8, E8 = BMPSs
#C9, D9, E9 = As
#C10, D10, E10 = Blanks

p1 <- read.csv('data/growth_run1_091224.csv',header=T)
p2 <- read.csv('data/growth_run2_091224.csv',header=T)
p3 <- read.csv('data/growth_run3_091224.csv',header=T)

p1$run <- 1
p2$run <- 2
p3$run <- 3

p12 <- merge(p1, p2, all = T)
p12 <- p12 %>%
  separate(well, c('row','col'))
p12 <- filter(p12, row == 'C' | row == 'D' | row == 'E')

p12$col <- as.numeric(p12$col)

p12 <- filter(p12, col >= 3)
p12 <- filter(p12, col <= 10)

p12$treat <- NA
p12$treat[p12$col == 3] <- 'B'
p12$treat[p12$col == 4] <- 'BM'
p12$treat[p12$col == 5] <- 'BPR'
p12$treat[p12$col == 6] <- 'BPS'
p12$treat[p12$col == 7] <- 'BMPR'
p12$treat[p12$col == 8] <- 'BMPS'
p12$treat[p12$col == 9] <- 'A'
p12$treat[p12$col == 10] <- 'blank'

#plate 3 - BMPS wells C8, E8 and F8. D8 is duplicate of C8. 

p3 <- p3 %>%
  separate(well, c('row','col'))

p3_f8 <- filter(p3, row == 'F' & col == 8)

p3 <- filter(p3, row == 'C' | row == 'D' | row == 'E')
p3 <- p3[-c(1160:1220),]
p3 <- merge(p3, p3_f8, all = T)
p3$col <- as.numeric(p3$col)
p3 <- filter(p3, col >= 3)
p3 <- filter(p3, col <= 10)

p3$treat <- NA
p3$treat[p3$col == 3] <- 'B'
p3$treat[p3$col == 4] <- 'BM'
p3$treat[p3$col == 5] <- 'BPR'
p3$treat[p3$col == 6] <- 'BPS'
p3$treat[p3$col == 7] <- 'BMPR'
p3$treat[p3$col == 8] <- 'BMPS'
p3$treat[p3$col == 9] <- 'A'
p3$treat[p3$col == 10] <- 'blank'

p123 <- merge(p12, p3, all = T)

p123$row[p123$row == 'F'] <- 'D'

p123$rep <- interaction(p123$treat, p123$row)

p123 <- mutate(p123,
              time_hr = time / 3600)

p123$od2 <- log(p123$od)

#visualise slopes

ggplot(p123, aes(x = time_hr, y = od2)) +
  geom_point(aes(col = run)) +
  facet_wrap(~treat)

p123 <- filter(p123, ! treat == 'blank')

ggplot(p123, aes(x = time_hr, y = od2)) +
  geom_point(aes(col = run)) +
  facet_wrap(~treat+run) 

bps <- filter(p123, treat == 'BPS') 

ggplot(bps, aes(x = time_hr, y = od2)) +
  geom_point(aes(col = run)) +
  facet_wrap(~run+rep) 

bpr <- filter(p123, treat == 'BPR') 

ggplot(bpr, aes(x = time_hr, y = od2)) +
  geom_point(aes(col = run)) +
  facet_wrap(~run+rep) 

#after 5 hrs the plots start to go a bit mad 

p123_2 <- p123

p123 <- filter(p123, time_hr < 5)

#rolling regression code

roll_regress <- function(x){
  temp <- data.frame(x)
  mod <- lm(temp)
  temp <- data.frame(slope = coef(mod)[[2]],
                     slope_lwr = confint(mod)[2, ][[1]],
                     slope_upr = confint(mod)[2, ][[2]],
                     intercept = coef(mod)[[1]],
                     rsq = summary(mod)$r.squared, stringsAsFactors = FALSE)
  return(temp)
}

num_points = ceiling(0.65*60/(60*0.167)) 

models <- p123 %>%
  group_by(rep, run, treat) %>%
  do(cbind(model = select(., od2, time_hr) %>% 
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = select(., time_hr),
           ln_od = select(., od2))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')

preds <- models %>%
  filter(., !is.na(slope)) %>%
  group_by(time_hr) %>%
  do(data.frame(time2 = c(.$time_hr - 2, .$time_hr + 2))) %>%
  left_join(., models) %>%
  mutate(pred = (slope*time2) + intercept)

preds <- models %>%
  filter(., !is.na(slope)) %>%
  left_join(., models) 

growth_rate <- filter(models, slope == max(slope, na.rm = TRUE))

grows <- group_by(growth_rate, rep, treat) %>%
  summarise(m.grow = mean(slope),
            slope_lwr = mean(slope_lwr),
            slope_upr = mean(slope_upr))

grows <- read.csv('est_growth_rates.csv', header = T)

b <- ggplot(grows, aes(x = treat2, y = m.grow, group = interaction(treat, res))) +
  MicrobioUoE::geom_pretty_boxplot(aes(col =res, fill =res)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8, colour = 'black'), axis.title = element_text(size = 10, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 10))  +
  ylab('Growth rate (hr-1)') +
  xlab('Evolution treatment') +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 2, labels = c('Susceptible', 'Resistant')) +
  palettetown::scale_color_poke(pokemon = 'articuno', spread = 2, labels = c('Susceptible', 'Resistant')) +
  labs(col = 'Phage resistance', fill = 'Phage resistance') +
  scale_x_discrete(labels = c('Ancestor', 'No phage\nNo macrophage', 'No phage\nMacrophage', 'Phage\nMacrophage', 'Phage\nNo macrophage')) +
  labs(title = '(b)')

m1 <- lm(m.grow ~ treat2 * res, data = grows)
m2 <- lm(m.grow ~ treat2 + res, data = grows)
anova(m1,m2)
m3 <- lm(m.grow ~ res, data = grows)
anova(m2,m3)
m4 <- lm(m.grow ~ treat2, data = grows)
anova(m2,m4)
m5 <- lm(m.grow ~ 1, data = grows)
anova(m3,m5)
