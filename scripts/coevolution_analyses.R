library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)

dat2 <- read.csv("data/coev_res_2.csv", header = T)

dat_means <- group_by(dat2, treat, time, name) %>%
  summarise(prop2 = mean(prop),
            se = sd(prop)/sqrt(length(prop)))

dat_means14 <- group_by(dat2, name) %>%
  summarise(prop2 = mean(prop),
            se = sd(prop)/sqrt(length(prop)))

#dat2_x <- filter(dat2, ! time == "t6")

#dat_means142 <- group_by(dat2_x, name) %>% summarise(prop2 = mean(prop), se = sd(prop)/sqrt(length(prop)))

dat2_noa <- filter(dat2, ! name == "X14")
dat2_noa <- filter(dat2_noa, ! name == "pnm")
dat_means_noa <- filter(dat_means, ! name == "pnm")
dat_means_noa <- filter(dat_means_noa, ! name == "X14")

dat_means_t <- group_by(dat2_noa, treat, time) %>%
  summarise(prop2 = mean(prop),
            se = sd(prop)/sqrt(length(prop)))

fac_labs <- c("(a) Macrophages present", "(b) Macrophages absent")
names(fac_labs) <- c("macro", "phage")

ggplot() +
  geom_point(data = dat2_noa, aes(x = time, y = prop, col = name, group = name), size = 1, alpha = 0.5, position = position_dodge(0.4)) +
  geom_point(data = dat_means_noa, aes(x = time, y = prop2, col = name, group = name), size = 2, position = position_dodge(0.4)) +
  facet_wrap(~treat, labeller = labeller(treat = fac_labs)) +
  geom_errorbar(data = dat_means_noa, aes(x = time, ymin = prop2-se, ymax = prop2+se, group = name, col = name), position = position_dodge(0.4), width = 0.1) + 
  theme_bw() +
  palettetown::scale_color_poke(pokemon = "nidoranm", spread = 3, labels = c("2", "4", "6")) +
  theme(axis.text = element_text(size = 13, colour = 'black'), axis.title = element_text(size = 15, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 15, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  labs(col = "Bacteriophage time point") +
  xlab("Bacteria time point") +
  scale_x_discrete(labels = c("2", "4", "6")) +
  ylab("Proportion bacteriophage resistant")

dat2a <- dat2[nchar(dat2$name) > 2, ]
datma <- dat_means[nchar(dat_means$name) > 2, ]

#plot of resistance to PNM and 14-1 at day 6
ggplot() +
  geom_point(data = dat2a, aes(x = name, y = prop, col = treat, group = treat), size = 1, alpha = 0.5, position = position_dodge(0.4)) +
  geom_point(data = datma, aes(x = name, y = prop2, col = treat, group = treat), size = 2, position = position_dodge(0.4)) +
  geom_errorbar(data = datma, aes(x = name, ymin = prop2-se, ymax = prop2+se, col = treat, group = treat), position = position_dodge(0.4), width = 0.1) + 
  theme_bw() +
  palettetown::scale_color_poke(pokemon = "nidoranm", spread = 3, labels = c('Macrophages present', 'Macrophages absent')) +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  xlab("Phage") +
  ylab("Proportion bacteriophage resistant") +
  scale_x_discrete(labels = c('PNM', '14-1')) +
  labs(col = 'Treatment')

#analysis of resistance to 14-1
m1 <- glm(cbind(res,sus) ~ treat, data = dat2a14, family = quasibinomial)
plot(m1) 
m2 <- glm(cbind(res,sus) ~ 1, data = dat2a14, family = quasibinomial)
anova(m1,m2, test = 'F')

#analysis of coevolution

dat2$rep2 <- interaction(dat2$treat, dat2$rep)

dat4 <- filter(dat2, ! name == "X14")
dat4 <- filter(dat4, ! name == "pnm")

m_dat <- filter(dat4, treat == "macro")
p_dat <- filter(dat4, treat == "phage")

p_dat2 <- group_by(p_dat, treat, rep, name) %>%
  summarise(m_prop = mean(prop),
            s_t = sum(total),
            s_s = sum(sus),
            s_r = sum(res))

m1 <- glmer(cbind(res,sus) ~ time * name + (1|rep2), family = binomial, data = m_dat) 
m2 <- glmer(cbind(res,sus) ~ time + name + (1|rep2), family = binomial, data = m_dat)
anova(m1,m2) # 0.9559
m3 <- glmer(cbind(res,sus) ~ name + (1|rep2), family = binomial, data = m_dat)
anova(m2,m3) #time = 2.933e-11 ***
m4 <- glmer(cbind(res,sus) ~ time + (1|rep2), family = binomial, data = m_dat)
anova(m2,m4) #name 0.2942
emmeans::emmeans(m4, pairwise ~ time) #increasing resistance through time between t2 and 4 and t2 and t6 to phages as a whole

m1 <- glmer(cbind(res,sus) ~ time * name + (1|rep2), family = binomial, data = p_dat)
m2 <- glmer(cbind(res,sus) ~ time + name + (1|rep2) , family = binomial, data = p_dat)
anova(m1,m2)
m3 <- glmer(cbind(res,sus) ~ name + (1|rep2), family = binomial, data = p_dat)
anova(m2,m3) #time = 0.09842 
m4 <- glmer(cbind(res,sus) ~ time + (1|rep2), family = binomial, data = p_dat)
anova(m2,m4) #name 0.007595

emmeans::emmeans(m3, pairwise ~ name) #bacteria more susceptible generally to phage from t6 then t2 

#look at contemporary resistance levels for both treatments

dat_con <- filter(dat4, time == "t2" & name == "X2" | time == "t4" & name == "X4"| time == "t6" & name == "X6")

m1 <- glmer(cbind(res,sus) ~ treat * time + (1|rep2), family = binomial, data = dat_con)
m2 <- glmer(cbind(res,sus) ~ treat + time + (1|rep2), family = binomial, data = dat_con) 
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ treat|time)

