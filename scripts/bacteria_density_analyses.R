#analysis of bacteria densities during coevolution

library(dplyr)
library(lme4)
library(emmeans)
library(ggplot2)
library(brms)

dat <- read.csv("data/coev_den_100223.csv",header=T)

dat <- mutate(dat,
              cfu = count * dil * dil_fact * froz,
              lcfu = log10(cfu))

dat2 <- na.omit(dat)
dat2$lcfu[dat2$lcfu == -Inf] <- 0
dat$lcfu[dat$lcfu == -Inf] <- 0

dat_means <- group_by(dat2, treat, time) %>%
  summarise(m_cfu = mean(lcfu),
            se = sd(lcfu)/sqrt(length(lcfu)))

dat_means2 <- group_by(dat2, treat) %>%
  summarise(m_cfu = mean(lcfu),
            se = sd(lcfu)/sqrt(length(lcfu)))

p1 <- ggplot() +
  geom_point(data = dat2, aes(x = time, y = lcfu, col = treat, group = interaction(rep,treat)), alpha = 0.3,size = 0.6) +
  geom_line(data = dat_means, aes(x = time, y = m_cfu, group = treat, col = treat),linewidth=1) +
  geom_line(data = dat2, aes(x = time, y = lcfu, col = treat, group = interaction(rep,treat)), alpha = 0.3) + 
  geom_errorbar(data = dat_means, aes(x = time, ymin = m_cfu-se, ymax = m_cfu+se), width=0.1) +
  geom_point(data = dat_means, aes(x = time, y = m_cfu), size = 1.5) +
  theme_bw() +
  palettetown::scale_color_poke(pokemon = "blastoise", spread = 4, labels = c("Bacteria", "Macrophages", "Bacteriophages", "Macrophages\nBacteriophages")) +
  ylab("Bacterial density\nlog10(CFU/mL)") +
  xlab("Time (days)") +
  theme(axis.text.y = element_text(size = 11, colour = "black"), axis.title = element_text(size = 12), axis.text.x = element_text(size = 11, colour = "black"), title = element_text(size = 10),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 11)) +
  labs(col = "Treatment", title = '(a) Bacteria density during coevolution') +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) 

dat$rep2 <- interaction(dat$treat, dat$rep)

dat$macrophage <- dat$treat
dat$phage <- dat$treat

dat$macrophage[dat$macrophage == "b" | dat$macrophage == "p"] <- "N"
dat$macrophage[dat$macrophage == "m" | dat$macrophage == "pm"] <- "Y"

dat$phage[dat$phage == "b" | dat$phage == "m"] <- "N"
dat$phage[dat$phage == "p" | dat$phage == "pm"] <- "Y"

m1 <- lmer(lcfu ~ time * macrophage * phage + (1|rep2), data = dat)
m2 <- lmer(lcfu ~ time + macrophage + phage + time:macrophage + time:phage + macrophage:phage + (1|rep2), data = dat)
anova(m1,m2) #sig three-way interaction
emmeans::emmeans(m1, pairwise ~ macrophage * phage | time)

#Make table of contrasts

mens <- data.frame(emmeans::emmeans(m1, pairwise ~ macrophage * phage | time)$contrasts)

mens$contrast <- as.character(mens$contrast)

mens$contrast[mens$contrast == "N N - Y N"] <- "Bacteria - Macrophages"
mens$contrast[mens$contrast == "N N - N Y"] <- "Bacteria - Bacteriophages"
mens$contrast[mens$contrast == "N N - Y Y"] <- "Bacteria - Both"
mens$contrast[mens$contrast == "Y N - N Y"] <- "Macrophages - Bacteriophages"
mens$contrast[mens$contrast == "Y N - Y Y"] <- "Macrophages - Both"
mens$contrast[mens$contrast == "N Y - Y Y"] <- "Bacteriophages - Both"

mens <- mens[,c(-4,-5)]

mens <- mutate(mens,
               estimate = round(estimate, digits = 3),
               t.ratio = round(t.ratio, digits = 3),
               p.value = round(p.value, digits = 3))

mens$p.value[mens$p.value == 0] <- "<0.001"
mens$time <- substring(mens$time, 2)

library(flextable)

men_flex <- flextable(mens) %>%
  set_header_labels(time = "Day", contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  hline(i = c(6,12,18,24,30,36), border = officer::fp_border(color="black"))%>%
  align(j = c(1:5), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  line_spacing(space = 0.75, part = "all")

macro1 <- filter(dat, treat == "m" | treat == "b")

m1 <- lmer(lcfu ~ time * macrophage + (1|rep2), data = macro1)
m2 <- lmer(lcfu ~ time + macrophage + (1|rep2), data = macro1)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ macrophage|time)

phag1 <- filter(dat, treat == "p" | treat == "b")

m1 <- lmer(lcfu ~ time * phage + (1|rep2), data = phag1)
m2 <- lmer(lcfu ~ time + phage + (1|rep2), data = phag1)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ phage|time)

mens2 <- data.frame(emmeans::emmeans(m1, pairwise ~ macrophage|time)$contrasts)

mens2 <- data.frame(emmeans::emmeans(m1, pairwise ~ phage|time)$contrasts)

mens2 <- mens2[,c(-4,-5)]
mens2 <- mutate(mens2,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

mens2$p.value[mens2$p.value == 0] <- "<0.001"
mens2$time <- substring(mens2$time, 2)

men_flex2 <- flextable(mens2) %>%
  set_header_labels(time = "Day", contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:5), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')
