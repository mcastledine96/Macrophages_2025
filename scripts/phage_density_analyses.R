library(dplyr)
library(lme4)
library(emmeans)
library(ggplot2)
library(brms)
library(Rmpfr)
library(tidyverse)

phag <- read.csv("data/phage_den_080423.csv",header=T)

phag <- na.omit(phag)

phag2 <- group_by(phag, macro, time) %>%
  summarise(m_den = mean(lden),
            se = sd(lden)/sqrt(length(lden)),
            m_den2 = mean(den),
            se2 = sd(den)/sqrt(length(den)))
p2 <- ggplot() +
  geom_point(data = phag, aes(x = time, y = lden, col = macro, group = interaction(rep,macro)), alpha = 0.3,size = 0.6) +
  geom_line(data = phag2, aes(x = time, y = m_den, group = macro, col = macro),linewidth=1) +
  geom_line(data = phag, aes(x = time, y = lden, col = macro, group = interaction(rep,macro)), alpha = 0.3) + 
  geom_errorbar(data = phag2, aes(x = time, ymin = m_den-se, ymax = m_den+se), width=0.1) +
  geom_point(data = phag2, aes(x = time, y = m_den), size = 1.5) +
  theme_bw() +
  ylab("Bacteriophage density\nlog10(PFU/mL)") +
  xlab("Time (days)") +
  theme(axis.text.y = element_text(size = 11, colour = "black"), axis.title = element_text(size = 12), axis.text.x = element_text(size = 11, colour = "black"),legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 11), title = element_text(size = 10)) +
  scale_x_discrete(labels = c("1","2","3","4","5","6")) +
  scale_color_manual(values = c("#484880","#e8a858"), labels = c("Absent", "Present")) +
  scale_fill_manual(values = c("#484880","#e8a858"), labels = c("Absent", "Present")) +
  labs(col = "Macrophages", title = '(b) Bacteriophage density during coevolution') 
phag$rep2 <- interaction(phag$rep, phag$macro)

pm1 <- lmer(lden ~ macro * time + (1|rep2), data = phag)
pm2 <- lmer(lden ~ macro + time + (1|rep2), data = phag)
anova(pm1,pm2)
pm3 <- lmer(lden ~ time + (1|rep2), data = phag)
anova(pm2,pm3)
pm4 <- lmer(lden ~ macro + (1|rep2), data = phag)
anova(pm2,pm4) #sig effect of time
pm5 <- lmer(lden ~ 1 + (1|rep2), data = phag)
anova(pm3,pm5)
emmeans::emmeans(pm3, pairwise ~ time) #only sig comparisons are to T6

mens2 <- data.frame(emmeans::emmeans(pm3, pairwise ~ time)$contrasts)

mens2 <- mens2[,c(-3,-4)]
mens2 <- mutate(mens2,
                estimate = round(estimate, digits = 3),
                t.ratio = round(t.ratio, digits = 3),
                p.value = round(p.value, digits = 3))

mens2$p.value[mens2$p.value == 0] <- "<0.001"
mens2$contrast <- gsub("T", "", mens2$contrast)

men_flex2 <- flextable(mens2) %>%
  set_header_labels(contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  line_spacing(space = 0.75, part = "all")
