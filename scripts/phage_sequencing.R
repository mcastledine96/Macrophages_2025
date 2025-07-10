library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(vcfR)
library(vegan)

seq_m <- read.csv('data/macro_phage_totals.csv', header =T)

seq_m <- mutate(seq_m, prop = n_reads/n_tot)

ggplot(seq_m, aes(x = time, y = prop)) +
  geom_point() +
  facet_wrap(~macrophage+phage)

pnm_dat <- filter(seq_m, phage == 'pnm')

pnm_dat_p <- filter(pnm_dat, macrophage == 'N')
#write.csv(pnm_dat_p, "data/pnm_propres.csv", row.names = F)

pnm_dat_p <- read.csv("data/pnm_propres.csv", header = T)

y <- log(pnm_dat_p$prop/(1-pnm_dat_p$prop))
m1 <- lmer(asin(sqrt(prop)) ~ time + (1|rep), data = pnm_dat_p)
m2 <- lmer(asin(sqrt(prop)) ~ 1 + (1|rep), data = pnm_dat_p)
anova(m1,m2)

y <- cbind(pnm_dat_p$s_r, pnm_dat_p$s_s)
m1 <- glmer(y ~ prop + (1|rep), data = pnm_dat_p, family = binomial)
m2 <- glmer(y ~ 1 + (1|rep), data = pnm_dat_p, family = binomial)
anova(m1,m2)

pnm_dat2 <- filter(pnm_dat, time == 'T6')

pnm_datm <- group_by(pnm_dat2, macrophage) %>%
  summarise(prop_s = sum(p_a)/6,
            sam = sum(p_a))

ggplot(pnm_datm, aes(x = macrophage, y = sam)) +
  theme_bw() +
  geom_bar(stat = 'identity')  +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 14, colour = 'black')) +
  ylab("Number of samples containing PNM") +
  xlab("Macrophages") +
  scale_x_discrete(labels = c("Absent", "Present"))

### variant analyses

p1 <- read_tsv("seq_vars/1-BMP1T2_14-1.tsv")
p1$Sample <- 'BMP1T2'
p2 <- read_tsv("seq_vars/2-BMP1T6_14-1.tsv")
p2$Sample <- 'BMP1T6'
p3 <- read_tsv("seq_vars/3-BMP2T2_14-1.tsv")
p3$Sample <- 'BMP2T2'
p4 <- read_tsv("seq_vars/4-BMP2T6_14-1.tsv")
p4$Sample <- 'BMP2T6'
p5 <- read_tsv("seq_vars/5-BMP3T2_14-1.tsv")
p5$Sample <- 'BMP3T2'
p6 <- read_tsv("seq_vars/6-BMP3T6_14-1.tsv")
p6$Sample <- 'BMP3T6'
p7 <- read_tsv("seq_vars/7-BMP4T2_14-1.tsv")
p7$Sample <- 'BMP4T2'
p8 <- read_tsv("seq_vars/8-BMP4T6_14-1.tsv")
p8$Sample <- 'BMP4T6'
p9 <- read_tsv("seq_vars/9-BMP5T2_14-1.tsv")
p9$Sample <- 'BMP5T2'
p10 <- read_tsv("seq_vars/10-BMP5T6_14-1.tsv")
p10$Sample <- 'BMP5T6'
p11 <- read_tsv("seq_vars/11-BMP6T2_14-1.tsv")
p11$Sample <- 'BMP6T2'
p12 <- read_tsv("seq_vars/12-BMP6T6_14-1.tsv")
p12$Sample <- 'BMP6T6'

p13 <- read_tsv("seq_vars/13-BP1T2_14-1.tsv")
p13$Sample <- 'BP1T2'
p14 <- read_tsv("seq_vars/14-BP1T6_14-1.tsv")
p14$Sample <- 'BP1T6'
p15 <- read_tsv("seq_vars/15-BP2T2_14-1.tsv")
p15$Sample <- 'BP2T2'
p16 <- read_tsv("seq_vars/16-BP2T6_14-1.tsv")
p16$Sample <- 'BP2T6'
p17 <- read_tsv("seq_vars/17-BP3T2_14-1.tsv")
p17$Sample <- 'BP3T2'
p18 <- read_tsv("seq_vars/18-BP3T6_14-1.tsv")
p18$Sample <- 'BP3T6'
p19 <- read_tsv("seq_vars/19-BP4T2_14-1.tsv")
p19$Sample <- 'BP4T2'
p20 <- read_tsv("seq_vars/20-BP4T6_14-1.tsv")
p20$Sample <- 'BP4T6'
p21 <- read_tsv("seq_vars/21-BP5T2_14-1.tsv")
p21$Sample <- 'BP5T2'
p22 <- read_tsv("seq_vars/22-BP5T6_14-1.tsv")
p22$Sample <- 'BP5T6'
p23 <- read_tsv("seq_vars/23-BP6T2_14-1.tsv")
p23$Sample <- 'BP6T2'
p24 <- read_tsv("seq_vars/24-BP6T6_14-1.tsv")
p24$Sample <- 'BP6T6'

m_dat <- merge(p1,p2,all=T)
m_dat <- merge(m_dat,p3,all=T)
m_dat <- merge(m_dat,p4,all=T)
m_dat <- merge(m_dat,p5,all=T)
m_dat <- merge(m_dat,p6,all=T)
m_dat <- merge(m_dat,p7,all=T)
m_dat <- merge(m_dat,p8,all=T)
m_dat <- merge(m_dat,p9,all=T)
m_dat <- merge(m_dat,p10,all=T)
m_dat <- merge(m_dat,p11,all=T)
m_dat <- merge(m_dat,p12,all=T)

m_dat$macro <- "Y"

p_dat <- merge(p13,p14,all=T)
p_dat <- merge(p_dat,p15,all=T)
p_dat <- merge(p_dat,p16,all=T)
p_dat <- merge(p_dat,p17,all=T)
p_dat <- merge(p_dat,p18,all=T)
p_dat <- merge(p_dat,p19,all=T)
p_dat <- merge(p_dat,p20,all=T)
p_dat <- merge(p_dat,p21,all=T)
p_dat <- merge(p_dat,p22,all=T)
p_dat <- merge(p_dat,p23,all=T)
p_dat <- merge(p_dat,p24,all=T)

p_dat$macro <- "N"

m_dat$time <- NA
p_dat$time <- NA

m_dat[grep("T2", m_dat$Sample), ]$time <- 'T2'
m_dat[grep("T6", m_dat$Sample), ]$time <- 'T6'
p_dat[grep("T2", p_dat$Sample), ]$time <- 'T2'
p_dat[grep("T6", p_dat$Sample), ]$time <- 'T6'

mp_dat <- merge(m_dat, p_dat, all = T)

mp_dat$Pos <- as.factor(mp_dat$Pos)
summary(mp_dat$Pos)
#19508 and 47349 occur in every sample - filter out

mp_dat <- filter(mp_dat, ! Pos == '19508')
mp_dat <- filter(mp_dat, ! Pos == '47349')

(summary(as.factor(mp_dat$macro))) #No macrophages = 22 variants. Macrophages = 29 variants

mp_dat2 <- mp_dat
mp_dat2 <- filter(mp_dat2, time == 'T6')
mp_dat2 <- mp_dat2[,c(1, 3, 8, 9, 10,11,18,21,23)]

mp_dat2 <- mp_dat2[,c(1,9,2:8)]

mp_dat2_flex2 <- flextable(mp_dat2) %>%
  set_header_labels(macro = "Macrophages", Pos = "Position", Var_frac = "Frequency", Ref_nt = "Reference", Var_nt = "Variant", snpEff_type = "SNP type") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

mp_dat6 <- filter(mp_dat, time == 'T6')

mp_dat6$Var_frac <- mp_dat6$Var_frac*100

mut3 <- group_by(mp_dat6, macro, Sample) %>%
  summarise(., dist = sum(Var_frac))

m1 <- lm(dist ~ macro, data = mut3)
m2 <- lm(dist ~ 1, data = mut3)
anova(m1,m2)
plot(m1)
emmeans::emmeans(m1, pairwise ~ macro)

summary(as.factor(mp_dat6$macro)) #No macrophages, 15 variants. Macrophages 24 variants

p2 <- ggplot(mut3, aes(x = macro, y = dist)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black')+
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 13, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 13, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 13), legend.title = element_text(size = 13)) +
  ylab('Distance from ancestor (original phage)') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Macrophages absent', 'Macrophages present')) +
  labs(title = '(b)')

library(viridis)

mp_dat6$nam <- interaction(mp_dat6$macro, mp_dat6$Sample)

labs <- c('1', '2', '3\n          Macrophages absent', '4', '5', '6', '1', '2', '3\n          Macrophages present', '4', '5', '6')

p1 <- ggplot(mp_dat6, aes(x = interaction(Sample, macro), y = interaction(Pos, Function), fill = Var_frac)) +
  geom_tile(width = 0.75) +
  scale_fill_viridis(option = "A") +
  theme_bw() + 
  ylab("Gene") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 8, hjust = 0), axis.text.y = element_text(size = 8, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), 
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 8), axis.text.x = element_text(size = 12, colour = 'black'), legend.position = 'bottom') +
  labs(fill = "Frequency") +
  scale_x_discrete(labels = labs) + 
  xlab('Sample') +
  geom_vline(aes(xintercept = 6.5), size = 0.5) +
  labs(title = '(a)')

p1 + p2 + patchwork::plot_layout()
