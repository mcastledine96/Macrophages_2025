library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(vcfR)
library(vegan)

p1 <- read_tsv("bact_vars/b1_vs.tsv")
p1$Sample <- 'B1'
p2 <- read_tsv("bact_vars/b2_vs.tsv")
p2$Sample <- 'B2'
p3 <- read_tsv("bact_vars/b3_vs.tsv")
p3$Sample <- 'B3'
p4 <- read_tsv("bact_vars/b4_vs.tsv")
p4$Sample <- 'B4'
p5 <- read_tsv("bact_vars/bm1_vs.tsv")
p5$Sample <- 'BM1'
p6 <- read_tsv("bact_vars/bm2_vs.tsv")
p6$Sample <- 'BM2'
p7 <- read_tsv("bact_vars/bm3_vs.tsv")
p7$Sample <- 'BM3'
p8 <- read_tsv("bact_vars/bp1_vs.tsv")
p8$Sample <- 'BP1'
p9 <- read_tsv("bact_vars/bp2_vs.tsv")
p9$Sample <- 'BP2'
p10 <- read_tsv("bact_vars/bp3_vs.tsv")
p10$Sample <- 'BP3'
p11 <- read_tsv("bact_vars/bp4_vs.tsv")
p11$Sample <- 'BP4'
p12 <- read_tsv("bact_vars/bp5_vs.tsv")
p12$Sample <- 'BP5'
p13 <- read_tsv("bact_vars/bp6_vs.tsv")
p13$Sample <- 'BP6'
p14 <- read_tsv("bact_vars/bmp1_vs.tsv")
p14$Sample <- 'BMP1'
p15 <- read_tsv("bact_vars/bmp2_vs.tsv")
p15$Sample <- 'BMP2'
p16 <- read_tsv("bact_vars/bmp3_vs.tsv")
p16$Sample <- 'BMP3'
p17 <- read_tsv("bact_vars/bmp4_vs.tsv")
p17$Sample <- 'BMP4'
p18 <- read_tsv("bact_vars/bmp5_vs.tsv")
p18$Sample <- 'BMP5'
p19 <- read_tsv("bact_vars/bmp6_vs.tsv")
p19$Sample <- 'BMP6'

bdat <- merge(p1,p2,all=T)
bdat <- merge(bdat, p3, all = T)
bdat <- merge(bdat, p4, all = T)
bdat$phage <- 'N'
bdat$macro <- 'N'

bmdat <- merge(p5,p6,all=T)
bmdat <- merge(bmdat,p7,all=T)
bmdat$phage <- 'N'
bmdat$macro <- 'Y'

bpdat <- merge(p8,p9,all=T)
bpdat2 <- merge(p10,p11,all=T)
bpdat3 <- merge(p12,p13,all=T)
bpdat <- merge(bpdat,bpdat2, all=T)
bpdat <- merge(bpdat,bpdat3, all=T)
bpdat$phage <- 'Y'
bpdat$macro <- 'N'

bmpdat <- merge(p14,p15,all=T)
bmpdat2 <- merge(p16,p17,all=T)
bmpdat3 <- merge(p18,p19,all=T)
bmpdat <- merge(bmpdat,bmpdat2, all=T)
bmpdat <- merge(bmpdat,bmpdat3, all=T)
bmpdat$phage <- 'Y'
bmpdat$macro <- 'Y'

bseq <- merge(bdat, bmdat, all = T)
bseq2 <- merge(bpdat, bmpdat, all = T)
bseq <- merge(bseq, bseq2, all = T)

summary(as.factor(bseq$Pos))

bseq_f <- bseq %>%
  group_by(Pos) %>% # Group by the column to check repeats
  filter(n() != 19) %>% # Keep groups that do not have exactly 19 occurrences
  ungroup()

bseq_f <- filter(bseq_f, ! Function == 'Rhs-family protein')

bseq_f$variant <- interaction(bseq_f$Pos, bseq_f$Function)

bseq_f <- bseq_f %>%
  mutate(y_axis_label = paste0(variant, "\n", Upstream_feature, "\n", Downstream_feature))

bseq_f <- bseq_f %>%
  mutate(y_axis_label2 = paste0(Function, "\n", Upstream_feature, "\n", Downstream_feature))

ggplot(bseq_f, aes(x = Sample, y = y_axis_label)) +
  geom_tile(aes(fill = Var_frac)) 

bseq_f$Sample2 <- interaction(bseq_f$Sample, bseq_f$macro)

#table of mutations
bseq_ft <- bseq_f

bseq_ft <- bseq_ft[,c(3,6,8,9,10,11,18:21,23,25)]

bseq_ft <- bseq_ft[,c(11,12,1:10)]

bseq_ft_flex <- flextable(bseq_ft) %>%
  set_header_labels(macro = "Macrophages", Pos = "Position", Var_frac = "Frequency", Ref_nt = "Reference", Var_nt = "Variant", snpEff_type = "SNP type") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 11, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#

new_values <- c("BP2.N", "BP5.N", "BP6.N", 'BMP1.Y', 'BMP3.Y', 'BMP6.Y')

library(tibble)
bseq_f2 <- bseq_f %>%
  add_row(Sample2 = new_values)

nams <- unique(bseq_f2$y_axis_label2)[1:4]

labs <- c('1', '2', '3\n          Macrophages present', '4', '5', '6', '1', '2', '3\n          Macrophages absent', '4', '5', '6')
library(viridis)
ggplot(bseq_f2, aes(x = Sample2, y = y_axis_label2)) +
  geom_tile(aes(fill = Var_frac)) +
  scale_y_discrete(limits = nams) +
  scale_fill_viridis(option = "A") +
  theme_bw()  +
  theme(axis.text = element_text(size = 10, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), 
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
  geom_vline(aes(xintercept = 6.5), size = 0.5) +
  scale_x_discrete(labels = labs) +
  ylab('Mutated gene\nUpstream gene\nDownstream gene') +
  xlab('Treatment replicate') +
  labs(fill = 'Frequency')

nams <- unique(bseq_f2$Function)[1:4]

ggplot(bseq_f2, aes(x = Sample2, y = Function)) +
  geom_tile(aes(fill = Var_frac)) +
  scale_y_discrete(limits = nams) +
  scale_fill_viridis(option = "A") +
  theme_bw()  +
  theme(axis.text = element_text(size = 13, colour = 'black'), axis.title = element_text(size = 15, colour = 'black'), 
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
  geom_vline(aes(xintercept = 6.5), size = 0.5) +
  scale_x_discrete(labels = labs) +
  ylab('Gene')
