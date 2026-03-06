this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(gtable)
library(lme4)
library(tidyverse)
library(lmerTest)
library(bootstrap)
library(ggpubr)
library(stringr)
library(brms)
library(BayesFactor)

`%notin%` <- Negate(`%in%`)
raw_data_path <- "thattracefull-trials.csv"
data<-read.csv(raw_data_path)
data_no_fills <- data %>%
  filter(condition %in% c("null", "that"))%>%
  mutate(block_number = (trial_sequence_total - 1) %/% 4 + 1)

data_summary <- data_no_fills %>%
  group_by(block_number, condition) %>%
  summarise(response = mean(response))

cbPalette = c("#e69d00", "#009e74","#d55e00",  "#cc79a7", "#0071b2")

satiation_plot <- ggplot(data_no_fills, aes(x = block_number, y=response, fill=condition)) +
  geom_point(data=data_summary,alpha=.9) +
  xlab("block number") +
  ylab("average acceptability")+
  geom_smooth(method=lm) +
  scale_fill_manual(values=cbPalette) +
  theme_bw()
satiation_plot

