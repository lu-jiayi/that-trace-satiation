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
raw_data_path <- "thattracenewfull-trials.csv"
data<-read.csv(raw_data_path)
data <- data %>%
  filter(!is.na(item)) %>%
  mutate(
    comp = case_when(
      sentence_id %/% 1000 == 1 & sentence_id %% 10 %in% c(1, 3) ~ "that",
      sentence_id %/% 1000 == 1 & sentence_id %% 10 %in% c(2, 4) ~ "null",
      TRUE ~ NA_character_
    ),
    
    gap = case_when(
      sentence_id %/% 1000 == 1 & sentence_id %% 10 %in% c(1, 2) ~ "subject",
      sentence_id %/% 1000 == 1 & sentence_id %% 10 %in% c(3, 4) ~ "object",
      TRUE ~ NA_character_
    )
  )
prepost <- data %>%
  filter(trial_sequence_total < 5 | trial_sequence_total > 28 ) %>%
  mutate(block = case_when(
    trial_sequence_total < 5 ~ "pre-exposure", 
    trial_sequence_total > 28 ~ "post-exposure"))
exposure <- data %>%
  filter(between(trial_sequence_total, 5, 28)) %>%
  arrange(workerid, trial_sequence_total) %>%
  group_by(workerid) %>%
  mutate(
    trial_order = cumsum(!is.na(gap))
  ) %>%
  filter(!is.na(gap)) %>%
  mutate(order = ((trial_order-1) %% 8 ) + 1)%>%
  ungroup()
exposure_summary <- exposure %>%
  group_by(order) %>%
  dplyr::summarise(response = mean(response))

exposure_plot <- ggplot(
  exposure,
  aes(x = order, y = response)
) +
  geom_point(data = exposure_summary, alpha = .9) +
  xlab("Presentation order") +
  ylab("Average acceptability") +
  geom_smooth(method = lm) +
  theme_bw() 

exposure_plot



prepost_summarize <- prepost %>%
  group_by(block, comp, gap) %>%
  dplyr::summarise(
    Mean = mean(response),
    CILow = ci.low(response),
    CIHigh = ci.high(response),
    .groups = "drop"
  ) %>%
  mutate(YMin = Mean - CILow, YMax = Mean + CIHigh)%>%
  mutate(
    block = factor(block, levels = c("pre-exposure", "post-exposure"))
  )
cbPalette = c("#e69d00", "#009e74","#d55e00",  "#cc79a7", "#0071b2")
bar_plot <- ggplot(
  prepost_summarize,
  aes(x = comp, y = Mean, fill = gap)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = YMin, ymax = YMax),
    position = position_dodge(width = 0.8),
    width = 0.3,
    linewidth = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~ block) +
  scale_fill_manual(values = cbPalette, name = "Gap position") +
  theme_bw() +
  xlab("Complementizer type") +
  ylab("Mean acceptability") +
  ylim(0, 1) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )


bar_plot

prepost$gap <- factor(prepost$gap, levels = c("object", "subject"))
prepost$comp <- factor(prepost$comp, levels = c("null", "that"))
prepost$block <- factor(prepost$block, levels = c("pre-exposure", "post-exposure"))
contrasts(prepost$gap) <- contr.sum(2)
contrasts(prepost$comp) <- contr.sum(2)
contrasts(prepost$block) <- contr.sum(2)
model <- lmer(response ~ gap*comp*block + 
                (1 + gap*comp*block - block - gap:comp:block|item) + 
                (1 + gap*comp*block|workerid), data = prepost)
summary(model)
