# Load libraries ----------------------------------------------------------
library(arrow)
library(tidyverse)
library(patchwork) # Not needed?
library(ggh4x)

BUCKET <- Sys.getenv('WORKSPACE_BUCKET')

# Load OR tables ----------------------------------------------------------

assay_or_df <- read_parquet(
  paste0(BUCKET, "/or-estimates/or-estimates-2025-12-20.parquet")
) %>%
  as_tibble()

predictor_or_df <- read_parquet(
  paste0(BUCKET, "/or-estimates/or-estimates-2025-12-11_predictors.parquet")
)

# Broad gene-phenotype classes

gene_groups_df <- tribble(
  ~Gene, ~Disease,
  "BAP1", "Cancer",
  "BARD1", "Cancer",
  "BRCA1", "Cancer",
  "BRCA2", "Cancer",
  "CALM3", "Cardiovascular",
  "CHEK2", "Cancer",
  "G6PD", "Metabolic",
  "GCK", "Metabolic",
  "KCNE1", "Cardiovascular",
  "KCNH2", "Cardiovascular",
  "KCNQ4", "Hearing loss", #Rare disease",
  "MSH2", "Cancer",
  "OTC", "Metabolic",
  "PALB2", "Cancer",
  "PTEN", "Cancer",
  "RAD51C", "Cancer",
  "RAD51D", "Cancer",
  "SCN5A", "Cardiovascular",
  "TARDBP", "Rare disease",
  "TP53", "Cancer",
  "TSC2", "Cancer"
)

# SGE figure --------------------------------------------------------------

sge_theme <-
  theme_linedraw() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.title = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 20, face = 'bold', color = 'black'),
    strip.background = element_blank()
  )

## Condensed version ----

sge_plot_together_df <- assay_or_df %>%
  filter(
    str_ends(Dataset, 'unpublished'),
    Classifier == 'StandardizedClass',
    `Cases with variants` > 0,
    Gene != 'G6PD', # G6PD does not have enough func abnormal data
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    Classification = factor(
      str_c('Functionally ', str_to_title(Classification)),
      levels = c('Functionally Normal', 'Functionally Abnormal')
    ),
    Gene = factor(
      case_when(
        Dataset == "TSC2_rapgap_unpublished" ~ "TSC2 RapGAP",
        .default = Gene
      ),
      levels = c('TSC2 RapGAP', 'BARD1', 'PALB2', 'RAD51D')
    )
  )

# Calc limits for small plots
limits_df <- sge_plot_together_df %>% 
  filter(Gene != "TSC2 RapGAP") %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  )

make_sge_plot <- function(plot_df, x_limits)
  ggplot(
      plot_df,
      aes(
        x=`Odds Ratio`,
        xmin=OR_LI,
        xmax=OR_UI,
        y=Classification
      )
    ) +
    facet_grid(cols = vars(Gene), scales = 'free_x') +
    geom_pointrange() +
    geom_errorbar(width = 0.2) +
    scale_x_log10(
      #labels = scales::label_number(accuracy = 0.1),
      minor_breaks = NULL,
      #minor_breaks = scales::minor_breaks_log(),
      guide = "axis_logticks"
    ) +
    geom_blank(aes(x = x_limits)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
    geom_vline(xintercept = 1, linetype = 'dashed')

sge_together_plot <- make_sge_plot(sge_plot_together_df, deframe(limits_df))

print(sge_together_plot + sge_theme)

ggsave('Fig3_O1_plot.pdf', sge_together_plot + sge_theme, width = 12, height = 3)


## Broken out version ----

sge_plot_breakdown_df <- assay_or_df %>%
  filter(
    str_ends(Dataset, 'unpublished'),
    str_starts(Classifier, 'StandardizedClass'),
    `Cases with variants` > 0,
    Gene != 'G6PD', # G6PD does not have enough func abnormal data
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    Classification = factor(
      str_c(
        'Functionally',
        str_to_title(Classification),
        case_match(
          Classifier,
          "StandardizedClass" ~ "(all)",
          "StandardizedClass (stop gain only)" ~ "(truncating)",
          "StandardizedClass (without stop gain)" ~ "(non-truncating)"
        ),
        sep = ' '
      ),
      levels = str_c(
        'Functionally ', c(
          'Normal (all)',
          'Abnormal (all)',
          'Abnormal (non-truncating)',
          'Abnormal (truncating)'
        )
      )
    ),
    Gene = factor(
      case_when(
        Dataset == "TSC2_rapgap_unpublished" ~ "TSC2 RapGAP",
        .default = Gene
      ),
      levels = c('TSC2 RapGAP', 'BARD1', 'PALB2', 'RAD51D')
    )
  )

sge_breakdown_plot_v1 <- make_sge_plot(
  sge_plot_breakdown_df %>% filter(
    Classification %in% str_c('Functionally ', c(
      'Normal (all)',
      'Abnormal (all)',
      'Abnormal (non-truncating)'
    ))
  ),
  deframe(limits_df)
)

print(sge_breakdown_plot_v1 + sge_theme)

sge_breakdown_plot_v2 <- make_sge_plot(
  sge_plot_breakdown_df %>% filter(
    Classification %in% str_c('Functionally ', c(
      'Normal (all)',
      'Abnormal (all)',
      'Abnormal (truncating)'
    ))
  ),
  deframe(limits_df)
)

print(sge_breakdown_plot_v2 + sge_theme)


sge_breakdown_plot_v3 <- make_sge_plot(
  sge_plot_breakdown_df %>% filter(
    Classification %in% str_c('Functionally ', c(
      'Normal (all)',
      'Abnormal (non-truncating)',
      'Abnormal (truncating)'
    ))
  ),
  deframe(limits_df)
)

print(sge_breakdown_plot_v3 + sge_theme)

ggsave('Fig3_O2_plot.pdf', sge_breakdown_plot_v1 + sge_theme, width = 12, height = 4)
ggsave('Fig3_O3_plot.pdf', sge_breakdown_plot_v2 + sge_theme, width = 12, height = 4)
ggsave('Fig3_O4_plot.pdf', sge_breakdown_plot_v3 + sge_theme, width = 12, height = 4)


# Predictor calibrations --------------------------------------------------

predictor_theme <-
  theme_linedraw() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.title = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 20, face = 'bold', color = 'black'),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

msh2_table <- predictor_or_df %>%
  filter(
    Dataset == 'REVEL',
    Gene == 'MSH2'
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI)
  )

# Calibration comparison

vep_levels = c(
  "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4"
)

predictor_plot_df <- predictor_or_df %>%
  filter(
    `Few samples` == FALSE,
    #`Cases with variants` > 0
    Classification %in% vep_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = vep_levels
    )
  ) %>%
  left_join(gene_groups_df)

make_predictor_plot <- function(in_df, scales = 'free_x')
  ggplot(
    in_df,
    aes(
      x=`Odds Ratio`,
      xmin=OR_LI,
      xmax=OR_UI,
      y=Classification,
      color = Classifier,
      shape = significance
    )
  ) +
  facet_grid(
    cols = vars(Gene),
    rows = vars(Dataset),
    scales = scales
  ) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks"
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

predictor_plot <- make_predictor_plot(predictor_plot_df)

print(predictor_plot + predictor_theme)

ggsave(
  'Predictor_plot.pdf',
  predictor_plot + predictor_theme,
  width = 40,
  height = 6,
  device = cairo_pdf)

## Condensed

condensed_predictor_plot <- predictor_plot_df %>%
  filter(!(Gene %in% c('G6PD', 'KCNH2'))) %>% # Poorly defined phenotypes
  inner_join(
    predictor_plot_df %>%
      group_by(Gene) %>%
      summarize(keep = any(significance == 'Significant at 95%')) %>%
      filter(keep) %>%
      select(Gene)
  ) %>%
  left_join(gene_groups_df) %>%
  make_predictor_plot() + geom_blank(aes(x=deframe(
    predictor_plot_df %>%
      filter(Gene == 'MSH2') %>%
      summarize(max(OR_UI))
  )))

print(condensed_predictor_plot + predictor_theme)

ggsave(
  'Condensed_predictor_plot.pdf',
  condensed_predictor_plot + predictor_theme,
  width = 15,
  height = 6,
  device = cairo_pdf)

# Assay OddsPath vs Zeiberg -----------------------------------------------

assay_classification_levels = c(
  "NORMAL", "ABNORMAL",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8"
)

assay_plot_df <-
  assay_or_df %>%
  filter(
    `Few samples` == FALSE,
    #Classifier %in% c('Calibrated (2025-12-08)', 'StandardizedClass'),
    #`Cases with variants` > 0
    Classification %in% assay_classification_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = assay_classification_levels
    )
  )

assay_plot <- ggplot(
  assay_plot_df,
  aes(
    x=`Odds Ratio`,
    xmin=OR_LI,
    xmax=OR_UI,
    y=Classification,
    shape = significance
  )
) +
  facet_nested(cols = vars(Gene, Dataset), rows = vars(Classifier), scales = 'free') +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    #labels = scales::label_number(accuracy = 0.1),
    guide = "axis_logticks"
  ) +
  scale_shape_manual(
    #limits = c('Not significant at 95%', 'Significant at 95%'),
    #minor_breaks = c('circle open'),
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

print(assay_plot)

# Condensed plot
condensed_assay_plot_df <-
  assay_plot_df %>%
  filter(
    Classifier %in% c(
      'Calibrated (2025-12-08)',
      'OP_points_18_25'
    ),
    Dataset %in% c(
      'BAP1_Waters_2024',
      'BARD1_unpublished',
      #'BRCA1_Adamovich_2022_Cisplatin',
      #'BRCA1_Adamovich_2022_HDR',
      'BRCA1_Findlay_2018',
      'BRCA2_Hu_2024',
      #'BRCA2_Sahu_2025_HDR',
      'KCNQ4_Zheng_2022_current_homozygous',
      'MSH2_Jia_2021',
      'PALB2_unpublished',
      'RAD51C_Olvera-León_2024_z_score_D4_D14',
      'RAD51D_unpublished',
      'TP53_Fayer_2021_meta'
      #'TP53_Boettcher_2019'
      #'TP53_Fortuno_2021_Kato_meta',
      #'TP53_Giacomelli_2018_combined_score'
      #'TSC2_rapgap_unpublished'
    ),
    Classification != "0"
  ) %>%
  left_join(gene_groups_df)

# Cancer only for v2
condensed_assay_plot2_df <-
  condensed_assay_plot_df %>%
  filter(Disease == 'Cancer')

# Limits for most panels
assay_plot_common_limits <- condensed_assay_plot_df %>% 
  filter(!(Gene %in% c("BRCA1", "MSH2"))) %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  ) %>% deframe()

condensed_assay_plot <- ggplot(
  condensed_assay_plot_df,
  aes(
    x=`Odds Ratio`,
    xmin=OR_LI,
    xmax=OR_UI,
    y=Classification,
    shape = significance
  )
) +
  facet_nested(
    cols = vars(Disease, Gene),
    scales = 'free',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    )
  ) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    #limits = c('Not significant at 95%', 'Significant at 95%'),
    #minor_breaks = c('circle open'),
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_blank(aes(x = assay_plot_common_limits)) +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

assay_theme = predictor_theme

print(condensed_assay_plot + assay_theme)

ggsave(
  'Condensed_assay_plot_v1.pdf',
  condensed_assay_plot +
    assay_theme +
    theme(ggh4x.facet.nestline = element_line(color="black")),
  width = 15,
  height = 6.5,
  device = cairo_pdf)

# Cancer only

condensed_assay_plot2 <- ggplot(
  condensed_assay_plot2_df,
  aes(
    x=`Odds Ratio`,
    xmin=OR_LI,
    xmax=OR_UI,
    y=Classification,
    shape = significance
  )
) +
  facet_grid(
    cols = vars(Gene),
    scales = 'free',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    )
  ) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    #limits = c('Not significant at 95%', 'Significant at 95%'),
    #minor_breaks = c('circle open'),
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_blank(aes(x = assay_plot_common_limits)) +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

print(condensed_assay_plot2 + assay_theme)

ggsave(
  'Cancer_assay_plot_v1.pdf',
  condensed_assay_plot2 + assay_theme,
  width = 15,
  height = 6,
  device = cairo_pdf)

# Combined points ---------------------------------------------------------

points_levels = c(
  "≤ -16", "≤ -15", "≤ -14", "≤ -13", "≤ -12", "≤ -11", "≤ -10", "≤ -9",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8",
  "≥ +9", "≥ +10", "≥ +11", "≥ +12", "≥ +13", "≥ +14", "≥ +15", "≥ +16"
)

combined_points_plot_df <-
  assay_or_df %>%
  filter(
    Dataset == 'Combined points',
    `Few samples` == FALSE #,
    #Classifier %in% c('Calibrated (2025-12-08)', 'StandardizedClass'),
    #`Cases with variants` > 0
    #Classification %in% assay_classification_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = points_levels
    )
  )

combined_points_plot <- ggplot(
  combined_points_plot_df,
  aes(
    x = `Odds Ratio`,
    xmin = OR_LI,
    xmax = OR_UI,
    y = Classification,
    shape = significance
  )
) +
  facet_grid(
    cols = vars(Gene),
    scales = 'free'
  ) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    #limits = c('Not significant at 95%', 'Significant at 95%'),
    #minor_breaks = c('circle open'),
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  #geom_blank(aes(x = assay_plot_common_limits)) +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

print(combined_points_plot)

combined_points_condensed_plot_df <-
  combined_points_plot_df %>%
  filter(
    !(Gene %in% c(
      "G6PD", "KCNH2", # Poor phenotype definition
      "GCK", "KCNE1", "PTEN", "SCN5A", "TSC2" # No sig. intervals
    )),
    Classification != "0"
  ) %>%
  left_join(gene_groups_df)

combined_points_condensed_plot <- ggplot(
  combined_points_condensed_plot_df,
  aes(
    x = `Odds Ratio`,
    xmin = OR_LI,
    xmax = OR_UI,
    y = Classification,
    shape = significance
  )
) +
  facet_nested(
    cols = vars(Disease, Gene),
    scales = 'free',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    )
  ) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    #limits = c('Not significant at 95%', 'Significant at 95%'),
    #minor_breaks = c('circle open'),
    values = c(
      'Not significant at 95%' = 'circle open',
      'Significant at 95%' = 'circle'),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  #geom_blank(aes(x = 0.25)) +
  geom_blank(aes(x = 6)) +
  guides(color = guide_legend(
    position = 'bottom',
    override.aes = aes(shape = 'circle open')
  ))

print(
  combined_points_condensed_plot +
    assay_theme +
    theme(ggh4x.facet.nestline = element_line(color="black"))
)

ggsave(
  'Combined_points_plot.pdf',
  combined_points_condensed_plot +
    assay_theme +
    theme(ggh4x.facet.nestline = element_line(color="black")),
  width = 15,
  height = 8,
  device = cairo_pdf)

# Fancy phenotype grouping

pheno_df <-
  combined_points_condensed_plot_df %>%
  select(Gene, `Case inclusion phenotypes`) %>%
  unique() %>%
  unnest(`Case inclusion phenotypes`) %>%
  cbind(Yes='X') %>%
  pivot_wider(
    names_from = `Case inclusion phenotypes`,
    values_from = Yes, values_fill = '.')

View(pheno_df)