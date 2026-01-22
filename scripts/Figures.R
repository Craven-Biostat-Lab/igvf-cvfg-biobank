# Install libraries
update.packages("ggplot2")
install.packages("patchwork")
install.packages("ggh4x")
install.packages("writexl")

# Load libraries ----------------------------------------------------------
library(arrow)
library(tidyverse)
library(patchwork)
library(ggh4x)

BUCKET <- Sys.getenv('WORKSPACE_BUCKET')
SAVE <- FALSE

# Load OR tables ----------------------------------------------------------

or_df <- read_parquet(
  paste0(BUCKET, "/or-estimates/publication-or-estimates-2026-01-20.parquet")
)

# Broad gene-phenotype classes

gene_groups_df <- tribble(
  ~Gene, ~Disease, ~D2,
  "BAP1", "Cancer", "Cancer",
  "BARD1", "Cancer", "Cancer",
  "BRCA1", "Cancer", "Cancer",
  "BRCA2", "Cancer", "Cancer",
  "CALM3", "Cardiovascular", "Cardio- vascular",
  "CHEK2", "Cancer", "Cancer",
  "G6PD", "Metabolic", "Meta- bolic",
  "GCK", "Metabolic", "Meta- bolic",
  "KCNE1", "Cardiovascular", "Cardio- vascular",
  "KCNH2", "Cardiovascular", "Cardio- vascular",
  "KCNQ4", "Hearing loss", "Hearing loss", #Rare disease",
  "MSH2", "Cancer", "Cancer",
  "OTC", "Metabolic", "Meta- bolic",
  "PALB2", "Cancer", "Cancer",
  "PTEN", "Cancer", "Cancer",
  "RAD51C", "Cancer", "Cancer",
  "RAD51D", "Cancer", "Cancer",
  "SCN5A", "Cardiovascular", "Cardio- vascular",
  "TARDBP", "Rare disease", "Rare disease",
  "TP53", "Cancer", "Cancer",
  "TSC2", "Cancer", "Cancer"
)

# Colors for predictor calibrations

color_mapping <- c(
  `single gene` = "#d62728", 
  `domain aggregate` = "#f5971e",
  `gene aggregate` = "#e9dfc6"
)

# Figure theme
point_in_mm = 0.3527778
nature_theme <- theme_linedraw() +
  theme(
    text = element_text(family = 'Helvetica', size = 7),
    line = element_line(linewidth = 1*point_in_mm),
    geom = element_geom(
      linewidth = 1*point_in_mm,
      borderwidth = 1*point_in_mm,
      pointsize = 3*point_in_mm
    ),
    #spacing = unit(3, 'points'),
    axis.ticks = element_line(linewidth = 0.5*point_in_mm),
    #axis.title = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = 'bold', color = 'black', margin = margin(3,3,3,3)),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    ggh4x.facet.nestline = element_line(linewidth = 1*point_in_mm, color="black"),
    plot.tag = element_text(face = 'bold')
  )

# SGE figure --------------------------------------------------------------

## Condensed version ----

sge_plot_together_df <- or_df %>%
  filter(
    str_ends(Dataset, 'unpublished'),
    Classifier == 'StandardizedClass',
    Powered,
    #`Cases with variants` > 0,
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
      levels = c('TSC2 RapGAP', 'BARD1', 'PALB2', 'RAD51D', 'XRCC2', 'CTCF', 'SFPQ')
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

print(sge_together_plot + nature_theme)

if (SAVE)
  ggsave('Fig3_O1_plot.pdf', sge_together_plot + sge_theme, width = 12, height = 3)


## Broken out version ----

sge_plot_breakdown_df <- or_df %>%
  filter(
    str_ends(Dataset, 'unpublished'),
    str_starts(Classifier, 'StandardizedClass'),
    Powered, #`Cases with variants` > 0,
    !(Gene %in% c('G6PD', 'CTCF')) # These genes do not have enough func abnormal data
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

fig2j_df <- sge_plot_breakdown_df %>% filter(
  Classification %in% str_c('Functionally ', c(
    'Normal (all)',
    'Abnormal (non-truncating)',
    'Abnormal (truncating)'
  ))
)

fig2j_plot <- make_sge_plot(fig2j_df, deframe(limits_df))

print(fig2j_plot + nature_theme)

if (SAVE) {
  ggsave('fig2j.pdf', fig2j_plot + nature_theme, width = 100, height = 40, units = 'mm')
}

# Predictor calibrations --------------------------------------------------

msh2_table <- or_df %>%
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

predictor_plot_df <- or_df %>%
  filter(
    Powered,
    #`Few samples` == FALSE,
    #`Cases with variants` > 0,
    Classifier %in% c('Gene-aggregated calibration', 'Gene-specific calibration'),
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

print(predictor_plot)

if (SAVE)
  ggsave(
    'Predictor_plot.pdf',
    predictor_plot + predictor_theme,
    width = 40,
    height = 6,
    device = cairo_pdf
  )

figure_predictor_plot_df <- predictor_plot_df %>%
  filter(
    Gene %in% c('BRCA1', 'BRCA2', 'MSH2', 'TP53'), # Have gene-specific calibration
  )

## Final figure
figure_predictor_plot <- figure_predictor_plot_df %>%
  make_predictor_plot() + geom_blank(aes(x=deframe(
    predictor_plot_df %>%
      filter(Gene == 'MSH2') %>%
      summarize(max(OR_UI))
  )))

print(figure_predictor_plot + nature_theme)

# Assay OddsPath vs Zeiberg -----------------------------------------------

assay_classification_levels = c(
  "NORMAL", "ABNORMAL",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8"
)

assay_plot_df <-
  or_df %>%
  filter(
    Powered, #`Few samples` == FALSE,
    Classifier %in% c('ExCALIBR', 'StandardizedClass'),
    `Cases with variants` > 0,
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


# Calibrated w/ significant intervals only

# Condensed plot
condensed_assay_plot_df <-
  assay_plot_df %>%
  filter(
    Classifier == 'ExCALIBR',
    Dataset %in% c(
      'BAP1_Waters_2024',
      'BARD1_unpublished',
      #'BRCA1_Adamovich_2022_Cisplatin',
      #'BRCA1_Adamovich_2022_HDR',
      'BRCA1_Findlay_2018',
      'BRCA2_Hu_2024',
      #'BRCA2_Sahu_2025_HDR',
      'GCK_Gersing_2023_complementation',
      'KCNH2_Jiang_2022',
      'KCNQ4_Zheng_2022_current_homozygous',
      'MSH2_Jia_2021',
      'PALB2_unpublished',
      'RAD51C_Olvera-León_2024_z_score_D4_D14',
      'RAD51D_unpublished',
      'SCN5A',
      'TP53_Fayer_2021_meta',
      #'TP53_Boettcher_2019'
      #'TP53_Fortuno_2021_Kato_meta',
      #'TP53_Giacomelli_2018_combined_score'
      'TSC2_rapgap_unpublished'
    ),
    Classification != "0"
  ) %>%
  left_join(gene_groups_df)

# Limits for most panels
assay_plot_common_limits <- condensed_assay_plot_df %>% 
  filter(!(Gene %in% c("BRCA1", "MSH2", "KCNH2", "TSC2", "GCK"))) %>%
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
  facet_nested_wrap(
    facets = vars(Disease, Gene),
    nrow = 3,
    scales = 'free_x',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    ),
    strip = strip_nested(text_x = element_text(margin = margin(0,0,0,0)))
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
    guide = 'none'
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_blank(aes(x = assay_plot_common_limits))

print(condensed_assay_plot)

# Supplemental figure together --------------------------------------------

fig_exd3 <- (figure_predictor_plot + nature_theme) /
  (condensed_assay_plot + nature_theme) +
  plot_annotation(tag_levels='a') +
  plot_layout(heights = c(4,9))

print(fig_exd3)

if(SAVE)
  ggsave(
    'Extended Data Figure 3.pdf',
    fig_exd3,
    width = 160, # Max 183
    height = 247,
    units = 'mm',
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
  or_df %>%
  filter(
    Dataset == 'Combined points',
    Powered,
    #`Few samples` == FALSE #,
    #Classifier %in% c('Calibrated (2025-12-08)', 'StandardizedClass'),
    `Cases with variants` > 0
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
      #"GCK", "KCNE1", "PTEN", "SCN5A", "TSC2" # No sig. intervals
      "SGCB" # Erroneously kept in spite of no sig. intervals
    )),
    Classification != "0"
  ) %>%
  inner_join(
    combined_points_plot_df %>%
      group_by(Gene) %>%
      summarize(keep = any(significance == 'Significant at 95%')) %>%
      filter(keep) %>%
      select(Gene)
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
    cols = vars(D2, Gene),
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
    nature_theme +
    theme(axis.text.x = element_text(angle = 30))
)
#    assay_theme +
#    theme(ggh4x.facet.nestline = element_line(color="black"))
#)

if(SAVE)
  ggsave(
    'Fig 6B 2026-01-16.pdf',
    combined_points_condensed_plot + nature_theme + theme(axis.text.x = element_text(angle = 45)),
    width = 183, # Max 183
    height = 90,
    units = 'mm',
    device = cairo_pdf)

ggsave(
  'Combined_points_plot.pdf',
  combined_points_condensed_plot +
    assay_theme +
    theme(ggh4x.facet.nestline = element_line(color="black")),
  width = 18,
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

# Build sup. tables --------
library(writexl)

sheets <- list(
  "Figure 2j" = fig2j_df %>%
    select(
      Gene,
      Classification,
      `Odds Ratio`,
      OR_LI,
      OR_UI,
      LogOR,
      LogOR_LI,
      LogOR_UI,
      `Wald p-value` = `p-value`
    ),
  "Figure 6b" = combined_points_condensed_plot_df %>%
    select(
      Gene,
      Classification,
      `Odds Ratio`,
      OR_LI,
      OR_UI,
      LogOR,
      LogOR_LI,
      LogOR_UI,
      `Wald p-value` = `p-value`
    ),
  "Extended Data Figure 3a" = condensed_assay_plot_df %>%
    select(
      Gene,
      Dataset,
      Classification,
      `Odds Ratio`,
      OR_LI,
      OR_UI,
      LogOR,
      LogOR_LI,
      LogOR_UI,
      `Wald p-value` = `p-value`
    ),
  "Extended Data Figure 3b" = figure_predictor_plot_df %>%
    select(
      Gene,
      Predictor = Dataset,
      Calibration = Classifier,
      Classification,
      `Odds Ratio`,
      OR_LI,
      OR_UI,
      LogOR,
      LogOR_LI,
      LogOR_UI,
      `Wald p-value` = `p-value`
    )
)

if(SAVE)
  write_xlsx(
    sheets,
    "Supplementary Table ORs.xlsx"
  )

# Notes
#fig2j_df # Fig 2j
#combined_points_condensed_plot_df # Fig 6b
#condensed_assay_plot_df # ExD Fig 3a
#figure_predictor_plot_df # Exd Fig 3b
