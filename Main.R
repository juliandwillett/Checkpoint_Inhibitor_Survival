library(survival)
library(survminer)
library(readxl)
library(writexl)
library(tidyverse)
library(viridis)
library(nph)
library(coxphf) # Added for Firth Penalized Likelihood
source('Functions.R')

# --- 1. Data Prep & Cleaning ---
surv_data <- read_excel("CaseRecurrenceSurvival.xlsx") %>%
  filter(Progression != "Unknown", !is.na(Progression), Case != 25) %>%
  mutate(
    `Date First Rx` = map_dbl(`Date First Rx`, to_excel_numeric),
    `Date Bx`       = map_dbl(`Date Bx`, to_excel_numeric),
    `Last FU`       = map_dbl(`Last FU`, to_excel_numeric),
    
    Prog_Num = ifelse(Progression == "None", NA, Progression),
    Prog_Num = map_dbl(Prog_Num, to_excel_numeric),
    
    Death_Num = ifelse(`Death from disease` == "Alive", NA, `Death from disease`),
    Death_Num = map_dbl(Death_Num, to_excel_numeric)
  ) %>%
  mutate(
    Progressed = ifelse(Progression == "None", 0, 1),
    Died = ifelse(`Death from disease` == "Alive", 0, 1),
    
    TTP = ifelse(Progressed == 1, 
                 Prog_Num - `Date First Rx`,
                 `Last FU` - `Date First Rx`),
    
    TTD = ifelse(Died == 1, 
                 Death_Num - `Date First Rx`,
                 `Last FU` - `Date First Rx`),
    
    Bx_to_Prog = ifelse(Progressed == 1 & !is.na(`Date Bx`), 
                        Prog_Num - `Date Bx`, 
                        NA),
    Bx_to_OS = ifelse(Died == 1 & !is.na(`Date Bx`), 
                      Death_Num - `Date Bx`, 
                      NA)
  ) %>%
  mutate(T4HS_Cat = ifelse(str_detect(Category,'T4HS'),'T4HS','NotT4HS'))

# Global Save Settings
IMG_WIDTH <- 8
IMG_HEIGHT <- 6
IMG_DPI <- 300

# ==============================================================================
# ANALYSIS 0.1: PFS (OVERALL COHORT)
# ==============================================================================
# Calculating overall progression-free survival for the entire dataset
fit_pfs_all <- survfit(Surv(TTP, Progressed) ~ 1, data = surv_data)

pfs_plot_all <- ggsurvplot(
  fit_pfs_all, 
  data = surv_data,
  pval = FALSE,       
  risk.table = TRUE,
  palette = c("#440154FF"), # Single color from viridis
  title = "Overall PFS: Entire Cohort",
  xlab = "Days from First Rx",
  legend = "none",    # Hide legend for a single curve
  ggtheme = theme_minimal(),
  conf.int = FALSE
)
print(pfs_plot_all)

# Save Overall PFS Plot
png("PFS_Entire_Cohort.png", width = IMG_WIDTH, height = IMG_HEIGHT, units = "in", res = IMG_DPI)
print(pfs_plot_all)
dev.off()

# ==============================================================================
# ANALYSIS 0.2: OS (OVERALL COHORT)
# ==============================================================================
# Calculating overall survival for the entire dataset
fit_os_all <- survfit(Surv(TTD, Died) ~ 1, data = surv_data)

os_plot_all <- ggsurvplot(
  fit_os_all, 
  data = surv_data,
  pval = FALSE,
  risk.table = TRUE,
  palette = c("#21908CFF"), # Distinct single color from viridis
  title = "Overall OS: Entire Cohort",
  xlab = "Days from First Rx",
  legend = "none",
  ggtheme = theme_minimal(),
  conf.int = FALSE
)
print(os_plot_all)

# Save Overall OS Plot
png("OS_Entire_Cohort.png", width = IMG_WIDTH, height = IMG_HEIGHT, units = "in", res = IMG_DPI)
print(os_plot_all)
dev.off()

# ==============================================================================
# ANALYSIS 1: PFS by T4HS
# ==============================================================================
fit_pfs <- survfit(Surv(TTP, Progressed) ~ T4HS_Cat, data = surv_data)
print(fit_pfs)

pfs_plot <- ggsurvplot(
  fit_pfs, 
  data = surv_data,
  pval = TRUE, 
  risk.table = TRUE,
  palette = "viridis",
  title = "PFS: Immunotherapy Start to Progression",
  xlab = "Days from First Rx",
  legend.labs = sort(unique(surv_data$T4HS_Cat)),
  ggtheme = theme_minimal()
)
print(pfs_plot)

# Save PFS Plot (using the survminer ggsave wrapper to catch the risk table)
png("PFS_T4HS_Comparison.png", width = IMG_WIDTH, height = IMG_HEIGHT, units = "in", res = IMG_DPI)
print(pfs_plot)
dev.off()

# Validate using Cox
cox_model_pfs <- coxph(Surv(TTP, Progressed) ~ T4HS_Cat, data = surv_data)
ph_test_pfs <- cox.zph(cox_model_pfs)
print("--- PROPORTIONAL HAZARDS TEST: PFS ---")
print(ph_test_pfs)

# Give late events more weight in regression, with T4HS mediated events anticipated to 
#.   occur later while T cell exhaustion is reversed.
fh_late_test <- logrank.test(
  time = surv_data$TTP,
  event = surv_data$Progressed,
  group = surv_data$T4HS_Cat,
  rho = 0,   # 0 means no extra weight on early events
  gamma = 1  # 1 means extra weight on late events
)
print(fh_late_test)

# Get hazards ratio
cox_model <- coxph(Surv(TTP, Progressed) ~ T4HS_Cat, data = surv_data)
summary(cox_model)

# ==============================================================================
# ANALYSIS 2: OS by T4HS cat
# ==============================================================================
fit_os <- survfit(Surv(TTD, Died) ~ T4HS_Cat, data = surv_data)

os_plot <- ggsurvplot(
  fit_os, 
  data = surv_data,
  pval = TRUE, 
  risk.table = TRUE,
  palette = "viridis",
  title = "OS: Immunotherapy Start to Death",
  xlab = "Days from First Rx",
  legend.labs = sort(unique(surv_data$T4HS_Cat)),
  ggtheme = theme_minimal()
)

print(os_plot)
png("OS_T4HS_Comparison.png", width = IMG_WIDTH, height = IMG_HEIGHT, units = "in", res = IMG_DPI)
print(os_plot)
dev.off()

# ==============================================================================
# ANALYSIS 3: BIOPSY DATE VS PROGRESSION DATE
# ==============================================================================
prog_only <- surv_data %>% filter(!is.na(Bx_to_Prog))

bx_prog_stats <- prog_only %>%
  group_by(T4HS_Cat) %>%
  summarise(
    N = n(),
    Median_Bx_to_Prog = median(Bx_to_Prog),
    IQR = IQR(Bx_to_Prog)
  )
print("--- BIOPSY TO PROGRESSION METRICS ---")
print(bx_prog_stats)

bx_prog_box <- ggplot(prog_only, aes(x = T4HS_Cat, y = Bx_to_Prog, fill = T4HS_Cat)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_viridis_d(option = "D") + 
  theme_minimal() +
  labs(
    title = "Interval: Biopsy Date to Progression Date",
    subtitle = "Positive = Bx before Prog; Negative = Bx after Prog",
    y = "Days Delta",
    x = "Morphologic Category"
  ) +
  theme(legend.position = "none")

print(bx_prog_box)

# Save Boxplot
ggsave("Biopsy_to_Progression_Boxplot.png", plot = bx_prog_box, width = IMG_WIDTH, height = IMG_HEIGHT, dpi = IMG_DPI)

kruskal_result <- kruskal.test(Bx_to_Prog ~ T4HS_Cat, data = prog_only)
print(kruskal_result)

# ==============================================================================
# ANALYSIS 4: BIOPSY DATE VS DEATH DATE (OS)
# ==============================================================================
os_only <- surv_data %>% filter(!is.na(Bx_to_OS))

bx_os_stats <- os_only %>%
  group_by(T4HS_Cat) %>%
  summarise(
    N = n(),
    Median_Bx_to_OS = median(Bx_to_OS),
    IQR = IQR(Bx_to_OS)
  )
print("--- BIOPSY TO DEATH (OS) METRICS ---")
print(bx_os_stats)

bx_os_box <- ggplot(os_only, aes(x = T4HS_Cat, y = Bx_to_OS, fill = T4HS_Cat)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_viridis_d(option = "C") + 
  theme_minimal() +
  labs(
    title = "Interval: Biopsy Date to Death Date",
    subtitle = "Positive = Bx before Death; Negative = Bx after Death",
    y = "Days Delta",
    x = "Morphologic Category"
  ) +
  theme(legend.position = "none")

print(bx_os_box)

# Save Boxplot
ggsave("Biopsy_to_OS_Boxplot.png", plot = bx_os_box, width = IMG_WIDTH, height = IMG_HEIGHT, dpi = IMG_DPI)

# Statistical Test
kruskal_os <- kruskal.test(Bx_to_OS ~ T4HS_Cat, data = os_only)
print(kruskal_os)

# ==============================================================================
# FINAL STEP: HIPAA COMPLIANT DATA EXPORT
# ==============================================================================

# Updated deidentify_cols to include all calculated delta columns
deidentify_cols <- c("MRN", "Date First Rx", "Date Bx", "Progression", 
                     "Death from disease", "Last FU", "Prog_Num", "Death_Num",
                     "Notes", "T4HS_Cat", "Bx_to_OS", "Bx_to_Prog", "Group_Status")

surv_data_deidentified <- surv_data %>%
  select(colnames(.)[colnames(.) %notin% deidentify_cols])

write_xlsx(surv_data_deidentified, "ST2.CaseRecurrence_Deidentified.xlsx")

print("Analysis complete. Cytotoxic-specific survival figures and stats exported.")
