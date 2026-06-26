# ============================================================
# Figure 11A. Opportunity decay after CGP
# 20260626
# Masaki Ishida
# Risk table: aligned under X-axis ticks
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(lubridate)
})

setwd("~/C-CAT_Drug Loss_Analysis")

fig_dir <- "figures/Figure11_OpportunityDecay"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

admin_censor_date <- as.Date("2025-12-31")
font_family <- "Arial"

landmark_days <- c(365)
landmark_labels <- c("12-month")
risk_months <- c(0, 12, 24, 36, 48)
risk_days <- round(risk_months * 365.25 / 12)

id_col <- "C-CAT調査結果-基本項目.ハッシュID"
panel_col <- "C-CAT調査結果-基本項目.パネル名"
cgp_date_col <- "C-CAT調査結果-基本項目.検体採取日"
domestic_approved_n <- "C-CAT調査結果-サマリ.国内承認薬数"
report_age_col <- "C-CAT調査結果-基本項目.年齢(年)"

first_token <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL", "nan", "不明", "不詳")] <- NA_character_
  ifelse(grepl(";", x), sub(";.*$", "", x), x)
}

parse_date <- function(x) {
  suppressWarnings(as.Date(parse_date_time(
    first_token(x),
    orders = c("Y/m/d", "Y-m-d"),
    quiet = TRUE
  )))
}

is_solid <- function(x) {
  x <- as.character(x)
  !is.na(x) & !grepl(
    "血液|リンパ|骨髄|白血病|リンパ腫|Hematologic|Lymphoid|Myeloid",
    x,
    ignore.case = TRUE
  )
}

fmt_rate_ci <- function(est, lo, hi) {
  sprintf("%.1f%% (%.1f-%.1f)", 100 * est, 100 * lo, 100 * hi)
}

theme_nm <- function() {
  theme_classic(base_size = 11, base_family = font_family) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.title = element_text(face = "plain", colour = "black", size = 11),
      axis.text = element_text(colour = "black", size = 10),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      plot.margin = margin(5, 10, 5, 5)
    )
}

build_survival_fields <- function(dt) {
  dt[, death_dt := as.Date(parse_date(death_date))]
  dt[, fu_dt := as.Date(parse_date(last_fu))]
  dt[, matched_tx_dt := as.Date(parse_date(matched_tx_date))]

  dt[, fu_end := fcase(
    outcome == "死亡" & !is.na(death_dt),
    pmin(death_dt, admin_censor_date, na.rm = TRUE),

    !is.na(fu_dt),
    pmin(fu_dt, admin_censor_date, na.rm = TRUE),

    default = as.Date(NA)
  )]

  dt[, tx_after_cgp := fifelse(
    !is.na(matched_tx_dt) & matched_tx_dt >= cgp_date,
    matched_tx_dt,
    as.Date(NA)
  )]

  dt[, death_after_cgp := fifelse(
    outcome == "死亡" & !is.na(death_dt) & death_dt >= cgp_date,
    death_dt,
    as.Date(NA)
  )]

  dt[, stop_dt := pmin(
    fifelse(!is.na(tx_after_cgp), tx_after_cgp, as.Date("9999-12-31")),
    fifelse(!is.na(death_after_cgp), death_after_cgp, as.Date("9999-12-31")),
    fu_end,
    admin_censor_date,
    na.rm = TRUE
  )]

  dt[, opp_event := as.integer(
    !is.na(death_after_cgp) &
      (is.na(tx_after_cgp) | death_after_cgp < tx_after_cgp)
  )]

  dt[, opp_time_days := as.numeric(stop_dt - cgp_date)]
  dt[opp_time_days < 0, opp_time_days := 0]

  dt
}

apply_exclusions <- function(dt) {
  n0 <- nrow(dt)

  dt[, flag_tx_before_cgp := !is.na(matched_tx_dt) & matched_tx_dt < cgp_date]
  dt[, flag_death_before_cgp := !is.na(death_dt) & death_dt < cgp_date]

  n_tx <- dt[flag_tx_before_cgp == TRUE, .N]
  n_death <- dt[flag_death_before_cgp == TRUE, .N]

  out <- dt[flag_tx_before_cgp == FALSE & flag_death_before_cgp == FALSE]

  n_miss_cgp <- out[is.na(cgp_date), .N]
  n_miss_fu <- out[is.na(fu_end), .N]

  out <- out[!is.na(cgp_date) & !is.na(fu_end)]

  message("Starting N: ", n0)
  message("Excluded matched therapy before CGP: ", n_tx)
  message("Excluded death before CGP: ", n_death)
  message("Excluded missing CGP date: ", n_miss_cgp)
  message("Excluded missing follow-up: ", n_miss_fu)
  message("Analyzable N: ", nrow(out))

  out
}

surv_step_df <- function(fit) {
  sm <- summary(fit)
  data.table(
    time_days = sm$time,
    time_months = sm$time / 30.4375,
    surv = sm$surv,
    lower = sm$lower,
    upper = sm$upper
  )
}

km_landmarks <- function(fit) {
  ss <- summary(fit, times = landmark_days, extend = TRUE)

  data.table(
    landmark = landmark_labels,
    time_days = landmark_days,
    rate = ss$surv,
    lower_95 = ss$lower,
    upper_95 = ss$upper,
    rate_ci = fmt_rate_ci(ss$surv, ss$lower, ss$upper)
  )
}

plot_km_panel_a <- function(fit, annotation_txt = NULL, x_max_months = 48) {
  sdf <- surv_step_df(fit)

  sdf[, surv_pct := 100 * surv]
  sdf[, lower_pct := 100 * lower]
  sdf[, upper_pct := 100 * upper]

  p_main <- ggplot(sdf, aes(x = time_months, y = surv_pct)) +
    geom_ribbon(
      aes(ymin = lower_pct, ymax = upper_pct),
      fill = "grey85",
      alpha = 0.65,
      colour = NA
    ) +
    geom_step(linewidth = 1.05, colour = "black") +
    scale_x_continuous(
      limits = c(0, x_max_months),
      breaks = risk_months,
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = c(0.01, 0)
    ) +
    labs(
      x = "Months from CGP",
      y = "Remaining genomic opportunity (%)"
    ) +
    theme_nm()

  if (!is.null(annotation_txt)) {
    p_main <- p_main +
      annotate(
        "text",
        x = 47.5,
        y = 98,
        hjust = 1,
        vjust = 1,
        label = annotation_txt,
        size = 3.1,
        lineheight = 0.95,
        family = font_family
      )
  }

  rs <- summary(fit, times = risk_days, extend = TRUE)

  risk_tbl <- data.table(
    months = risk_months,
    n_risk = rs$n.risk
  )

  risk_plot <- ggplot(risk_tbl, aes(x = months)) +
    geom_text(aes(y = 1, label = months), size = 3.0, family = font_family) +
    geom_text(aes(y = 0, label = format(n_risk, big.mark = ",")), size = 3.0, family = font_family) +
    annotate("text", x = -3.0, y = 1, label = "Months from CGP", hjust = 0, size = 3.0, family = font_family) +
    annotate("text", x = -3.0, y = 0, label = "No. at risk", hjust = 0, size = 3.0, family = font_family) +
    scale_x_continuous(
      limits = c(-6, x_max_months),
      breaks = risk_months,
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(-0.4, 1.4)) +
    theme_void(base_family = font_family) +
    theme(plot.margin = margin(0, 10, 0, 5))

  plot_grid(
    p_main,
    risk_plot,
    ncol = 1,
    rel_heights = c(1, 0.16),
    align = "v"
  )
}

# ============================================================
# Load data
# ============================================================

patient_case <- readRDS("output/patient_case.rds")
setDT(patient_case)

report_all <- readRDS("output/report_all.rds")
setDT(report_all)

report_cols <- c(id_col, panel_col, cgp_date_col, domestic_approved_n, report_age_col)
report_cols <- report_cols[report_cols %in% names(report_all)]
report_all <- report_all[, ..report_cols]

panels <- report_all[, .(
  panel_name = {
    x <- unique(na.omit(as.character(get(panel_col))))
    if (length(x) == 0) NA_character_ else x[1]
  },
  cgp_date_report = {
    x <- unique(na.omit(as.character(get(cgp_date_col))))
    if (length(x) == 0) NA_character_ else x[1]
  },
  age_report = suppressWarnings(as.numeric(get(report_age_col))[1])
), by = id_col]

actionable <- report_all[
  ,
  .(max_domestic_approved = suppressWarnings(max(as.numeric(get(domestic_approved_n)), na.rm = TRUE))),
  by = id_col
]

actionable[is.infinite(max_domestic_approved), max_domestic_approved := 0]
actionable[, actionable := max_domestic_approved >= 1]

rm(report_all)
gc()

cache_tx <- "output/matched_tx_dates_actionable_F1CDx.rds"

if (!file.exists(cache_tx)) {
  stop("Missing output/matched_tx_dates_actionable_F1CDx.rds. Run Figure4 opportunity decay script first.")
}

tx_dates <- readRDS(cache_tx)
setDT(tx_dates)

patient <- merge(patient_case, panels, by = id_col, all.x = TRUE)
patient <- merge(patient, actionable[, c(id_col, "actionable"), with = FALSE], by = id_col, all.x = TRUE)

patient[, hash_id := get(id_col)]

patient <- merge(patient, tx_dates, by = "hash_id", all.x = TRUE)

patient[, f1cdx := grepl("FoundationOne.*CDx", panel_name, ignore.case = TRUE)]
patient[, solid := is_solid(cancer_l1)]
patient[, cgp_date := parse_date(cgp_date_report)]
patient[is.na(cgp_date) & !is.na(cgp_date_case), cgp_date := parse_date(cgp_date_case)]
patient[is.na(actionable), actionable := FALSE]

base <- patient[f1cdx == TRUE & solid == TRUE & actionable == TRUE]

base <- build_survival_fields(base)
analysis <- apply_exclusions(copy(base))

message("Events: ", analysis[, sum(opp_event)])
message("Censored: ", analysis[, .N - sum(opp_event)])
message("Matched therapy censored: ", analysis[opp_event == 0 & !is.na(tx_after_cgp), .N])

fit_a <- survfit(Surv(opp_time_days, opp_event) ~ 1, data = analysis, conf.int = 0.95)

lm_a <- km_landmarks(fit_a)

anno_a <- paste(
  "Median opportunity: not reached",
  paste0(lm_a$landmark, ": ", lm_a$rate_ci),
  sep = "\n"
)

p_a <- plot_km_panel_a(fit_a, annotation_txt = anno_a)

print(p_a)

ggsave(
  file.path(fig_dir, "Figure11A_OpportunityDecay_Overall_Months.pdf"),
  p_a,
  width = 7.2,
  height = 5.8,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Figure11A_OpportunityDecay_Overall_Months.png"),
  p_a,
  width = 7.2,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

fwrite(lm_a, file.path(fig_dir, "Figure11A_summary_table_months.csv"))

message("Figure 11A completed.")
message("Saved to: ", fig_dir)
