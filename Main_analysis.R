################################################################################
# Main_analysis.R
# Complete Post-Docking Analytics Pipeline for Vina + SwissADME + Efficiency Metrics
# Cleaned, unified, executable, GitHub-ready
################################################################################

############################
# 1) INITIAL CONFIGURATION #
############################

VINA_WORKSPACE <- path.expand("~/Vina_Workspace")
adme_file <- file.path(VINA_WORKSPACE, "compounds", "adme_features_for_efficiency.csv")

if (!file.exists(adme_file)) {
  warning("ADME file not found at: ", adme_file)
}

receptors_dir <- file.path(VINA_WORKSPACE, "receptors")
if (!dir.exists(receptors_dir)) stop("Receptors folder not found.")

# Auto-detect proteins
protein_dirs <- list.dirs(receptors_dir, recursive = FALSE, full.names = FALSE)
protein_dirs <- Filter(function(x) {
  d <- file.path(receptors_dir, x, "dockings")
  dir.exists(d) && length(list.files(d, "\\.log$")) > 0
}, protein_dirs)

if (!length(protein_dirs)) stop("No proteins with .log docking files found.")

cat("Detected proteins:\n"); print(protein_dirs)

###########################
# 2) REQUIRED PACKAGES    #
###########################

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(plotly)
library(viridis)
library(grid)
library(htmlwidgets)
library(readr)
library(janitor)
library(rcdk)
library(forcats)
library(scales)

#############################################################
# 3) UNIVERSAL FUNCTION TO PARSE AUTODOCK VINA LOG FILES    #
#############################################################

process_vina_log <- function(fpath) {
  name <- tools::file_path_sans_ext(basename(fpath))
  m <- str_match(name, "^(.*)_pocket(\\d+)$")

  if (!is.na(m[1, 1])) {
    prefix <- m[1, 2]
    parts  <- str_split(prefix, "_")[[1]]
    drug   <- parts[length(parts)]
    protein <- paste(parts[-length(parts)], collapse = "_")
    pocket  <- paste0("pocket", m[1, 3])
  } else {
    parts <- str_split(name, "_")[[1]]
    protein <- parts[1]
    drug    <- ifelse(length(parts) >= 2, parts[2], NA)
    pocket  <- ifelse(length(parts) >= 3, parts[3], NA)
  }

  lines <- readLines(fpath, warn = FALSE)
  idx <- grep("^\\s*1\\s+", lines)[1]

  affinity <- if (!is.na(idx))
    as.numeric(str_split(str_trim(lines[idx]), "\\s+")[[1]][2])
  else NA_real_

  tibble(protein, drug, pocket, affinity)
}

#############################################################
# 4) MAIN LOOP — PROCESS EACH PROTEIN INDEPENDENTLY         #
#############################################################

for (PROTEIN_NAME in protein_dirs) {

  cat("\n======================================\n")
  cat(" Processing:", PROTEIN_NAME, "\n")
  cat("======================================\n")

  log_dir <- file.path(receptors_dir, PROTEIN_NAME, "dockings")
  setwd(file.path(receptors_dir, PROTEIN_NAME))

  files <- list.files(log_dir, "\\.log$", full.names = TRUE)
  df_logs <- map_df(files, process_vina_log)

  ###############################################
  # 5) BEST AFFINITY PER PROTEIN × LIGAND       #
  ###############################################

  df_min <- df_logs %>%
    group_by(protein, drug) %>%
    summarise(min_affinity = min(affinity, na.rm = TRUE), .groups = "drop")

  df_filtered <- df_min %>%
    filter(min_affinity > -20 & min_affinity < 0)

  ###############################################
  # 6) HEATMAPS (ComplexHeatmap)                #
  ###############################################

  mat <- df_filtered %>%
    pivot_wider(names_from = drug, values_from = min_affinity) %>%
    column_to_rownames("protein") %>%
    as.matrix()

  col_fun <- colorRamp2(
    c(min(mat, na.rm=TRUE), mean(mat, na.rm=TRUE), max(mat, na.rm=TRUE)),
    c("navy", "white", "firebrick")
  )

  ht1 <- Heatmap(
    mat,
    name = "Affinity (kcal/mol)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_names_rot = 45
  )

  png(paste0("Heatmap_Protein_vs_Compounds_", PROTEIN_NAME, ".png"),
      width=4000, height=2000, res=150)
  draw(ht1)
  dev.off()

  ###############################################
  # 7) 3D AFFINITY LANDSCAPE (plotly)           #
  ###############################################

  df_plot <- df_logs %>%
    mutate(
      pocket_num = parse_number(pocket),
      protein_f = factor(protein)
    ) %>%
    filter(!is.na(affinity))

  p3d <- plot_ly(
    df_plot,
    x = ~protein_f,
    y = ~pocket_num,
    z = ~affinity,
    color = ~drug,
    colors = viridis(length(unique(df_plot$drug))),
    type = "scatter3d",
    mode = "markers"
  ) %>%
  layout(scene=list(
    xaxis=list(title="Protein"),
    yaxis=list(title="Pocket"),
    zaxis=list(title="Affinity (kcal/mol)")
  ))

  saveWidget(p3d, paste0("Affinity_3D_", PROTEIN_NAME, ".html"), selfcontained=TRUE)

  ###############################################
  # 8) ADME IMPORT AND NORMALIZATION            #
  ###############################################

  adme <- read.csv(adme_file, stringsAsFactors = FALSE) %>% clean_names()

  norm_key <- function(x)
    x |> tolower() |> stringi::stri_trans_general("Latin-ASCII") |> gsub("[^a-z0-9]", "", x)

  adme <- adme %>% mutate(drug_key = norm_key(drug))
  df_min <- df_min %>% mutate(drug_key = norm_key(drug))

  adme_u <- adme %>% group_by(drug_key) %>%
    summarise(across(everything(), ~if(is.numeric(.x)) median(.x, na.rm=TRUE) else first(.x)),
              .groups="drop")

  df_eff <- df_min %>%
    left_join(adme_u, by="drug_key")

  ###################################################
  # 9) LIGAND EFFICIENCY METRICS (LE, LLE, FQ, pKd) #
  ###################################################

  R <- 1.9872036e-3
  T <- 310.15

  df_eff <- df_eff %>%
    mutate(
      Kd  = exp(min_affinity/(R*T)),
      pKd = -log10(Kd),
      LE  = ifelse(hac > 0, -min_affinity / hac, NA),
      logP_sel = coalesce(consensus_log_p, xlogp3, wlogp, mlogp, i_logp),
      LLE = pKd - logP_sel
    )

  fq_den <- function(ha) 0.0715 + 7.5328/ha + 25.7079/ha^2 - 361.4722/ha^3
  df_eff <- df_eff %>% mutate(FQ = ifelse(hac > 0, (pKd/hac)/fq_den(hac), NA))

  write.csv(df_eff, paste0("Dock_Efficiency_", PROTEIN_NAME, ".csv"), row.names=FALSE)

  ###############################################
  # 10) BOILED-EGG MODEL                        #
  ###############################################

  adme_u <- adme_u %>%
    mutate(
      egg_logp = coalesce(wlogp, logP_sel),
      egg_hia  = (tpsa <= 131.6) & between(egg_logp, -0.7, 6),
      egg_bbb  = (tpsa <= 90)    & between(egg_logp, -0.7, 6)
    )

  ###############################################
  # 11) SUMMARY PLOTS                           #
  ###############################################

  p_egg <- ggplot(adme_u, aes(x=egg_logp, y=tpsa)) +
    geom_point(aes(color=egg_bbb, shape=egg_hia), size=2) +
    scale_color_manual(values=c("FALSE"="blue","TRUE"="red")) +
    theme_minimal() +
    labs(title="BOILED-Egg Model", x="logP", y="TPSA")

  ggsave(paste0("Boiled_Egg_", PROTEIN_NAME, ".png"), p_egg, width=6, height=6)

  ###############################################
  # 12) TOP 20 BEST LIGANDS                     #
  ###############################################

  top20 <- df_eff %>% arrange(desc(LLE)) %>% slice_head(n=20)
  write.csv(top20, paste0("Top20_", PROTEIN_NAME, ".csv"), row.names=FALSE)

  cat("\nFinished:", PROTEIN_NAME, "\n")
}

cat("\n======================================\n")
cat(" All proteins processed successfully!\n")
cat("======================================\n")
