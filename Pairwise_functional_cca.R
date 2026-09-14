###################################################################
# Pairwise Functional CCA for NRR Curves
#
# Goal:
#   Compare 4 sets / relationships of NRR curves:
#     1. Fundus vs OCT
#     2. Fundus vs Fused
#     3. Fundus vs Landmark registration
#     4. Fused vs Landmark registration
#
# Main question:
#   Which fused curve set is more strongly associated with Fundus?
#
# Main comparison:
#   CCA(Fundus, Fused)
#   CCA(Fundus, Landmark registration)
#
# Added baseline/reference:
#   CCA(Fundus, OCT)
#
# Interpretation:
#   Higher canonical correlations suggest stronger shared variation.
#
# Final plot conventions:
#   Fundus vs OCT                   : black, circle
#   Fundus vs Fused                 : blue, square
#   Fundus vs Landmark registration : red, triangle
#   Fused vs Landmark registration  : green, star
#
# Notes:
#   All plotted CCA curves use solid lines.
#   Legend is placed at bottom left.
###################################################################

#####################################
#####     Set the libraries     #####
#####################################

library.set <- c(
  "fda",
  "pracma",
  "tidyverse"
)

new.packages <- library.set[
  !(library.set %in% installed.packages()[, "Package"])
]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE)
}

lapply(library.set, library, character.only = TRUE)

###################################################################
# 1. Set paths
###################################################################

output_dir <- "/Users"

raw_data_file <- file.path(
  output_dir,
  "OCT_Fundus_NRR.csv"
)

fused_file <- file.path(
  output_dir,
  "NRR_Normal_Norm_FusedCurve.csv"
)

landmark_registration_file <- file.path(
  output_dir,
  "NRR_Normal_Norm_landmark_registered_fused_curves.csv"
)

cca_summary_file <- file.path(
  output_dir,
  "NRR_pairwise_functional_CCA_summary_with_Fundus_OCT_final.csv"
)

cca_correlations_file <- file.path(
  output_dir,
  "NRR_pairwise_functional_CCA_correlations_with_Fundus_OCT_final.csv"
)

cca_decision_file <- file.path(
  output_dir,
  "NRR_pairwise_functional_CCA_decision_with_Fundus_OCT_final.csv"
)

# New file name so it does NOT overwrite previous CCA figures
cca_plot_file <- file.path(
  output_dir,
  "NRR_pairwise_functional_CCA_comparison_with_Fundus_OCT_final_symbols.pdf"
)

###################################################################
# 2. Read files
###################################################################

if(!file.exists(raw_data_file)) {
  stop("Raw OCT/Fundus file not found: ", raw_data_file)
}

if(!file.exists(fused_file)) {
  stop("Fused file not found: ", fused_file)
}

if(!file.exists(landmark_registration_file)) {
  stop("Landmark registration file not found: ", landmark_registration_file)
}

fullOCTdata_normal_NRR <- read.csv(
  raw_data_file,
  stringsAsFactors = FALSE
)

Fused_Data <- read.csv(
  fused_file,
  stringsAsFactors = FALSE
)

Landmark_Registration_Data <- read.csv(
  landmark_registration_file,
  stringsAsFactors = FALSE
)

fullOCTdata_normal_NRR <- unique(fullOCTdata_normal_NRR)
fullOCTdata_normal_NRR <- na.omit(fullOCTdata_normal_NRR)

###################################################################
# 3. Helper functions
###################################################################

get_id_column <- function(data) {
  
  possible_id_cols <- c(
    "PatientID_EYE",
    "Subject_Eye",
    "V1",
    "ID"
  )
  
  id_col <- possible_id_cols[possible_id_cols %in% names(data)][1]
  
  if(is.na(id_col)) {
    stop("Could not find an ID column. Please check metadata column names.")
  }
  
  return(id_col)
}

extract_curve_matrix <- function(data, prefix, expected_ncol = 180) {
  
  curve_cols <- grep(prefix, names(data), value = TRUE)
  
  if(length(curve_cols) != expected_ncol) {
    stop(
      "Expected ",
      expected_ncol,
      " curve columns for prefix ",
      prefix,
      ", but found ",
      length(curve_cols)
    )
  }
  
  curve_mat <- as.matrix(data[, curve_cols])
  storage.mode(curve_mat) <- "numeric"
  
  return(curve_mat)
}

extract_curve_matrix_by_candidates <- function(
    data,
    prefixes,
    expected_ncol = 180,
    label = "curve"
) {
  
  for(prefix in prefixes) {
    
    curve_cols <- grep(prefix, names(data), value = TRUE)
    
    if(length(curve_cols) == expected_ncol) {
      
      curve_mat <- as.matrix(data[, curve_cols])
      storage.mode(curve_mat) <- "numeric"
      
      message(label, " columns detected using prefix: ", prefix)
      
      return(curve_mat)
    }
  }
  
  numeric_cols <- which(sapply(data, is.numeric))
  
  if(length(numeric_cols) >= expected_ncol) {
    
    curve_cols <- names(data)[tail(numeric_cols, expected_ncol)]
    
    curve_mat <- as.matrix(data[, curve_cols])
    storage.mode(curve_mat) <- "numeric"
    
    message(label, " columns detected using last ", expected_ncol, " numeric columns.")
    
    return(curve_mat)
  }
  
  stop("Could not detect ", expected_ncol, " curve columns for ", label)
}

extract_fused_matrix <- function(data) {
  
  extract_curve_matrix_by_candidates(
    data = data,
    prefixes = c(
      "^Fused_NRR_",
      "^Original_Fused_NRR_",
      "^Norm_Fused_NRR_",
      "^NRR_Normal_Norm_FusedCurve_",
      "^FusedCurve_",
      "^NRR_",
      "^V[0-9]+$"
    ),
    expected_ncol = 180,
    label = "Fused"
  )
}

extract_landmark_registration_matrix <- function(data) {
  
  extract_curve_matrix_by_candidates(
    data = data,
    prefixes = c(
      "^Landmark_Fused_NRR_"
    ),
    expected_ncol = 180,
    label = "Landmark registration"
  )
}

normalize.curve <- function(curves, method = "Kron") {
  
  curves <- as.numeric(curves)
  m <- length(curves)
  
  curve_fn <- function(t) curves[t]
  
  X_og_area <- integral(
    curve_fn,
    xmin = 1,
    xmax = m,
    method = method,
    no_intervals = m,
    random = FALSE,
    reltol = 1e-8,
    abstol = 0
  )
  
  X_target <- curves / X_og_area
  
  return(X_target)
}

stratified.curve.norm <- function(temp_OG_curves) {
  
  temp_OG_curves <- as.matrix(temp_OG_curves)
  
  nrow_temp_curves <- nrow(temp_OG_curves)
  ncol_temp_curves <- ncol(temp_OG_curves)
  
  temp_Norm_curves <- matrix(
    NA_real_,
    nrow = nrow_temp_curves,
    ncol = ncol_temp_curves
  )
  
  for(i in seq_len(nrow_temp_curves)) {
    
    temp_Norm_curves[i, ] <- normalize.curve(
      temp_OG_curves[i, ]
    )
  }
  
  return(temp_Norm_curves)
}

###################################################################
# 4. Rebuild OCT and Fundus curves exactly like registration script
###################################################################

NRR_OCT_Normal_OG <- fullOCTdata_normal_NRR[
  ,
  c(c(1:8), grep("OCT", names(fullOCTdata_normal_NRR)))
]

NRR_Fundus_Normal_OG <- fullOCTdata_normal_NRR[
  ,
  c(c(1:8), grep("Fundus", names(fullOCTdata_normal_NRR)))
]

OCT_Norm <- NRR_OCT_Normal_OG
Fundus_Norm <- NRR_Fundus_Normal_OG

OCT_Norm[, 9:ncol(OCT_Norm)] <-
  stratified.curve.norm(
    NRR_OCT_Normal_OG[, -c(1:8)]
  )

Fundus_Norm[, 9:ncol(Fundus_Norm)] <-
  stratified.curve.norm(
    NRR_Fundus_Normal_OG[, -c(1:8)]
  )

d_NRR <- 1:180

NRR_basis_15 <- create.fourier.basis(
  rangeval = c(1, 180),
  nbasis = 15
)

fd_OCT_p15 <- Data2fd(
  argvals = d_NRR,
  y = t(as.matrix(OCT_Norm[, -c(1:8)])),
  basisobj = NRR_basis_15
)

fd_Fundus_p15 <- Data2fd(
  argvals = d_NRR,
  y = t(as.matrix(Fundus_Norm[, -c(1:8)])),
  basisobj = NRR_basis_15
)

OCT_p15_evaluated <- t(
  eval.fd(d_NRR, fd_OCT_p15)
)

Fundus_p15_evaluated <- t(
  eval.fd(d_NRR, fd_Fundus_p15)
)

id_col_raw <- get_id_column(fullOCTdata_normal_NRR)

raw_ids <- fullOCTdata_normal_NRR[[id_col_raw]]
oct_ids <- raw_ids
fundus_ids <- raw_ids

###################################################################
# 5. Extract Fused and Landmark registration matrices
###################################################################

Fused_mat <- extract_fused_matrix(
  Fused_Data
)

Landmark_registration_mat <- extract_landmark_registration_matrix(
  Landmark_Registration_Data
)

id_col_fused <- get_id_column(Fused_Data)
id_col_landmark <- get_id_column(Landmark_Registration_Data)

fused_ids <- Fused_Data[[id_col_fused]]
landmark_ids <- Landmark_Registration_Data[[id_col_landmark]]

###################################################################
# 6. Match subjects across all four datasets
###################################################################

common_ids <- Reduce(
  intersect,
  list(
    oct_ids,
    fundus_ids,
    fused_ids,
    landmark_ids
  )
)

if(length(common_ids) < 5) {
  stop("Too few common subjects across OCT, Fundus, Fused, and Landmark registration datasets.")
}

oct_idx <- match(common_ids, oct_ids)
fundus_idx <- match(common_ids, fundus_ids)
fused_idx <- match(common_ids, fused_ids)
landmark_idx <- match(common_ids, landmark_ids)

OCT_mat <- OCT_p15_evaluated[
  oct_idx,
  ,
  drop = FALSE
]

Fundus_mat <- Fundus_p15_evaluated[
  fundus_idx,
  ,
  drop = FALSE
]

Fused_mat <- Fused_mat[
  fused_idx,
  ,
  drop = FALSE
]

Landmark_registration_mat <- Landmark_registration_mat[
  landmark_idx,
  ,
  drop = FALSE
]

complete_idx <- complete.cases(OCT_mat) &
  complete.cases(Fundus_mat) &
  complete.cases(Fused_mat) &
  complete.cases(Landmark_registration_mat)

OCT_mat <- OCT_mat[complete_idx, , drop = FALSE]
Fundus_mat <- Fundus_mat[complete_idx, , drop = FALSE]
Fused_mat <- Fused_mat[complete_idx, , drop = FALSE]
Landmark_registration_mat <- Landmark_registration_mat[complete_idx, , drop = FALSE]
common_ids <- common_ids[complete_idx]

message("Number of matched complete subjects: ", length(common_ids))

###################################################################
# 7. Convert matrices to functional data objects
###################################################################
# Matrix format:
#   rows = subjects
#   columns = angular grid points
#
# Data2fd requires:
#   rows = grid points
#   columns = subjects
###################################################################

make_fd_from_matrix <- function(
    mat,
    basisobj,
    argvals = 1:180
) {
  
  Data2fd(
    argvals = argvals,
    y = t(mat),
    basisobj = basisobj
  )
}

fd_OCT <- make_fd_from_matrix(
  OCT_mat,
  basisobj = NRR_basis_15,
  argvals = d_NRR
)

fd_Fundus <- make_fd_from_matrix(
  Fundus_mat,
  basisobj = NRR_basis_15,
  argvals = d_NRR
)

fd_Fused <- make_fd_from_matrix(
  Fused_mat,
  basisobj = NRR_basis_15,
  argvals = d_NRR
)

fd_Landmark_Registration <- make_fd_from_matrix(
  Landmark_registration_mat,
  basisobj = NRR_basis_15,
  argvals = d_NRR
)

###################################################################
# 8. Run pairwise functional CCA
###################################################################
# Pairwise comparisons:
#   A. Fundus vs OCT
#   B. Fundus vs Fused
#   C. Fundus vs Landmark registration
#   D. Fused vs Landmark registration
#
# Main comparison:
#   B vs C
#
# Fundus vs OCT is included as a no-fusion baseline/reference.
###################################################################

run_pairwise_cca <- function(
    fd1,
    fd2,
    name1,
    name2,
    ncan = 5,
    lambda = 1e-4
) {
  
  fdPar1 <- fdPar(
    fdobj = fd1,
    Lfdobj = 2,
    lambda = lambda
  )
  
  fdPar2 <- fdPar(
    fdobj = fd2,
    Lfdobj = 2,
    lambda = lambda
  )
  
  cca_fit <- cca.fd(
    fdobj1 = fd1,
    fdobj2 = fd2,
    ncan = ncan,
    ccafdPar1 = fdPar1,
    ccafdPar2 = fdPar2,
    centerfns = TRUE
  )
  
  cca_cor <- NULL
  
  possible_cor_names <- c(
    "ccacorr",
    "corr",
    "corrs",
    "canonical.corr",
    "canonical_correlations"
  )
  
  for(nm in possible_cor_names) {
    if(nm %in% names(cca_fit)) {
      cca_cor <- cca_fit[[nm]]
      break
    }
  }
  
  if(is.null(cca_cor)) {
    stop(
      "Could not automatically find canonical correlations in cca.fd output. ",
      "Run names(cca_fit) to inspect the object."
    )
  }
  
  cca_cor <- as.numeric(cca_cor)
  
  out <- data.frame(
    comparison = paste(name1, "vs", name2),
    curve_set_1 = name1,
    curve_set_2 = name2,
    component = seq_along(cca_cor),
    canonical_correlation = cca_cor
  )
  
  return(
    list(
      fit = cca_fit,
      correlations = out
    )
  )
}

ncan_use <- 5
lambda_use <- 1e-4

cca_Fundus_OCT <- run_pairwise_cca(
  fd1 = fd_Fundus,
  fd2 = fd_OCT,
  name1 = "Fundus",
  name2 = "OCT",
  ncan = ncan_use,
  lambda = lambda_use
)

cca_Fundus_Fused <- run_pairwise_cca(
  fd1 = fd_Fundus,
  fd2 = fd_Fused,
  name1 = "Fundus",
  name2 = "Fused",
  ncan = ncan_use,
  lambda = lambda_use
)

cca_Fundus_Landmark <- run_pairwise_cca(
  fd1 = fd_Fundus,
  fd2 = fd_Landmark_Registration,
  name1 = "Fundus",
  name2 = "Landmark registration",
  ncan = ncan_use,
  lambda = lambda_use
)

cca_Fused_Landmark <- run_pairwise_cca(
  fd1 = fd_Fused,
  fd2 = fd_Landmark_Registration,
  name1 = "Fused",
  name2 = "Landmark registration",
  ncan = ncan_use,
  lambda = lambda_use
)

CCA_All_Correlations <- bind_rows(
  cca_Fundus_OCT$correlations,
  cca_Fundus_Fused$correlations,
  cca_Fundus_Landmark$correlations,
  cca_Fused_Landmark$correlations
)

print(CCA_All_Correlations)

###################################################################
# 9. Summarize which curve set is better relative to Fundus
###################################################################

summarize_cca <- function(cor_df) {
  
  cor_df %>%
    group_by(comparison) %>%
    summarise(
      first_canonical_correlation = canonical_correlation[component == 1][1],
      mean_first_3 = mean(canonical_correlation[component <= 3], na.rm = TRUE),
      mean_all = mean(canonical_correlation, na.rm = TRUE),
      min_all = min(canonical_correlation, na.rm = TRUE),
      .groups = "drop"
    )
}

CCA_Summary <- summarize_cca(
  CCA_All_Correlations
)

print(CCA_Summary)

fundus_fused <- CCA_Summary %>%
  filter(comparison == "Fundus vs Fused")

fundus_landmark <- CCA_Summary %>%
  filter(comparison == "Fundus vs Landmark registration")

comparison_decision <- data.frame(
  criterion = c(
    "First canonical correlation",
    "Mean of first 3 canonical correlations",
    "Mean of all reported canonical correlations"
  ),
  fused = c(
    fundus_fused$first_canonical_correlation,
    fundus_fused$mean_first_3,
    fundus_fused$mean_all
  ),
  landmark_registration = c(
    fundus_landmark$first_canonical_correlation,
    fundus_landmark$mean_first_3,
    fundus_landmark$mean_all
  )
)

comparison_decision$better_method <- ifelse(
  comparison_decision$landmark_registration >
    comparison_decision$fused,
  "Landmark registration",
  "Fused"
)

print(comparison_decision)

###################################################################
# 10. Save CCA outputs
###################################################################

write.csv(
  CCA_All_Correlations,
  cca_correlations_file,
  row.names = FALSE
)

write.csv(
  CCA_Summary,
  cca_summary_file,
  row.names = FALSE
)

write.csv(
  comparison_decision,
  cca_decision_file,
  row.names = FALSE
)

message("Saved CCA correlations to: ", cca_correlations_file)
message("Saved CCA summary to: ", cca_summary_file)
message("Saved CCA decision table to: ", cca_decision_file)

###################################################################
# 11. Plot canonical correlations
###################################################################
# Final requested plot styling:
#   a) Different point symbols:
#        Fundus vs OCT                   = circle
#        Fundus vs Fused                 = square
#        Fundus vs Landmark registration = triangle
#        Fused vs Landmark registration  = star
#
#   b) Solid lines for all comparisons
#
#   c) Color mapping:
#        Fundus vs OCT                   = black
#        Fundus vs Fused                 = blue
#        Fundus vs Landmark registration = red
#        Fused vs Landmark registration  = green
#
#   d) Legend at bottom left
###################################################################

cols <- c(
  "Fundus vs OCT" = "black",
  "Fundus vs Fused" = "blue",
  "Fundus vs Landmark registration" = "red3",
  "Fused vs Landmark registration" = "darkgreen"
)

ltys <- c(
  "Fundus vs OCT" = 1,
  "Fundus vs Fused" = 1,
  "Fundus vs Landmark registration" = 1,
  "Fused vs Landmark registration" = 1
)

pchs <- c(
  "Fundus vs OCT" = 19,                    # circle
  "Fundus vs Fused" = 15,                  # square
  "Fundus vs Landmark registration" = 17,  # triangle
  "Fused vs Landmark registration" = 8     # star
)

pdf(
  cca_plot_file,
  width = 9,
  height = 6
)

par(mar = c(5, 5, 4, 2))

plot(
  NULL,
  xlim = c(1, ncan_use),
  ylim = c(0, 1),
  xlab = "Canonical component",
  ylab = "Canonical correlation",
  main = "",
  cex.lab = 2,
  cex.axis = 1.5
)

abline(
  h = 0.5,
  col = "gray50",
  lty = 2,
  lwd = 1.5
)

for(comp_name in names(cols)) {
  
  temp <- CCA_All_Correlations[
    CCA_All_Correlations$comparison == comp_name,
  ]
  
  if(nrow(temp) > 0) {
    
    lines(
      temp$component,
      temp$canonical_correlation,
      type = "b",
      pch = pchs[comp_name],
      lwd = 2,
      cex = 1.5,
      col = cols[comp_name],
      lty = ltys[comp_name]
    )
  }
}

legend(
  "bottomleft",
  legend = names(cols),
  col = cols,
  lty = ltys,
  lwd = rep(2, length(cols)),
  pch = pchs,
  bty = "n",
  cex = 1.05
)

dev.off()

message("Saved CCA comparison plot to: ", cca_plot_file)

###################################################################
# 12. Optional PNG version
###################################################################

cca_png_file <- sub("\\.pdf$", ".png", cca_plot_file)

png(
  cca_png_file,
  width = 2700,
  height = 1800,
  res = 300
)

par(mar = c(5, 5, 4, 2))

plot(
  NULL,
  xlim = c(1, ncan_use),
  ylim = c(0, 1),
  xlab = "Canonical component",
  ylab = "Canonical correlation",
  main = "",
  cex.lab = 2,
  cex.axis = 1.5
)

abline(
  h = 0.5,
  col = "gray50",
  lty = 2,
  lwd = 1.5
)

for(comp_name in names(cols)) {
  
  temp <- CCA_All_Correlations[
    CCA_All_Correlations$comparison == comp_name,
  ]
  
  if(nrow(temp) > 0) {
    
    lines(
      temp$component,
      temp$canonical_correlation,
      type = "b",
      pch = pchs[comp_name],
      lwd = 2,
      cex = 1.5,
      col = cols[comp_name],
      lty = ltys[comp_name]
    )
  }
}

legend(
  "bottomleft",
  legend = names(cols),
  col = cols,
  lty = ltys,
  lwd = rep(2, length(cols)),
  pch = pchs,
  bty = "n",
  cex = 1.05
)

dev.off()

message("Saved PNG CCA comparison plot to: ", cca_png_file)

###################################################################
# 13. Optional: inspect CCA object structure
###################################################################
# If the code errors while extracting canonical correlations,
# run one of these:
#
# names(cca_Fundus_OCT$fit)
# names(cca_Fundus_Fused$fit)
# str(cca_Fundus_OCT$fit)
# str(cca_Fundus_Fused$fit)
#
# The correlations are usually stored in:
#   cca_fit$ccacorr
#
###################################################################

###################################################################
# End of Pairwise Functional CCA Script
###################################################################