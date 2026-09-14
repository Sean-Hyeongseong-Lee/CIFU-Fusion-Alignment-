###################################################################
# Module: Landmark Registration for OCT-Fundus NRR Curves
# Dataset: Non-glaucoma eyes only
#
# Goal:
#   For each eye:
#     1. Use the OCT curve as the target/reference curve
#     2. Detect two OCT landmarks and two Fundus landmarks
#     3. Register the Fundus curve to the OCT curve using landmark registration
#     4. Average OCT + registered Fundus to generate one fused curve
#
# This script uses direct landmark-based warping rather than time_warping()
# or fda::landmarkreg().
###################################################################

#####################################
#####     Set the libraries     #####
#####################################
library.set <- c("fda", "pracma")

new.packages <- library.set[
  !(library.set %in% installed.packages()[, "Package"])
]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE)
}

lapply(library.set, library, character.only = TRUE)

###################################################################
# 1. Read in the data
###################################################################

localPath_normal_NRR <- "/Users/OCT_Fundus_NRR.csv"

fullOCTdata_normal_NRR <- read.csv(
  localPath_normal_NRR,
  stringsAsFactors = TRUE
)

# Remove duplicate rows and rows with missing values
fullOCTdata_normal_NRR <- unique(fullOCTdata_normal_NRR)
fullOCTdata_normal_NRR <- na.omit(fullOCTdata_normal_NRR)

###################################################################
# 2. Parse the data into OCT and Fundus
###################################################################

NRR_OCT_Normal_OG <- fullOCTdata_normal_NRR[
  , c(c(1:8), grep("OCT", x = names(fullOCTdata_normal_NRR)))
]

NRR_Fundus_Normal_OG <- fullOCTdata_normal_NRR[
  , c(c(1:8), grep("Fundus", x = names(fullOCTdata_normal_NRR)))
]

###################################################################
# 3. Normalize the curves so area under each curve = 1
###################################################################

normalize.curve <- function(curves, method = "Kron") {
  
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
  
  nrow_temp_curves <- dim(temp_OG_curves)[1]
  ncol_temp_curves <- dim(temp_OG_curves)[2]
  
  temp_Norm_curves <- matrix(
    NA,
    nrow = nrow_temp_curves,
    ncol = ncol_temp_curves,
    byrow = TRUE
  )
  
  for(i in 1:nrow_temp_curves) {
    
    temp_Norm_curves_i <- temp_OG_curves[i, ]
    
    temp_Norm_curves[i, ] <- normalize.curve(
      t(temp_Norm_curves_i)
    )
  }
  
  temp_Norm_curves <- apply(temp_Norm_curves, 2, as.numeric)
  
  return(temp_Norm_curves)
}

# Empty matrices for normalized curves
NRR_OCT_Normal_Norm <- matrix(
  NA,
  nrow = nrow(NRR_OCT_Normal_OG),
  ncol = ncol(NRR_OCT_Normal_OG),
  byrow = TRUE
)

NRR_Fundus_Normal_Norm <- matrix(
  NA,
  nrow = nrow(NRR_Fundus_Normal_OG),
  ncol = ncol(NRR_Fundus_Normal_OG),
  byrow = TRUE
)

# Normalize only the 180 curve values
NRR_OCT_Normal_Norm[, 9:ncol(NRR_OCT_Normal_OG)] <-
  stratified.curve.norm(NRR_OCT_Normal_OG[, -c(1:8)])

NRR_Fundus_Normal_Norm[, 9:ncol(NRR_Fundus_Normal_OG)] <-
  stratified.curve.norm(NRR_Fundus_Normal_OG[, -c(1:8)])

# Convert to data frames and restore metadata
NRR_OCT_Normal_Norm <- as.data.frame(NRR_OCT_Normal_Norm)
NRR_Fundus_Normal_Norm <- as.data.frame(NRR_Fundus_Normal_Norm)

NRR_OCT_Normal_Norm[, 1:8] <- fullOCTdata_normal_NRR[, 1:8]
NRR_Fundus_Normal_Norm[, 1:8] <- fullOCTdata_normal_NRR[, 1:8]

###################################################################
# 4. Create p = 15 Fourier-smoothed functional curves
###################################################################

d_NRR <- seq(from = 1, to = 180)

NRR_basis_15 <- create.fourier.basis(
  rangeval = c(1, 180),
  nbasis = 15
)

fdNRR_OCT_Normal_Norm_p15 <- Data2fd(
  argvals = d_NRR,
  y = t(NRR_OCT_Normal_Norm[, -c(1:8)]),
  basisobj = NRR_basis_15
)

fdNRR_Fundus_Normal_Norm_p15 <- Data2fd(
  argvals = d_NRR,
  y = t(NRR_Fundus_Normal_Norm[, -c(1:8)]),
  basisobj = NRR_basis_15
)

# Evaluate smoothed curves back on the 180-point grid
NRR_OCT_Normal_Norm_p15_evaluated <- t(
  eval.fd(d_NRR, fdNRR_OCT_Normal_Norm_p15)
)

NRR_Fundus_Normal_Norm_p15_evaluated <- t(
  eval.fd(d_NRR, fdNRR_Fundus_Normal_Norm_p15)
)

###################################################################
# 5. Detect landmarks
###################################################################
# We use:
#   Landmark 1 = dominant peak in first half of the curve
#   Landmark 2 = dominant peak in second half of the curve
#
# The OCT landmarks are the target landmarks.
# The Fundus curve will be warped so its landmarks align to OCT landmarks.
###################################################################

detect.two.landmarks <- function(curve) {
  
  landmark_1 <- which.max(curve[2:90]) + 1
  landmark_2 <- which.max(curve[91:179]) + 90
  
  return(c(landmark_1, landmark_2))
}

OCT_landmarks <- t(
  apply(
    NRR_OCT_Normal_Norm_p15_evaluated,
    1,
    detect.two.landmarks
  )
)

Fundus_landmarks <- t(
  apply(
    NRR_Fundus_Normal_Norm_p15_evaluated,
    1,
    detect.two.landmarks
  )
)

colnames(OCT_landmarks) <- c(
  "OCT_landmark_1",
  "OCT_landmark_2"
)

colnames(Fundus_landmarks) <- c(
  "Fundus_landmark_1",
  "Fundus_landmark_2"
)

###################################################################
# 6. Landmark-registration warping function
###################################################################
# For each eye:
#
#   target grid = OCT time grid
#   source grid = Fundus time grid
#
# We construct a monotone piecewise-linear warping function h(t)
# such that:
#
#   h(1)                    = 1
#   h(OCT landmark 1)       = Fundus landmark 1
#   h(OCT landmark 2)       = Fundus landmark 2
#   h(180)                  = 180
#
# Then:
#
#   registered_Fundus(t) = Fundus( h(t) )
#
# This makes the registered Fundus landmarks occur at the OCT landmark
# locations.
###################################################################

landmark.register.one.curve <- function(
    fundus_curve,
    oct_landmarks,
    fundus_landmarks,
    grid = 1:180
) {
  
  # Target landmark positions: OCT
  target_times <- c(
    1,
    oct_landmarks[1],
    oct_landmarks[2],
    180
  )
  
  # Source landmark positions: Fundus
  source_times <- c(
    1,
    fundus_landmarks[1],
    fundus_landmarks[2],
    180
  )
  
  # Safety check: landmark times must remain ordered
  if(
    any(diff(target_times) <= 0) ||
    any(diff(source_times) <= 0)
  ) {
    stop("Landmarks are not strictly increasing.")
  }
  
  # Piecewise-linear monotone warping function:
  # maps OCT-target grid locations -> Fundus-source grid locations
  warped_grid <- approx(
    x = target_times,
    y = source_times,
    xout = grid,
    method = "linear",
    ties = "ordered"
  )$y
  
  # Evaluate the Fundus curve at the warped positions
  registered_fundus <- approx(
    x = grid,
    y = fundus_curve,
    xout = warped_grid,
    method = "linear",
    rule = 2
  )$y
  
  return(
    list(
      registered_fundus = registered_fundus,
      warped_grid = warped_grid
    )
  )
}

###################################################################
# 7. Register each Fundus curve to its paired OCT curve
###################################################################

n_eye <- nrow(fullOCTdata_normal_NRR)

Fundus_registered <- matrix(
  NA_real_,
  nrow = n_eye,
  ncol = length(d_NRR)
)

Fused_landmark_registered <- matrix(
  NA_real_,
  nrow = n_eye,
  ncol = length(d_NRR)
)

Warped_grids <- matrix(
  NA_real_,
  nrow = n_eye,
  ncol = length(d_NRR)
)

registration_status <- rep("ok", n_eye)

for(i in seq_len(n_eye)) {
  
  tryCatch({
    
    reg_i <- landmark.register.one.curve(
      fundus_curve = NRR_Fundus_Normal_Norm_p15_evaluated[i, ],
      oct_landmarks = OCT_landmarks[i, ],
      fundus_landmarks = Fundus_landmarks[i, ],
      grid = d_NRR
    )
    
    Fundus_registered[i, ] <- reg_i$registered_fundus
    
    Warped_grids[i, ] <- reg_i$warped_grid
    
    Fused_landmark_registered[i, ] <- 0.5 * (
      NRR_OCT_Normal_Norm_p15_evaluated[i, ] +
        Fundus_registered[i, ]
    )
    
  }, error = function(e) {
    
    registration_status[i] <<- paste0(
      "failed: ",
      conditionMessage(e)
    )
    
  })
  
  if(i %% 50 == 0) {
    message("Processed ", i, " of ", n_eye, " eyes")
  }
}

###################################################################
# 8. Check registration success
###################################################################

keep <- registration_status == "ok"

message("Successful registrations: ", sum(keep), " / ", n_eye)
message("Failed registrations: ", sum(!keep), " / ", n_eye)

print(table(registration_status))

###################################################################
# 9. Build output data frames
###################################################################

########################################
# 9A. Landmark-registered fused curves
########################################

Fused_Landmark_Registered_Data <- data.frame(
  fullOCTdata_normal_NRR[keep, 1:8],
  Fused_landmark_registered[keep, , drop = FALSE]
)

names(Fused_Landmark_Registered_Data)[
  9:ncol(Fused_Landmark_Registered_Data)
] <- paste0("Landmark_Fused_NRR_", 0:179)

########################################
# 9B. Registered Fundus curves
########################################

Registered_Fundus_Data <- data.frame(
  fullOCTdata_normal_NRR[keep, 1:8],
  Fundus_registered[keep, , drop = FALSE]
)

names(Registered_Fundus_Data)[
  9:ncol(Registered_Fundus_Data)
] <- paste0("Landmark_Registered_Fundus_NRR_", 0:179)

########################################
# 9C. Warping grids
########################################

Warping_Grid_Data <- data.frame(
  fullOCTdata_normal_NRR[keep, 1:8],
  Warped_grids[keep, , drop = FALSE]
)

names(Warping_Grid_Data)[
  9:ncol(Warping_Grid_Data)
] <- paste0("Warped_Grid_", 0:179)

########################################
# 9D. QC file
########################################

Landmark_Registration_QC <- data.frame(
  fullOCTdata_normal_NRR[, 1:8],
  OCT_landmarks,
  Fundus_landmarks,
  OCT_minus_Fundus_landmark_1 =
    OCT_landmarks[, 1] - Fundus_landmarks[, 1],
  OCT_minus_Fundus_landmark_2 =
    OCT_landmarks[, 2] - Fundus_landmarks[, 2],
  registration_status = registration_status
)

###################################################################
# 10. Save CSV outputs
###################################################################

output_dir <- "/Users"

write.csv(
  Fused_Landmark_Registered_Data,
  file = file.path(
    output_dir,
    "NRR_Normal_Norm_landmark_registered_fused_curves.csv"
  ),
  row.names = FALSE
)

write.csv(
  Registered_Fundus_Data,
  file = file.path(
    output_dir,
    "NRR_Normal_Norm_landmark_registered_fundus_curves.csv"
  ),
  row.names = FALSE
)

write.csv(
  Warping_Grid_Data,
  file = file.path(
    output_dir,
    "NRR_Normal_Norm_landmark_registration_warping_grids.csv"
  ),
  row.names = FALSE
)

write.csv(
  Landmark_Registration_QC,
  file = file.path(
    output_dir,
    "NRR_Normal_Norm_landmark_registration_qc.csv"
  ),
  row.names = FALSE
)

###################################################################
# 11. Example visual check for one eye
###################################################################

example_eye <- 1

if(registration_status[example_eye] == "ok") {
  
  plot(
    d_NRR,
    NRR_OCT_Normal_Norm_p15_evaluated[example_eye, ],
    type = "l",
    lwd = 2,
    xlab = "Angular Index",
    ylab = "Normalized NRR",
    main = paste(
      "Landmark Registration Example:",
      fullOCTdata_normal_NRR$PatientID_EYE[example_eye]
    )
  )
  
  lines(
    d_NRR,
    NRR_Fundus_Normal_Norm_p15_evaluated[example_eye, ],
    lty = 2,
    lwd = 2
  )
  
  lines(
    d_NRR,
    Fundus_registered[example_eye, ],
    lwd = 2
  )
  
  lines(
    d_NRR,
    Fused_landmark_registered[example_eye, ],
    lwd = 3
  )
  
  points(
    OCT_landmarks[example_eye, ],
    NRR_OCT_Normal_Norm_p15_evaluated[
      example_eye,
      OCT_landmarks[example_eye, ]
    ],
    pch = 19
  )
  
  points(
    Fundus_landmarks[example_eye, ],
    NRR_Fundus_Normal_Norm_p15_evaluated[
      example_eye,
      Fundus_landmarks[example_eye, ]
    ],
    pch = 17
  )
  
  legend(
    "topright",
    legend = c(
      "OCT target curve",
      "Fundus before registration",
      "Fundus after landmark registration",
      "Landmark-registered fused curve",
      "OCT landmarks",
      "Fundus landmarks"
    ),
    lty = c(1, 2, 1, 1, NA, NA),
    lwd = c(2, 2, 2, 3, NA, NA),
    pch = c(NA, NA, NA, NA, 19, 17),
    bty = "n"
  )
}

###################################################################
# 12. Optional: plot the warping function for the same eye
###################################################################

if(registration_status[example_eye] == "ok") {
  
  plot(
    d_NRR,
    Warped_grids[example_eye, ],
    type = "l",
    lwd = 2,
    xlab = "OCT target grid",
    ylab = "Fundus source grid",
    main = paste(
      "Landmark Warping Function:",
      fullOCTdata_normal_NRR$PatientID_EYE[example_eye]
    )
  )
  
  abline(0, 1, lty = 2)
  
  points(
    c(1, OCT_landmarks[example_eye, ], 180),
    c(1, Fundus_landmarks[example_eye, ], 180),
    pch = 19
  )
}

###################################################################
# End of Landmark Registration Pipeline
###################################################################