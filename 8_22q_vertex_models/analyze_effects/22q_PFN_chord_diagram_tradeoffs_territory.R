library(readr)
library(circlize)

#==========================================================
# PATHS
#==========================================================

#TCA
dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_TCA/Chord_diag"
cov_mod <- "_TCA"

#no TCA
# dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_no_TCA/Chord_diag"
# cov_mod <- "_no_TCA"

filename <- paste0("stat_22q_tradeoff_selected_networks_within90pct_all17",cov_mod,".csv")
tradeoff_file <- file.path(
  dirout,
  filename
)

dat <- read_csv(
  tradeoff_file,
  show_col_types = FALSE
)
#==========================================================
# NETWORK DEFINITIONS
# Order corresponds to PFN1 ... PFN17
#==========================================================

PFN_names <- c(
  "DM-B",       # PFN1
  "SM-foot",    # PFN2
  "FP-A",       # PFN3
  "SM-face",    # PFN4
  "DA-post",    # PFN5
  "VS-peri",    # PFN6
  "Sal",        # PFN7
  "Lang",       # PFN8
  "CingOp",     # PFN9
  "VS-cent",    # PFN10
  "SM-lh",      # PFN11
  "DM-A",       # PFN12
  "SM-rh",      # PFN13
  "DA-ant",     # PFN14
  "ParM",       # PFN15
  "Aud",        # PFN16
  "FP-B"        # PFN17
)

PFN_ids <- paste0(
  "PFN",
  1:17
)

# PFN# -> readable name
pfn_name_lookup <- setNames(
  PFN_names,
  PFN_ids
)

#==========================================================
# 17-NETWORK COLOR PALETTE
#==========================================================

PFN_colors <- rgb(
  r = c(
    220, 173, 244, 73, 65, 170, 255, 226, 235,
    102, 33, 170, 7, 0, 216, 78, 204
  ),
  g = c(
    95, 216, 197, 143, 171, 105, 175, 57, 70,
    5, 113, 12, 69, 109, 144, 49, 109
  ),
  b = c(
    135, 230, 115, 191, 93, 190, 205, 93, 150,
    122, 181, 61, 132, 44, 72, 168, 14
  ),
  maxColorValue = 255
)

network_color_lookup <- setNames(
  PFN_colors,
  PFN_names
)

#==========================================================
# DESIRED ORDER AROUND CIRCLE
#==========================================================

network_order <- c(
  "VS-peri",
  "VS-cent",
  "SM-foot",
  "SM-face",
  "SM-lh",
  "SM-rh",
  "Aud",
  "DA-post",
  "DA-ant",
  "Sal",
  "CingOp",
  "FP-A",
  "FP-B",
  "ParM",
  "Lang",
  "DM-B",
  "DM-A"
)

#==========================================================
# SET 22q NETWORK-LEVEL RESULT ARROWS
#==========================================================

show_22q_arrows <- "TRUE"

func_area_arrow <- setNames(
  rep("", length(PFN_names)),
  PFN_names
)

# # Lose area effect = down arrow:
# func_area_arrow["VS-peri"] <- "\u2193"
# 
# # Gain area effect = up arrow:
# func_area_arrow["CingOp"] <- "\u2191"

# Lose area effect = down triangle:
func_area_arrow["VS-peri"]  <- "\u25BC"

# Gain area effect = up triangle:
func_area_arrow["CingOp"] <- "\u25B2"

#==========================================================
# CLEAN INPUT
#==========================================================

dat$Reference_PFN <- as.character(
  dat$Reference_PFN
)

dat$Selected_Tradeoff_PFN <- as.character(
  dat$Selected_Tradeoff_PFN
)

#----------------------------------------------------------
# Identify reference-mask vertex count column
#----------------------------------------------------------

ref_mask_col <- "N_Ref_Directional_Sig_Vertices"
overlap_col  <- "N_Overlap"

if (!ref_mask_col %in% names(dat)) {
  stop(
    paste0(
      "Could not find column: ",
      ref_mask_col,
      "\nAvailable columns are:\n",
      paste(names(dat), collapse = ", ")
    )
  )
}

if (!overlap_col %in% names(dat)) {
  stop(
    paste0(
      "Could not find column: ",
      overlap_col,
      "\nAvailable columns are:\n",
      paste(names(dat), collapse = ", ")
    )
  )
}

#----------------------------------------------------------
# Remove incomplete rows
#----------------------------------------------------------

dat <- dat[
  !is.na(dat$Reference_PFN) &
    !is.na(dat$Selected_Tradeoff_PFN) &
    !is.na(dat$Selected_Tradeoff_Mean_Beta) &
    !is.na(dat$Mean_Beta_Among_All_Sig_Ref_Vertices) &
    !is.na(dat[[ref_mask_col]]) &
    !is.na(dat[[overlap_col]]),
]
#==========================================================
# BUILD REFERENCE-LEVEL METADATA BEFORE LINK FILTERING
#
# This table determines which networks themselves have a
# sufficiently strong directional reference effect.
#
# IMPORTANT:
# This is intentionally created BEFORE filtering on the
# selected tradeoff target beta. Thus:
#
#   bold/sign = property of the reference network itself
#   ribbon    = property of a qualifying network pair
#==========================================================

reference_meta <- unique(
  dat[, c(
    "Reference_PFN",
    "Reference_Direction",
    "Mean_Beta_Among_All_Sig_Ref_Vertices",
    ref_mask_col
  )]
)

# After collapsing identical rows, each reference PFN should
# have only one reference-level direction/effect specification.
if (any(duplicated(reference_meta$Reference_PFN))) {
  stop(
    paste0(
      "Some reference PFNs have inconsistent reference-level ",
      "metadata across tradeoff rows."
    )
  )
}

#==========================================================
# EFFECT-SIZE AND MASK-SIZE THRESHOLDS
#
# A reference PFN can only generate tradeoff links if:
#   1. |reference mean beta| >= 0.025
#   2. |selected tradeoff mean beta| >= 0.025
#   3. reference mask contains at least 25 vertices
#
# Note:
# A weak network can still appear as the TARGET of a
# stronger reference network only if it passes the selected
# tradeoff beta threshold for that directional relationship.
#==========================================================

# min_ref_beta <- median(
#   abs(reference_meta$Mean_Beta_Among_All_Sig_Ref_Vertices),
#   na.rm = TRUE
# )
min_ref_beta <- 0.025

min_tradeoff_beta <- min_ref_beta
min_ref_mask_vertices <- 25

# Reference-only criteria:
#   1. sufficiently strong reference effect
#   2. sufficiently large directional reference mask
#
reference_meta <- reference_meta[
  abs(reference_meta$Mean_Beta_Among_All_Sig_Ref_Vertices) >=
    min_ref_beta &
    reference_meta[[ref_mask_col]] >=
    min_ref_mask_vertices,
]

#==========================================================
# APPLY TRADEOFF LINK FILTERS
#
# A directional link is retained only if:
#   1. the reference effect is sufficiently strong
#   2. the selected target/tradeoff effect is sufficiently strong
#   3. the directional reference mask is sufficiently large
#
# These are PAIRWISE display criteria and do not determine
# whether the reference network itself receives a bold/sign label.
#==========================================================

dat <- dat[
  abs(dat$Mean_Beta_Among_All_Sig_Ref_Vertices) >=
    min_ref_beta &
    abs(dat$Selected_Tradeoff_Mean_Beta) >=
    min_tradeoff_beta &
    dat[[ref_mask_col]] >=
    min_ref_mask_vertices,
]

# Optional summary of retained directional relationships
cat(
  "\nDirectional tradeoff links retained with:\n",
  "  |reference mean beta| >=",
  min_ref_beta,
  "\n  |tradeoff mean beta| >=",
  min_tradeoff_beta,
  "\n  reference mask vertices >=",
  min_ref_mask_vertices,
  "\n\n"
)

print(
  dat[, c(
    "Reference_PFN",
    "Reference_Direction",
    "Selected_Tradeoff_PFN",
    "Mean_Beta_Among_All_Sig_Ref_Vertices",
    "Selected_Tradeoff_Mean_Beta",
    ref_mask_col,
    overlap_col
  )]
)


#==========================================================
# CONVERT PFN# TO READABLE NETWORK NAMES
#
# If the CSV already contains readable names,
# they are retained.
#==========================================================

dat$Reference_Network <- ifelse(
  dat$Reference_PFN %in% PFN_ids,
  unname(
    pfn_name_lookup[
      dat$Reference_PFN
    ]
  ),
  dat$Reference_PFN
)

dat$Tradeoff_Network <- ifelse(
  dat$Selected_Tradeoff_PFN %in% PFN_ids,
  unname(
    pfn_name_lookup[
      dat$Selected_Tradeoff_PFN
    ]
  ),
  dat$Selected_Tradeoff_PFN
)

# Convert reference-level metadata separately.
reference_meta$Reference_Network <- ifelse(
  reference_meta$Reference_PFN %in% PFN_ids,
  unname(
    pfn_name_lookup[
      reference_meta$Reference_PFN
    ]
  ),
  reference_meta$Reference_PFN
)


#==========================================================
# CHECK NETWORK NAMES
#==========================================================

bad_names <- unique(
  c(
    dat$Reference_Network[
      !dat$Reference_Network %in% network_order
    ],
    dat$Tradeoff_Network[
      !dat$Tradeoff_Network %in% network_order
    ],
    reference_meta$Reference_Network[
      !reference_meta$Reference_Network %in% network_order
    ]
  )
)

bad_names <- bad_names[
  !is.na(bad_names)
]

if (length(bad_names) > 0) {
  
  stop(
    paste(
      "Unrecognized network names:",
      paste(
        bad_names,
        collapse = ", "
      )
    )
  )
}


#==========================================================
# DETERMINE GAINING NETWORK
#
# Reference positive:
#   reference network is gaining
#
# Reference negative:
#   target/tradeoff network is gaining
#==========================================================

dat$Gaining_Network <- ifelse(
  dat$Reference_Direction == "positive",
  dat$Reference_Network,
  dat$Tradeoff_Network
)


#==========================================================
# STRENGTH USED FOR RIBBON WIDTH
#
# Aggregate tradeoff burden:
#
#   Tradeoff_Beta_Sum =
#     mean target/tradeoff beta
#     x number of overlapping vertices
#
# This combines:
#   - target effect magnitude
#   - spatial extent over which that target effect was evaluated
#
# Sign is retained in Tradeoff_Beta_Sum.
# Absolute value is used for ribbon width.
#==========================================================

dat$Tradeoff_Beta_Sum <- (
  dat$Selected_Tradeoff_Mean_Beta *
    dat[[overlap_col]]
)

dat$Tradeoff_Strength <- abs(
  dat$Tradeoff_Beta_Sum
)

dat$Tradeoff_Strength_Type <-
  "abs_mean_beta_times_overlapping_vertices"


#==========================================================
# CREATE UNDIRECTED NETWORK-PAIR ID
#
# A -> B and B -> A should map to the same pair.
#==========================================================

network_index <- setNames(
  seq_along(network_order),
  network_order
)

ref_index <- unname(
  network_index[
    dat$Reference_Network
  ]
)

target_index <- unname(
  network_index[
    dat$Tradeoff_Network
  ]
)

dat$Pair_A <- ifelse(
  ref_index < target_index,
  dat$Reference_Network,
  dat$Tradeoff_Network
)

dat$Pair_B <- ifelse(
  ref_index < target_index,
  dat$Tradeoff_Network,
  dat$Reference_Network
)

dat$Pair_ID <- paste(
  dat$Pair_A,
  dat$Pair_B,
  sep = "__"
)


#==========================================================
# COLLAPSE RECIPROCAL RELATIONSHIPS
#
# One ribbon per unique pair.
#
# Reciprocal:
#   both A -> B and B -> A were selected.
#
# Ribbon width:
#   mean aggregate beta burden across directional analyses.
#
# Aggregate beta burden =
#   abs(selected tradeoff mean beta * overlapping vertices)
#
# Ribbon color:
#   gaining network.
#==========================================================

pair_ids <- unique(
  dat$Pair_ID
)

collapsed_list <- vector(
  "list",
  length(pair_ids)
)

for (i in seq_along(pair_ids)) {
  
  this_pair <- dat[
    dat$Pair_ID == pair_ids[i],
  ]
  
  pair_a <- unname(
    this_pair$Pair_A[1]
  )
  
  pair_b <- unname(
    this_pair$Pair_B[1]
  )
  
  
  #--------------------------------------------------------
  # Is relationship reciprocal?
  #--------------------------------------------------------
  
  has_a_to_b <- any(
    this_pair$Reference_Network == pair_a &
      this_pair$Tradeoff_Network == pair_b
  )
  
  has_b_to_a <- any(
    this_pair$Reference_Network == pair_b &
      this_pair$Tradeoff_Network == pair_a
  )
  
  reciprocal <- (
    has_a_to_b &
      has_b_to_a
  )
  
  
  #--------------------------------------------------------
  # Identify gaining network
  #--------------------------------------------------------
  
  gaining_networks <- unique(
    this_pair$Gaining_Network
  )
  
  gaining_networks <- gaining_networks[
    !is.na(gaining_networks)
  ]
  
  if (length(gaining_networks) == 1) {
    
    gaining_network <- gaining_networks[1]
    
  } else {
    
    # If there is an unexpected discrepancy,
    # use the gaining network from the strongest effect.
    strongest_row <- which.max(
      this_pair$Tradeoff_Strength
    )
    
    gaining_network <-
      this_pair$Gaining_Network[
        strongest_row
      ]
    
    warning(
      paste0(
        "Pair ",
        pair_a,
        " - ",
        pair_b,
        " has inconsistent gaining-network assignments. ",
        "Using ",
        gaining_network,
        " from the strongest directional result."
      )
    )
  }
  
  
  #--------------------------------------------------------
  # Combine strength
  #
  # Mean instead of sum:
  # reciprocal relationships should not automatically
  # become twice as thick.
  #--------------------------------------------------------
  
  combined_strength <- mean(
    this_pair$Tradeoff_Strength,
    na.rm = TRUE
  )
  
  collapsed_list[[i]] <- data.frame(
    from = pair_a,
    to = pair_b,
    value = combined_strength,
    
    Mean_Directional_Tradeoff_Beta_Sum =
      mean(
        this_pair$Tradeoff_Beta_Sum,
        na.rm = TRUE
      ),
    
    Mean_Directional_Tradeoff_Beta_SumAbs =
      mean(
        abs(this_pair$Tradeoff_Beta_Sum),
        na.rm = TRUE
      ),
    
    Max_Directional_Tradeoff_Beta_SumAbs =
      max(
        abs(this_pair$Tradeoff_Beta_Sum),
        na.rm = TRUE
      ),
    
    Mean_Directional_Tradeoff_MeanBetaAbs =
      mean(
        abs(this_pair$Selected_Tradeoff_Mean_Beta),
        na.rm = TRUE
      ),
    
    Mean_Overlap_Mask_Vertices =
      mean(
        this_pair[[overlap_col]],
        na.rm = TRUE
      ),
    
    Gaining_Network = gaining_network,
    Reciprocal = reciprocal,
    N_Directional_Links = nrow(this_pair),
    
    Tradeoff_Strength_Type =
      "abs_mean_beta_times_overlapping_vertices",
    
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

links <- do.call(
  rbind,
  collapsed_list
)

links <- links[
  order(
    links$Reciprocal,  # FALSE first, TRUE last
    links$value        # weaker first, stronger last
  ),
]


#==========================================================
# DOMINANT DIRECTION + DISPLAY LABELS
#
# Bold/sign status is a property of the REFERENCE NETWORK,
# not of whether that network ultimately generates a
# displayed tradeoff ribbon.
#
# Networks passing the reference-only criteria:
#   - show (+) or (-)
#   - displayed in bold
#
# Networks not passing:
#   - show network name only
#   - displayed in regular font
#
# Arrows indicate network-level functional area 22q results
#==========================================================

direction_lookup <- setNames(
  rep(NA_character_, length(PFN_names)),
  PFN_names
)

# Reference-network status comes from reference_meta,
# which was defined BEFORE target/link filtering.
passes_threshold <- setNames(
  PFN_names %in%
    unique(reference_meta$Reference_Network),
  PFN_names
)


# Get dominant direction from the reference-level metadata,
# rather than from the filtered link table.
for (network in PFN_names) {
  
  vals <- unique(
    reference_meta$Reference_Direction[
      reference_meta$Reference_Network == network
    ]
  )
  
  vals <- vals[
    !is.na(vals)
  ]
  
  if (length(vals) > 1) {
    stop(
      paste0(
        "Reference network ",
        network,
        " has multiple directional assignments."
      )
    )
  }
  
  if (length(vals) == 1) {
    direction_lookup[network] <- vals
  }
}


# Direction symbols
# Networks without a retained reference-level effect
# receive an empty string, not NA.
direction_symbol <- setNames(
  rep("", length(PFN_names)),
  PFN_names
)

direction_symbol[
  direction_lookup == "positive" &
    !is.na(direction_lookup)
] <- "+"

direction_symbol[
  direction_lookup == "negative" &
    !is.na(direction_lookup)
] <- "-"


# Build labels
display_labels <- setNames(
  PFN_names,
  PFN_names
)

for (network in PFN_names) {
  
  #------------------------------------------
  # Reference-direction label
  #------------------------------------------
  
  ref_dir <- ""
  
  if (passes_threshold[network]) {
    
    ref_dir <- paste0(
      " (",
      direction_symbol[network],
      ")"
    )
  }
  
  #------------------------------------------
  # Network-level 22q effect arrow
  #------------------------------------------
  
  area_arrow <- ""
  
  if (show_22q_arrows) {
    area_arrow <- func_area_arrow[network]
    if (is.na(area_arrow)) area_arrow <- ""
  }
  
  
  #------------------------------------------
  # Final label
  #------------------------------------------
  
  display_labels[network] <- paste0(
    ifelse(
      area_arrow != "",
      paste0(area_arrow, " "),
      ""
    ),
    network,
    ref_dir
  )
}


# Font:
# 2 = bold
# 1 = regular
label_font <- ifelse(
  passes_threshold,
  2,
  1
)

names(label_font) <- PFN_names
#==========================================================
# ALLOCATE RIBBON POSITIONS WITHIN EACH SECTOR
#
# All sectors will have the SAME x-axis range.
#
# We first calculate the total link weight attached to
# each network. The largest of those totals determines
# the common sector range.
#
# Therefore:
#   - all sectors have equal angular width
#   - ribbon width remains globally proportional to aggregate tradeoff burden
#     computed as abs(mean beta * overlap-mask vertex count)
#==========================================================

incident_weight <- setNames(
  rep(
    0,
    length(network_order)
  ),
  network_order
)

for (i in seq_len(nrow(links))) {
  
  incident_weight[
    links$from[i]
  ] <- incident_weight[
    links$from[i]
  ] + links$value[i]
  
  incident_weight[
    links$to[i]
  ] <- incident_weight[
    links$to[i]
  ] + links$value[i]
}

# standardize sector width across scaling and 22q
common_sector_width <- 40
#   max(
#   incident_weight
# ) * 1.15

#----------------------------------------------------------
# Allocate each ribbon a non-overlapping interval
# within each network sector
#----------------------------------------------------------

links$from_start <- NA_real_
links$from_end <- NA_real_

links$to_start <- NA_real_
links$to_end <- NA_real_

current_position <- setNames(
  rep(
    0,
    length(network_order)
  ),
  network_order
)

for (i in seq_len(nrow(links))) {
  
  a <- links$from[i]
  b <- links$to[i]
  
  w <- links$value[i]
  
  # From end
  links$from_start[i] <-
    current_position[a]
  
  links$from_end[i] <-
    current_position[a] + w
  
  # To end
  links$to_start[i] <-
    current_position[b]
  
  links$to_end[i] <-
    current_position[b] + w
  
  # Advance current positions
  current_position[a] <-
    current_position[a] + w
  
  current_position[b] <-
    current_position[b] + w
}

#==========================================================
# COLORS
#==========================================================

# Base ribbon:
# gaining-network hue with moderate transparency
link_colors <- adjustcolor(
  unname(
    network_color_lookup[
      links$Gaining_Network
    ]
  ),
  alpha.f = 0.8
)

# Reciprocal ribbons receive a dark outline
link_borders <- ifelse(
  links$Reciprocal,
  "grey20",
  NA
)

link_lwd <- ifelse(
  links$Reciprocal,
  2.2,
  1
)

#==========================================================
# SAVE COLLAPSED TABLE
#==========================================================

filename <- paste0("22q_tradeoff_chord_collapsed_links_territory",cov_mod,".csv")
write_csv(
  links,
  file.path(
    dirout,
    filename
  )
)

#==========================================================
# DRAW FIGURE
#==========================================================

filename <- paste0("22q_tradeoff_chord_even_networks_overlap_territory_normed_0.025_beta_sw_40",cov_mod,".pdf")
cairo_pdf(
  filename = file.path(
    dirout,
    filename
  ),
  width = 11,
  height = 11,
  family = "Arial Unicode MS"
)

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 2.5,
  track.margin = c(
    0.01,
    0.01
  ),
  cell.padding = c(
    0,
    0,
    0,
    0
  ),
  canvas.xlim = c(
    -1.30,
    1.30
  ),
  canvas.ylim = c(
    -1.30,
    1.30
  )
)

#==========================================================
# INITIALIZE EQUAL-WIDTH NETWORK SECTORS
#
# Every network gets the exact same x-axis range.
#==========================================================

sector_xlim <- matrix(
  rep(
    c(
      0,
      common_sector_width
    ),
    length(network_order)
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(
  sector_xlim
) <- network_order

circos.initialize(
  factors = factor(
    network_order,
    levels = network_order
  ),
  xlim = sector_xlim
)

#==========================================================
# OUTER NETWORK RING
#==========================================================

circos.trackPlotRegion(
  ylim = c(
    0,
    1
  ),
  track.height = 0.10,
  
  panel.fun = function(x, y) {
    
    sector_name <- get.cell.meta.data(
      "sector.index"
    )
    
    xlim <- get.cell.meta.data(
      "xlim"
    )
    
    # Network-colored outer sector
    circos.rect(
      xleft = xlim[1],
      ybottom = 0,
      xright = xlim[2],
      ytop = 1,
      col = network_color_lookup[
        sector_name
      ],
      border = "white",
      lwd = 1
    )
    
    # Network label + dominant direction
    circos.text(
      x = mean(xlim),
      y = 1.2,
      
      labels = display_labels[
        sector_name
      ],
      
      facing = "clockwise",
      niceFacing = TRUE,
      
      adj = c(
        0,
        0.5
      ),
      
      cex = 1.5,
      
      font = label_font[
        sector_name
      ],
      family = "Arial Unicode MS"
    )
  },
  
  bg.border = NA
)

#==========================================================
# DRAW RIBBONS
#==========================================================

for (i in seq_len(nrow(links))) {
  
  #--------------------------------------------------------
  # Main semi-transparent ribbon
  #
  # Ribbon width = aggregate tradeoff burden
  #                abs(mean beta * overlap-mask vertex count)
  # Reciprocal relationships receive a dark outline
  #--------------------------------------------------------
  
  circos.link(
    sector.index1 = links$from[i],
    
    point1 = c(
      links$from_start[i],
      links$from_end[i]
    ),
    
    sector.index2 = links$to[i],
    
    point2 = c(
      links$to_start[i],
      links$to_end[i]
    ),
    
    col = link_colors[i],
    
    border = link_borders[i],
    
    lwd = link_lwd[i]
  )
}

circos.clear()

dev.off()