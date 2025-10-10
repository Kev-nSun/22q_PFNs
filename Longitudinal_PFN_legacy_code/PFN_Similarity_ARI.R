
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
library(ggplot2)
library(dplyr)
library(fs)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(mclust)

# Root directories
new_root <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/Longitudinal_PFNs/19_PFNs_INV10-_082325/Personalized_FN"      # has e.g. sub-XXXX+ses-baseline/, sub-XXXX+ses-2/
old_root <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/Longitudinal_PFNs/19_PFNs_INV10-_082325/Old_pipeline"         # has e.g. sub-XXXX/

# Sessions (labels are the tail after 'ses-')
baseline_sessions <- c("baselineYear1Arm1")
year2_sessions    <- c("2YearFollowUpYArm1")

# File pattern for PFN mats
new_pfn_glob <- "**/FN.mat"
old_pfn_glob <- "**/final_UV.mat"

get_subject_id <- function(path) {
  nm <- fs::path_file(path)
  m  <- stringr::str_match(nm, "^sub-([A-Za-z0-9]+)(?:\\+.*)?$")  # '+...' optional
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}
get_session <- function(path) {
  nm <- path_file(path)
  m <- str_match(nm, "ses-([A-Za-z0-9]+)$")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

list_sub_sess_dirs <- function(root, recurse = Inf) {
  # list everything under root; fail="never" avoids errors on permissions/placeholders
  all_dirs <- dir_ls(root, type = "directory", recurse = recurse, fail = FALSE)
  # keep dirs whose *basename* contains sub-<ID>+ses-<LABEL>
  keep(all_dirs, ~ grepl("sub-[A-Za-z0-9]+\\+ses-[A-Za-z0-9]+", path_file(.x)))
}
list_sub_dirs <- function(root, recurse = Inf) {
  # list everything under root; fail="never" avoids errors on permissions/placeholders
  all_dirs <- dir_ls(root, type = "directory", recurse = recurse, fail = FALSE)
  # keep dirs whose *basename* contains sub-<ID>
  keep(all_dirs, ~ grepl("sub-[A-Za-z0-9]+$", path_file(.x)))
}

list_new_pfn_files <- function(dir) {
  # Return character vector of matching files under a given sub+sess directory
  if (!dir_exists(dir)) character(0) else
    dir_ls(dir, type = "file", recurse = TRUE, glob = new_pfn_glob)
}
list_old_pfn_files <- function(dir) {
  # Return character vector of matching files under a given sub+sess directory
  if (!dir_exists(dir)) character(0) else
    dir_ls(dir, type = "file", recurse = TRUE, glob = old_pfn_glob)
}

# Needs matrices of shape ~ 59412 x 17; transpose if 17 x 59412
.pick_matrix <- function(obj, expected = c(59412L, 17L)) {
  # convert any 2D numeric array/list entry to a matrix; keep names
  as_mat_2d <- function(x) {
    if (is.matrix(x) && is.numeric(x)) return(x)
    if (is.array(x)  && is.numeric(x) && length(dim(x)) == 2L) return(as.matrix(x))
    NULL
  }
  mats <- lapply(obj, as_mat_2d)
  keep <- vapply(mats, function(x) !is.null(x), logical(1))
  mats <- mats[keep]
  if (length(mats) == 0) {
    stop("No 2D numeric matrices found in .mat object.", call. = FALSE)
  }
  
  nm <- names(mats)
  prior <- c("FN","FINAL_UV","U","UV","V")  # prefer these names if multiple candidates
  
  # helper to choose best index among a set, preferring prior names if present
  pick_best <- function(idx) {
    if (length(idx) == 1L) return(idx)
    if (is.null(nm)) return(idx[1L])
    ord <- order(match(toupper(nm[idx]), prior, nomatch = length(prior) + 1L))
    idx[ord][1L]
  }
  
  dims <- lapply(mats, dim)
  exact_idx <- which(vapply(dims, function(d) identical(d, expected), logical(1)))
  rev_idx   <- which(vapply(dims, function(d) identical(d, rev(expected)), logical(1)))
  
  if (length(exact_idx) >= 1L) {
    i <- pick_best(exact_idx)
    return(mats[[i]])                      # already 59412 x 17
  }
  if (length(rev_idx) >= 1L) {
    i <- pick_best(rev_idx)
    return(t(mats[[i]]))                   # if 17 x 59412 -> transpose to 59412 x 17
  }
  
  # Nothing acceptable; report what we did see
  dim_str <- paste(unique(vapply(dims, function(d) paste(d, collapse = "x"), character(1))), collapse = ", ")
  stop(sprintf("No matrix with dims 59412x17 (or 17x59412) found. Candidates: %s", if (nzchar(dim_str)) dim_str else "<none>"),
       call. = FALSE)
}

# ---- Robust MAT reader: 2 methods ----
read_pfn_mat <- function(path) {
  # 1) Try R.matlab (works for v5/v6)
  out <- tryCatch(R.matlab::readMat(path), error = function(e) e)
  if (!inherits(out, "error")) {
    M <- .pick_matrix(out)
    if (!is.null(M)) {
      d <- dim(M)
      if (identical(d, c(17L, 59412L))) M <- t(M)
      return(M)
    }
  }
  
  # 2) rhdf5 fallback
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    listing <- tryCatch(rhdf5::h5ls(path), error = function(e) NULL)
    if (!is.null(listing)) {
      candidates <- subset(listing, otype == "H5I_DATASET")
      pref <- c("FN","final_UV","U","UV","V")
      idx  <- match(tolower(candidates$name), tolower(pref))
      ord  <- order(ifelse(is.na(idx), length(pref) + 1L, idx))
      for (i in ord) {
        dset <- file.path(candidates$group[i], candidates$name[i])
        M <- tryCatch(rhdf5::h5read(path, dset), error = function(e) NULL)
        if (is.numeric(M) && length(dim(M)) == 2L) {
          d <- dim(M)
          if (identical(d, c(59412L, 17L))) return(as.matrix(M))
          if (identical(d, c(17L, 59412L))) return(t(as.matrix(M)))
        }
      }
      stop("HDF5 read succeeded, but no dataset had dims 59412x17 or 17x59412.", call. = FALSE)
    }
  }
}

softparc_to_hardparc <- function(M, nrow = 59412L, ncol = 17L) {
  if (!is.matrix(M)) M <- as.matrix(M)
  d <- dim(M)
  if (identical(d, c(nrow, ncol))) A <- M
  else if (identical(d, c(ncol, nrow))) A <- t(M)
  else stop(sprintf("Unexpected dims: %sx%s", d[1], d[2]), call. = FALSE)
  
  A[is.na(A)] <- -Inf                        # so NA never wins argmax
  zero_row <- rowSums(is.finite(A) & A != 0) == 0
  PFN_label <- max.col(A, ties.method = "first")
  PFN_label[zero_row] <- NA_integer_
  as.integer(PFN_label)                           # length 59412
}

#1. List subject IDs of those with baseline directories in new PFN pipeline directory
#2. Cycle through those baseline directories to pull those PFN mats- get hard parcellation
#3. Cycle through these subjects to see if they have directories in the old PFN pipeline directory, if so pull those old PFN mats, get hard parcellation
#4. Cycle through these subjects to see if they have year 2 directories in the new PFN pipeline directory, if so pull those year 2 PFN mats, get hard parcellation

# ---- STEP 1: subjects with NEW baseline sessions ----
new_all_sub_sess <- list_sub_sess_dirs(new_root, recurse = Inf)

# Create a full table of all new pipeline sub+sess dirs
new_tbl <- tibble::tibble(
  dir        = new_all_sub_sess,
  subject_id = purrr::map_chr(new_all_sub_sess, get_subject_id),
  session    = purrr::map_chr(new_all_sub_sess, get_session)
) %>%
  dplyr::filter(!is.na(subject_id), !is.na(session))

# Parse subject/session and keep only baseline sessions
new_baseline_tbl <- new_tbl %>% 
  filter(session %in% baseline_sessions)

# Paths to the new-baseline directories (what you'll iterate over to pull PFN mats)
subjects_with_new_baseline <- new_baseline_tbl$dir 

# Unique subject IDs represented in those baseline dirs
new_baseline_subject_ids <- new_baseline_tbl %>% 
  distinct(subject_id) %>% 
  pull(subject_id)


# ---- STEP 2: pull PFN mats from NEW baseline ----
new_baseline_PFNs <- map_dfr(subjects_with_new_baseline, function(sub_sess_dir) {
  files <- list_new_pfn_files(sub_sess_dir)
  tibble(
    subject_id = get_subject_id(sub_sess_dir),
    pipeline   = "new",
    session    = get_session(sub_sess_dir),
    filepath   = files
  )
}) %>%
  filter(!is.na(.data$filepath) & nzchar(.data$filepath)) %>%
  mutate(
    mat = purrr::map(.data$filepath, ~ read_pfn_mat(.x)),
    PFN_labels = purrr::map(mat, ~ softparc_to_hardparc(.x))
  ) %>%
  # keep only successful loads with expected size
  filter(!map_lgl(PFN_labels, is.null)) %>%
  # drop the full matrices to save memory:
  dplyr::select(-mat)


# ---- STEP 3: for those subjects, pull PFN mats from OLD pipeline ----
old_all_sub <- list_sub_dirs(old_root, recurse = Inf)

old_tbl <- tibble(
  dir        = old_all_sub,
  subject_id = map_chr(old_all_sub, get_subject_id),
) %>%
  filter(!is.na(subject_id)) %>%
  semi_join(tibble(subject_id = new_baseline_subject_ids), by = "subject_id")

old_baseline_PFNs <- map_dfr(old_tbl$dir, function(sub_dir) {
  files <- list_old_pfn_files(sub_dir)
  tibble(
    subject_id = get_subject_id(sub_dir),
    pipeline   = "old",
    session    = "baseline",
    filepath   = files
  )
}) %>%
  filter(!is.na(.data$filepath) & nzchar(.data$filepath)) %>%
  mutate(
    mat = purrr::map(.data$filepath, ~ read_pfn_mat(.x)),
    PFN_labels = purrr::map(mat, ~ softparc_to_hardparc(.x))
  ) %>%
  # keep only successful loads with expected size
  filter(!map_lgl(PFN_labels, is.null)) %>%
  # drop the full matrices to save memory:
  dplyr::select(-mat)



# ---- STEP 4: for those subjects, pull NEW year-2 PFN mats ----
new_year2_tbl <- new_tbl %>%
  filter(subject_id %in% new_baseline_subject_ids, session %in% year2_sessions)

new_year2_PFNs <- map_dfr(new_year2_tbl$dir, function(sub_sess_dir) {
  files <- list_new_pfn_files(sub_sess_dir)
  tibble(
    subject_id = get_subject_id(sub_sess_dir),
    pipeline   = "new",
    session    = get_session(sub_sess_dir),
    filepath   = files
  )
}) %>%
  filter(!is.na(.data$filepath) & nzchar(.data$filepath)) %>%
  mutate(
    mat = purrr::map(.data$filepath, ~ read_pfn_mat(.x)),
    PFN_labels = purrr::map(mat, ~ softparc_to_hardparc(.x))
  ) %>%
  # keep only successful loads with expected size
  filter(!map_lgl(PFN_labels, is.null)) %>%
  # drop the full matrices to save memory:
  dplyr::select(-mat)




safe_ARI <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (!any(ok)) return(NA_real_)
  adjustedRandIndex(a[ok], b[ok])  # mclust
}

# ---- STEP 5: Intra-subject ARI (adjusted rand index) ----

# Old vs New (Baseline)
newX <- new_baseline_PFNs %>% mutate(new_idx = row_number()) %>% dplyr::select(subject_id, new_idx, v_new = PFN_labels)
oldX <- old_baseline_PFNs %>% mutate(old_idx = row_number()) %>% dplyr::select(subject_id, old_idx, v_old = PFN_labels)

both_oldnew <- inner_join(newX, oldX, by = "subject_id", relationship = "many-to-many")

pairwise_oldnew <- both_oldnew %>%
  mutate(ARI = map2_dbl(v_new, v_old, ~ safe_ARI(.x, .y)))

ARI_old_new <- pairwise_oldnew %>%
  group_by(subject_id) %>%
  slice_max(ARI, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(subject_id, ARI, comp = "Old vs New (Baseline)")

write.csv(ARI_old_new, file = "./19_PFNs_INV10-_082325/ARI_old_new.csv")

# Baseline vs Year 2 (within new)
y2X <- new_year2_PFNs %>% mutate(y2_idx = row_number()) %>% dplyr::select(subject_id, y2_idx, v_y2 = PFN_labels)

both_bly2 <- inner_join(newX, y2X, by = "subject_id", relationship = "many-to-many")

pairwise_bly2 <- both_bly2 %>%
  mutate(ARI = map2_dbl(v_new, v_y2, ~ safe_ARI(.x, .y)))

ARI_bl_y2 <- pairwise_bly2 %>%
  group_by(subject_id) %>%
  slice_max(ARI, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(subject_id, ARI, comp = "Baseline vs Year 2")

write.csv(ARI_bl_y2, file = "./19_PFNs_INV10-_082325/ARI_bl_y2.csv")


# ---- STEP 6: Plotting intra-subject ARI (adjusted rand index) ----
df_long <- bind_rows(ARI_old_new, ARI_bl_y2) %>% filter(!is.na(ARI))

# x positions
df_long <- df_long %>%
  mutate(x_num = ifelse(comp == "Old vs New (Baseline)", 1, 2))

# lines need wide form
df_wide <- df_long %>%
  dplyr::select(subject_id, comp, ARI) %>%
  pivot_wider(names_from = comp, values_from = ARI) %>%
  filter(!is.na(`Old vs New (Baseline)`) & !is.na(`Baseline vs Year 2`)) %>%
  mutate(x1 = 1, y1 = `Old vs New (Baseline)`,
         x2 = 2, y2 = `Baseline vs Year 2`)

# “zoom” y-range to central 1–99% (tweak probs as you like)
yr <- quantile(df_long$ARI, c(0, 1), na.rm = TRUE)


# ---- With box plots (lines optional)----
p <- ggplot() +
  geom_boxplot(data = df_long, aes(x = x_num, y = ARI, group = x_num),color="black", width=.2, alpha = 0) + #boxplot
  geom_segment(data = df_wide,
               aes(x = x1, xend = x2, y = y1, yend = y2, group = subject_id),
               alpha = 0.3, linewidth = 0.1) +
  geom_point(data = df_long,
             aes(x = x_num, y = ARI),
             #width=.04, 
             size = 1.5, alpha = 0.3) +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Baseline, Old vs pNet pipeline", "Baseline vs Year 2, pNet"),
                     limits = c(0.5, 2.5)) +
  coord_cartesian(ylim =c(0.3,0.7)) + 
  labs(x = NULL, y = "Adjusted Rand Index (ARI)", color = NULL,
       title = "Within-subject PFN Similarity") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16,hjust = 0.5),
    axis.text.x = element_text(size = 10, margin = margin(t = 6)),
    axis.title=element_text(size = 14, margin = margin(t = 6)),
    axis.text.y = element_text(size = 12, margin = margin(t = 6))
  )

print(p)
ggsave('19_PFNs_INV10-_082325/PFN_ARI_19_0.3-0.7_boxplot+lines.png', width = 12, height = 14, dpi = 600, units = "cm")


# ---- STEP 7: Inter-Subject ARI (adjusted rand index) ----
# Pairwise ARI within a group; returns long table of i<j pairs
pairwise_ARI <- function(df_subject) {
  ids <- df_subject$subject_id
  label   <- df_subject$PFN_labels
  n   <- length(ids)
  if (n < 2) {
    return(tibble(subject_i = character(), subject_j = character(), ARI = numeric()))
  }
  
  out_i <- character()
  out_j <- character()
  out_ARI <- numeric()
  
  for (i in seq_len(n - 1L)) {
    label_i <- label[[i]]
    for (j in (i + 1L):n) {
      label_j <- label[[j]]
      ARI <- safe_ARI(label_i,label_j)
      out_i <- c(out_i, ids[i])
      out_j <- c(out_j, ids[j])
      out_ARI <- c(out_ARI, ARI)
    }
  }
  tibble(subject_i = out_i, subject_j = out_j, ARI = out_ARI)
}

# Compute pairwise ARI within each group
ARI_new_baseline <- pairwise_ARI(new_baseline_PFNs) %>% mutate(group = "Baseline, pNet")
ARI_old_baseline <- pairwise_ARI(old_baseline_PFNs) %>% mutate(group = "Baseline, Old")
ARI_new_year2    <- pairwise_ARI(new_year2_PFNs)     %>% mutate(group = "Year 2, pNet")


# ---- STEP 8: Plotting Inter-Subject ARI (adjusted rand index) ----
df_long <- bind_rows(
  ARI_old_baseline %>% mutate(group = "Baseline, Old"),
  ARI_new_baseline %>% mutate(group = "Baseline, pNet"),
  ARI_new_year2    %>% mutate(group = "Year 2, pNet")
)


# x positions- choose the left→right order you want on the x-axis
x_labels <- c("Baseline, Old", "Baseline, pNet", "Year 2, pNet")

df_long <- df_long %>%
  dplyr::mutate(
    group  = factor(group, levels = x_labels),
    x_num  = as.integer(group)   # 1, 2, 3
  )

# make longitudinal pair lines
# subject-pair key function so (i,j) == (j,i)
make_pair_key <- function(a, b) ifelse(a <= b,
                                       paste(a, b, sep = "__"),
                                       paste(b, a, sep = "__"))

df_pairs <- df_long %>%
  mutate(
    pair  = make_pair_key(subject_i, subject_j),
    group = factor(group, levels = x_labels),
  )

df_lines <- df_pairs %>%
  arrange(pair, group) %>%
  group_by(pair) %>%
  mutate(x2 = lead(as.integer(group)),
         y2 = lead(ARI),
         x1 = as.integer(group), #1, 2, 3
         y1 = ARI) %>%
  ungroup() %>%
  filter(!is.na(x2), !is.na(y2))  # only keep pairs present in consecutive groups

# ---- With box plots (lines optional)----
p <- ggplot() +
  geom_boxplot(data = df_long, aes(x = x_num, y = ARI, group = x_num),color="black", width=.25, alpha = 0) + #boxplot
  geom_segment(data = df_lines, #longitudinal lines
               aes(x = x1, xend = x2, y = y1, yend = y2, group = pair),
               alpha = 0.2, linewidth = 0.1) +
  geom_point(data = df_long,
             aes(x = x_num, y = ARI),
             #width = .12, 
             size = 1, alpha = 0.1) +
  scale_x_continuous(breaks = c(1:3),
                     labels = x_labels,
                     limits = c(0.5, 3.5)) +
  coord_cartesian(ylim = c(0.3,0.7)) +
  labs(x = NULL, y = "Adjusted Rand Index (ARI)", color = NULL,
       title = "Between-subject PFN Similarity") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16,hjust = 0.5),
    axis.text.x = element_text(size = 10, margin = margin(t = 6)),
    axis.title=element_text(size = 14, margin = margin(t = 6)),
    axis.text.y = element_text(size = 12, margin = margin(t = 6))
  )

print(p)
ggsave('19_PFNs_INV10-_082325/PFN_ARI_inter_subject_19_0.3-0.7_boxplot+lines.png', width = 14, height = 14, dpi = 600, units = "cm")

