#' Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Any Numbers of Factors (>=2)
#'
#' @param levels Range of levels,L is F<=L<=choose(F+2,2), where, F is number of factors.
#' @param factors Any number of factors, F >=2.
#' @param w1 Weight of maximum absolute correlation. Between 0 to 1. So that w1+w2+w3=1.
#' @param w2 Weight of Phi_p criterion. Between 0 to 1. So that w1+w2+w3=1.
#' @param w3 Weight of Maxpro criterion. Between 0 to 1. So that w1+w2+w3=1.
#' @param pop_size Default population size is 30. 
#' @param generations Default number of generations is 100.
#' @param mut_prob Mutation probability, by default it is 1/(F-1).
#' @returns Generates Latin hypercube designs for a given number of factor-level combinations.
#' @export
#' @examples
#' \dontrun{
#' library(CompExpDes)
#' wtLHDs(9,3,0.5,0.5,0)
#' }
wtLHDs <- function(levels, factors, w1, w2, w3,
                   pop_size = 30, generations = 100, mut_prob = 1/(factors-1)) {
  
  k = factors
  n = levels
  p = 2 * k + 1
  if (levels > choose(p, 2) || levels < factors || factors < 
      2) {
    return(message("Factors, F (>2) is any integer and Levels, L should be in the range from F to choose(2F+1,2)"))
  }
  t0 <- Sys.time()
  
  # --------------- basic params & checks ----------------
  k <- as.integer(factors)
  n <- as.integer(levels)
  if (k < 2) stop("factors must be >= 2")
  if (n < 1) stop("levels must be >= 1")
  p <- k
  if (sum(c(w1, w2, w3)) != 1 || any(c(w1, w2, w3) < 0) || any(c(w1, w2, w3) > 1)) {
    stop("Weights must be in [0,1] and sum to 1")
  }
  # if (levels > factors^2 || levels < factors) {
  #   stop("Levels L must satisfy F <= L <= choo")
  # }
  
  mut_ops <- c("swap", "inversion", "scramble", "insertion")
  mut_op_probs_init <- c(0.35, 0.25, 0.20, 0.20)
  crossover_ops <- c("uniform", "single_point", "two_point", "matrix_swap", "geometric")
  crossover_op_probs_init <- c(0.3, 0.15, 0.15, 0.2, 0.2)
  adapt_ops <- TRUE; adapt_interval <- 10; adapt_alpha <- 0.15
  tournament_k <- 3; elite_keep <- 2; verbose <- TRUE
  scale_for_metrics <- TRUE
  
  # ----------------- small helpers -----------------------
  safe_sample <- function(x, size = NULL, replace = FALSE, prob = NULL) {
    if (length(x) == 1 && is.numeric(x) && is.null(dim(x))) {
      vec <- seq_len(as.integer(x))
    } else {
      if (is.matrix(x) || is.data.frame(x)) x <- as.vector(x)
      vec <- x
    }
    if (length(vec) == 0 || all(is.na(vec))) stop("safe_sample: input vector is empty or all NA")
    if (is.null(size)) size <- length(vec)
    if (!replace && size > length(vec)) replace <- TRUE
    sample(vec, size = size, replace = replace, prob = prob)
  }
  force_numeric_matrix <- function(x) {
    m <- tryCatch(as.matrix(x), error = function(e) stop("cannot coerce to matrix"))
    if (!is.numeric(m)) {
      m2 <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m)))
      orig_na <- is.na(as.vector(m)); coerced_na <- is.na(as.vector(m2))
      if (any(coerced_na & !orig_na)) stop("matrix contains non-numeric entries that cannot be coerced")
      m <- m2
    }
    storage.mode(m) <- "numeric"
    m
  }
  
  # ----------------- rotation / triangular builder ------------
  rotation_mat <- function(matrix_in) {
    a <- as.matrix(matrix_in)
    ini <- seq_len(ncol(a))
    final <- NULL
    for (i in 1:(ncol(a) - 1)) final <- rbind(final, ini - i)
    final <- rbind(ini, final)
    final <- final %% ncol(a); final[final == 0] <- ncol(a)
    fmat <- NULL
    for (j in 1:ncol(a)) {
      out <- NULL
      for (k2 in 1:ncol(a)) out <- cbind(out, a[, final[j, k2]])
      fmat <- rbind(fmat, out)
    }
    fmat
  }
  
  # require (p-1) even as your partition uses (p-1)/2
  if ((p - 1) %% 2 != 0) stop("This implementation requires (factors - 1) to be even (i.e. factors odd).")
  
  seq1 <- seq_len(p - 1)
  seq2_template <- seq(p - 1, 2)
  tri_list <- list()
  for (j in seq1) {
    a <- j; b <- j; i <- 1
    seq2 <- seq2_template
    while (i <= length(seq2)) {
      b <- c(b, a + seq2[i])
      a <- b[length(b)]; i <- i + 1
    }
    tri_list <- append(tri_list, list(b))
    seq2_template <- seq2_template[-length(seq2_template)]
  }
  Tri1 <- tri_list; Tri2 <- Tri1[length(Tri1):1]
  array1 <- NULL
  for (i in seq_along(Tri1)) array1 <- rbind(array1, c(Tri1[[i]], Tri2[[i]]))
  final_array <- t(array1)
  partition1 <- final_array[, 1:((p - 1)/2), drop = FALSE]
  des <- rotation_mat(partition1)
  
  nr_des <- nrow(des)
  if (nr_des %% p != 0) stop("Unexpected design size from rotation step")
  reg_store <- list()
  for (a_idx in seq_len(nr_des / p)) reg_store <- append(reg_store, list(des[((a_idx * p) - (p - 1)):(a_idx * p), , drop = FALSE]))
  full_design <- do.call(rbind, reg_store)
  LHD2 <- full_design
  
  # ----------------- mark > n as NA (safe) -------------------
  if (nrow(LHD2) > n) {
    deleted_numbers <- seq(n + 1, nrow(LHD2))
    LHD2[LHD2 %in% deleted_numbers] <- NA
  }
  
  # ----------------- form LHD2_revised robustly -----------------
  # collect non-NA entries per column and pad shorter columns with NA
  cols_list <- lapply(seq_len(ncol(LHD2)), function(j) {
    vals <- LHD2[, j]
    vals[!is.na(vals)]
  })
  if (length(cols_list) == 0) {
    LHD2_revised <- matrix(NA_integer_, nrow = 0, ncol = max(1, k))
  } else {
    maxlen <- max(sapply(cols_list, length))
    cols_padded <- lapply(cols_list, function(v) {
      if (length(v) < maxlen) c(v, rep(NA_integer_, maxlen - length(v))) else v
    })
    LHD2_revised <- do.call(cbind, cols_padded)
    storage.mode(LHD2_revised) <- "numeric"
  }
  # pad columns if fewer than k (so indexing up to k is safe)
  if (ncol(LHD2_revised) < k) {
    extra_cols <- k - ncol(LHD2_revised)
    pad <- matrix(NA_integer_, nrow = nrow(LHD2_revised), ncol = extra_cols)
    LHD2_revised <- cbind(LHD2_revised, pad)
  }
  
  # ----------------- build initial_lhd safely -----------------
  initial_lhd <- matrix(NA_integer_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    rawcol <- integer(0)
    if (nrow(LHD2_revised) > 0) rawcol <- as.integer(na.omit(LHD2_revised[, j]))
    used <- integer(0); result_col <- rep(NA_integer_, n)
    pos <- 1
    if (length(rawcol) > 0) {
      for (v in rawcol) {
        if (pos > n) break
        if (!is.na(v) && v >= 1L && v <= n && !(v %in% used)) {
          used <- c(used, v)
          result_col[pos] <- v
          pos <- pos + 1
        }
      }
    }
    remaining_vals <- setdiff(seq_len(n), used)
    na_positions <- which(is.na(result_col))
    if (length(na_positions) > 0) {
      if (length(remaining_vals) >= length(na_positions)) sampled <- sample(remaining_vals, length(na_positions), replace = FALSE)
      else sampled <- sample(remaining_vals, length(na_positions), replace = TRUE)
      result_col[na_positions] <- sampled
    }
    if (!all(sort(result_col) == seq_len(n))) result_col <- sample(seq_len(n), n)
    initial_lhd[, j] <- as.integer(result_col)
  }
  storage.mode(initial_lhd) <- "numeric"
  
  # ----------------- fitness helpers -----------------------
  fitness_cache <- new.env(hash = TRUE)
  is_valid_LHD <- function(mat, n_check, k_check) {
    ok <- tryCatch({
      m <- as.matrix(mat); if (any(is.na(m))) return(FALSE)
      if (nrow(m) != n_check || ncol(m) != k_check) return(FALSE)
      for (j in seq_len(ncol(m))) if (!all(sort(as.integer(m[, j])) == seq_len(n_check))) return(FALSE)
      TRUE
    }, error = function(e) FALSE)
    ok
  }
  repair_LHD <- function(m, n_check) {
    m <- as.matrix(m)
    if (nrow(m) != n_check) {
      m <- do.call(cbind, replicate(ncol(m), sample(seq_len(n_check)), simplify = FALSE))
      return(force_numeric_matrix(m))
    }
    k_local <- ncol(m)
    for (j in seq_len(k_local)) {
      col <- as.integer(round(m[, j])); col[ !(col %in% seq_len(n_check)) ] <- NA_integer_
      seen <- integer(n_check); keep <- logical(length(col))
      for (i in seq_along(col)) {
        v <- col[i]; if (is.na(v)) next
        if (seen[v] == 0L) { keep[i] <- TRUE; seen[v] <- 1L }
      }
      missing <- which(seen == 0L); replace_idx <- which(!keep)
      if (length(replace_idx) != length(missing)) { m[, j] <- sample(seq_len(n_check)); next }
      if (length(missing) > 0L) m[replace_idx, j] <- sample(missing, length(missing)) else m[, j] <- col
    }
    storage.mode(m) <- "numeric"; m
  }
  
  safe_metric_call_inf <- function(fun, mat, fun_name) {
    out <- tryCatch({
      mm <- force_numeric_matrix(mat)
      if (scale_for_metrics && n > 1) mm2 <- apply(mm, 2, function(col) (as.numeric(col) - 1) / (n - 1)) else mm2 <- mm
      res <- fun(mm2)
      if (!is.numeric(res) || length(res) != 1 || !is.finite(res)) stop(paste0(fun_name, " returned non-finite or non-scalar"))
      as.numeric(res)
    }, error = function(e) {
      if (verbose) cat("Metric", fun_name, "failed (treated as Inf):", conditionMessage(e), "\n"); Inf
    })
    out
  }
  
  calculate_fitness <- function(mat, w1, w2, w3) {
    m <- force_numeric_matrix(mat)
    if (any(is.na(m))) stop("calculate_fitness: matrix contains NA")
    if (nrow(m) != n) stop("calculate_fitness: input row count mismatch")
    m1 <- safe_metric_call_inf(MAC, m, "MAC")
    m2 <- safe_metric_call_inf(function(x) PhipMeasure(x, q = 2), m, "PhipMeasure")
    m3 <- safe_metric_call_inf(Maxpro_Measure, m, "Maxpro_Measure")
    if (!is.finite(m1) || !is.finite(m2) || !is.finite(m3)) return(Inf)
    w1 * m1 + w2 * m2 + w3 * m3
  }
  
  get_fitness <- function(mat, w1, w2, w3) {
    key <- tryCatch({ mm <- force_numeric_matrix(mat); paste(as.numeric(as.vector(t(mm))), collapse = ",") }, error = function(e) paste0("badkey_", sample(1e8, 1)))
    if (exists(key, envir = fitness_cache)) return(get(key, envir = fitness_cache))
    val <- tryCatch({ f <- calculate_fitness(mat, w1, w2, w3); assign(key, f, envir = fitness_cache); f }, error = function(e) Inf)
    val
  }
  
  # ----------------- crossover + mutation functions ----------------
  uniform_crossover <- function(p1, p2) { child <- p1; for (j in seq_len(ncol(p1))) if (runif(1) < 0.5) child[, j] <- p2[, j]; storage.mode(child) <- "numeric"; child }
  single_point_crossover <- function(p1, p2) { nc <- ncol(p1); if (nc <= 1) return(uniform_crossover(p1, p2)); cp <- safe_sample(nc - 1, 1); child <- p1; child[, (cp + 1):nc] <- p2[, (cp + 1):nc]; storage.mode(child) <- "numeric"; child }
  two_point_crossover <- function(p1, p2) { nc <- ncol(p1); if (nc <= 2) return(single_point_crossover(p1, p2)); points <- sort(safe_sample(nc - 1, 2)); child <- p1; child[, (points[1] + 1):points[2]] <- p2[, (points[1] + 1):points[2]]; storage.mode(child) <- "numeric"; child }
  matrix_swap_crossover <- function(p1, p2) { nc <- ncol(p1); if (nc <= 1) return(uniform_crossover(p1, p2)); block_len <- sample(1:max(1, floor(nc/2)), 1); start <- sample(seq_len(max(1, nc - block_len + 1)), 1); cols <- start:(start + block_len - 1); child <- p1; child[, cols] <- p2[, cols]; storage.mode(child) <- "numeric"; child }
  geometric_crossover <- function(p1, p2) { uniform_crossover(p1, p2) }
  safe_select_crossover <- function(p1, p2, op_name) {
    switch(op_name, uniform = uniform_crossover(p1, p2), single_point = single_point_crossover(p1, p2),
           two_point = two_point_crossover(p1, p2), matrix_swap = matrix_swap_crossover(p1, p2),
           geometric = geometric_crossover(p1, p2), uniform_crossover(p1, p2))
  }
  
  mutate_swap <- function(vec) { L <- length(vec); if (L <= 1) return(vec); idx <- safe_sample(L, 2); tmp <- vec[idx[1]]; vec[idx[1]] <- vec[idx[2]]; vec[idx[2]] <- tmp; vec }
  mutate_inversion <- function(vec) { L <- length(vec); if (L <= 2) return(vec); idx <- sort(safe_sample(L, 2)); vec[idx[1]:idx[2]] <- rev(vec[idx[1]:idx[2]]); vec }
  mutate_scramble <- function(vec) { L <- length(vec); if (L <= 2) return(vec); idx <- sort(safe_sample(L, 2)); sub <- vec[idx[1]:idx[2]]; vec[idx[1]:idx[2]] <- sample(sub, length(sub)); vec }
  mutate_insertion <- function(vec) { L <- length(vec); if (L <= 2) return(vec); i <- safe_sample(L, 2); from <- i[1]; to <- i[2]; val <- vec[from]; vec2 <- vec[-from]; if (to == 1) vec2 <- c(val, vec2) else if (to > length(vec2)) vec2 <- c(vec2, val) else vec2 <- append(vec2, val, after = to - 1); vec2 }
  apply_mut_op <- function(mat, op_name, mut_prob) {
    m <- as.matrix(mat)
    for (j in seq_len(ncol(m))) if (runif(1) < mut_prob) m[, j] <- switch(op_name, swap = mutate_swap(m[, j]), inversion = mutate_inversion(m[, j]), scramble = mutate_scramble(m[, j]), insertion = mutate_insertion(m[, j]), mutate_swap(m[, j]))
    storage.mode(m) <- "numeric"; m
  }
  
  # ----------------- operator probs ----------------------------
  mut_op_probs <- setNames(mut_op_probs_init / sum(mut_op_probs_init), mut_ops)
  mut_success_counts <- setNames(rep(1e-6, length(mut_ops)), names(mut_op_probs))
  crossover_op_probs <- setNames(crossover_op_probs_init / sum(crossover_op_probs_init), crossover_ops)
  crossover_success_counts <- setNames(rep(1e-6, length(crossover_ops)), crossover_ops)
  if (w1 == 1 && w2 == 0 && w3 == 0) {
    crossover_op_probs <- setNames(1, "uniform")
    mut_op_probs <- setNames(c(0.5, 0.4, 0.1), c("swap", "inversion", "scramble"))
    mut_success_counts <- setNames(rep(1e-6, length(mut_op_probs)), names(mut_op_probs))
  }
  
  # ----------------- initial population (seeded) ----------------
  population <- replicate(pop_size, apply(matrix(1:n, nrow = n, ncol = k), 2, sample), simplify = FALSE)
  if (is_valid_LHD(initial_lhd, n, k)) population[[1]] <- initial_lhd else population[[1]] <- repair_LHD(initial_lhd, n)
  fitness_scores <- sapply(population, function(mat) get_fitness(mat, w1, w2, w3))
  
  best_score <- Inf; best_design <- NULL
  tournament_select_idx <- function(fitness_vec, k_tour) {
    L <- length(fitness_vec); if (L <= 1) return(1L)
    k_use <- min(max(1L, as.integer(k_tour)), L)
    ids <- safe_sample(L, size = k_use, replace = (k_use > L))
    fin <- fitness_vec; fin[!is.finite(fin)] <- .Machine$double.xmax / 4
    ids[which.min(fin[ids])]
  }
  
  history <- numeric(generations)
  if (verbose) pb <- txtProgressBar(min = 0, max = generations, style = 3)
  
  # ----------------- GA main loop -----------------------------
  for (gen in seq_len(generations)) {
    ord <- order(fitness_scores, na.last = TRUE); population <- population[ord]; fitness_scores <- fitness_scores[ord]
    history[gen] <- fitness_scores[1]
    if (is.finite(fitness_scores[1]) && fitness_scores[1] < best_score - 1e-12) { best_score <- fitness_scores[1]; best_design <- population[[1]] }
    
    new_population <- list()
    nkeep <- min(elite_keep, length(population))
    if (nkeep > 0) for (i in seq_len(nkeep)) new_population[[i]] <- population[[i]]
    
    while (length(new_population) < pop_size) {
      p1_idx <- tournament_select_idx(fitness_scores, tournament_k)
      p2_idx <- tournament_select_idx(fitness_scores, tournament_k)
      if (p2_idx == p1_idx && length(population) > 1) p2_idx <- safe_sample(setdiff(seq_along(population), p1_idx), 1)
      parent1 <- population[[p1_idx]]; parent2 <- population[[p2_idx]]
      crossover_choice <- sample(names(crossover_op_probs), 1, prob = as.numeric(crossover_op_probs))
      child <- safe_select_crossover(parent1, parent2, crossover_choice)
      mut_choice <- sample(names(mut_op_probs), 1, prob = as.numeric(mut_op_probs))
      child_mut <- apply_mut_op(child, mut_choice, mut_prob)
      child_mut <- repair_LHD(child_mut, n)
      child_fit <- tryCatch({ if (!is_valid_LHD(child_mut, n, k)) Inf else get_fitness(child_mut, w1, w2, w3) }, error = function(e) Inf)
      sel_vals <- fitness_scores[c(p1_idx, p2_idx)]
      if (all(is.na(sel_vals) | !is.finite(sel_vals))) parent_baseline <- Inf else parent_baseline <- min(sel_vals, na.rm = TRUE)
      if (is.finite(child_fit) && child_fit < parent_baseline - 1e-12) { crossover_success_counts[crossover_choice] <- crossover_success_counts[crossover_choice] + 1; mut_success_counts[mut_choice] <- mut_success_counts[mut_choice] + 1 }
      if (!is_valid_LHD(child_mut, n, k)) child_mut <- apply(matrix(1:n, nrow = n, ncol = k), 2, sample)
      new_population[[length(new_population) + 1]] <- child_mut
    }
    
    population <- new_population
    fitness_scores <- sapply(population, function(mat) get_fitness(mat, w1, w2, w3))
    
    if (adapt_ops && (gen %% adapt_interval == 0)) {
      total_cx_success <- sum(crossover_success_counts)
      if (total_cx_success > 0) {
        cx_success_prop <- crossover_success_counts / total_cx_success
        crossover_op_probs <- (1 - adapt_alpha) * crossover_op_probs + adapt_alpha * cx_success_prop
        crossover_op_probs <- crossover_op_probs / sum(crossover_op_probs)
      }
      crossover_success_counts[] <- 1e-6
      total_mut_success <- sum(mut_success_counts)
      if (total_mut_success > 0) {
        mut_success_prop <- mut_success_counts / total_mut_success
        mut_op_probs <- (1 - adapt_alpha) * mut_op_probs + adapt_alpha * mut_success_prop
        mut_op_probs <- mut_op_probs / sum(mut_op_probs)
      }
      mut_success_counts[] <- 1e-6
    }
    
    if (verbose) setTxtProgressBar(pb, gen)
  }
  
  if (verbose) { close(pb); cat("\n") }
  
  if (is.null(best_design)) { ord <- order(fitness_scores, na.last = TRUE); best_design <- population[[ord[1]]] }
  final_mat <- repair_LHD(best_design, n)
  final_fit <- get_fitness(final_mat, w1, w2, w3)
  final_metrics <- c(MAC = safe_metric_call_inf(MAC, final_mat, "MAC"),
                     Phip = safe_metric_call_inf(function(x) PhipMeasure(x, q = 2), final_mat, "PhipMeasure"),
                     MaxPro = safe_metric_call_inf(Maxpro_Measure, final_mat, "Maxpro_Measure"))
  if (verbose) cat("Time elapsed:", round(difftime(Sys.time(), t0, units = "secs"), 2), "secs\n")
  
  result <- list(design = final_mat,
                 fitness = final_fit,
                 metrics = final_metrics,
                 parameters = list(levels = levels, factors = factors, weights = c(w1 = w1, w2 = w2, w3 = w3), generations_completed = generations),
                 initial_lhd = initial_lhd,
                 history = history)
  class(result) <- "wtLHD"
  return(result$design)
}