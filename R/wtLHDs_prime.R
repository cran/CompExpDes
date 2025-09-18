#' Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Prime Numbers
#'
#' @param levels Range of levels,L is F<=L<=F^2, where, F is number of factors.
#' @param factors Any number of prime factors, F >=3.
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
#' wtLHDs_prime(9,3,1,0,0)
#' }
wtLHDs_prime <- function(levels, factors, w1, w2, w3,
                        pop_size = 30, generations = 100, mut_prob = 1/(factors-1)) {
  if(sum(w1,w2,w3)!=1 || any(c(w1,w2,w3)<0) || any(c(w1,w2,w3)>1)){
    return(message(" sum of w_i should be 1 and 0<= w_i <= 1."))
  }
  t0 <- Sys.time()
  
  # ------------------------
  # Internal fixed settings
  # ------------------------
  mut_ops <- c("swap", "inversion", "scramble", "insertion")
  mut_op_probs_init <- c(0.35, 0.25, 0.20, 0.20)
  crossover_ops <- c("uniform", "single_point", "two_point", "matrix_swap", "geometric")
  crossover_op_probs_init <- c(0.3, 0.15, 0.15, 0.2, 0.2)
  
  adapt_ops <- TRUE; adapt_interval <- 10; adapt_alpha <- 0.15
  tournament_k <- 3; elite_keep <- 2; verbose <- TRUE
  
  init_mode <- "safe"        # simple init that randomly fills missing entries
  generate_mode <- "random"  # generate individuals by random column permutations
  repair_mode <- "conservative"
  
  # ------------------------
  # Lightweight helpers
  # ------------------------
  is.prime <- function(value) {
    if (value <= 1) return(FALSE)
    if (value <= 3) return(TRUE)
    if (value %% 2 == 0 || value %% 3 == 0) return(FALSE)
    i <- 5
    while (i * i <= value) {
      if (value %% i == 0 || value %% (i + 2) == 0) return(FALSE)
      i <- i + 6
    }
    TRUE
  }
  
  circshift <- function(v, shift) {
    L <- length(v); shift <- ((shift %% L) + L) %% L
    if (shift == 0) return(v)
    c(tail(v, shift), head(v, L - shift))
  }
  
  fast_sample_int <- function(n, size = 1, replace = FALSE) sample.int(n, size = size, replace = replace)
  
  # ------------------------
  # Input checks
  # ------------------------
  #if (!requireNamespace("CompExpDes", quietly = TRUE)) warning("CompExpDes not found; ensure MAC/PhipMeasure/Maxpro_Measure available")
  #if (length(levels) != 1 || length(factors) != 1) stop("levels and factors must be scalars")
  #if (any(c(w1, w2, w3) < 0)) stop("Weights must be non-negative")
  #if (abs(w1 + w2 + w3 - 1) !=0) return(message("Weights don't sum to 1."))
  
  k <- as.integer(factors); n <- as.integer(levels)
  if (levels > factors^2 || levels < factors || !is.prime(factors) || factors < 3) {
    return(message("Factors, F (>2) must be prime and Levels, L in [F, F^2]"))
  }
  
  # ------------------------
  # Build base submatrices (minimal version)
  # ------------------------
  s <- factors
  row_wise_add <- s * (1:(s - 1)) + 1
  column_wise_add <- (s^2 + 1) - row_wise_add
  basic <- matrix(1:(s^2), nrow = s, ncol = s, byrow = TRUE)
  
  list_of_first_rows <- lapply(1:(s - 1), function(kidx) {
    base_row <- 1
    for (l in 1:(s - 1)) base_row <- c(base_row, base_row[length(base_row)] + row_wise_add[kidx])
    base_row <- base_row %% (s^2); base_row[base_row == 0] <- s^2; base_row
  })
  
  list_of_arrays <- lapply(1:(s - 1), function(m) {
    arr <- t(as.matrix(list_of_first_rows[[m]]))
    for (nrow_i in 1:(s - 1)) arr <- rbind(arr, arr[nrow(arr), ] + column_wise_add[m])
    arr
  })
  
  list_of_arrays1 <- lapply(list_of_arrays, function(mat) mat %% s^2)
  list_of_arrays2 <- lapply(list_of_arrays1, function(mat) { mat[mat == 0] <- s^2; mat })
  list_of_arrays3 <- append(list(basic), list_of_arrays2)
  
  new <- outer(0:(s - 1), 1:s, "+") %% s
  new[new == 0] <- s
  
  Final_list <- lapply(1:s, function(i) list_of_arrays3[[i]][, new[, i]])
  deleted_values <- tail(1:(factors^2), factors^2 - levels)
  LHD_with_NA <- do.call(rbind, Final_list)
  LHD_with_NA[LHD_with_NA %in% deleted_values] <- NA
  
  # ------------------------
  # Simple initial_lhd: replace NA by sampling missing values (random)
  # ------------------------
  initial_lhd <- matrix(NA_integer_, nrow = n, ncol = k)
  
  for (j in seq_len(k)) {
    rawcol <- LHD_with_NA[, j]
    
    # preserve first valid occurrences (if any)
    used <- integer(0)
    result_col <- rep(NA_integer_, n)
    for (i in seq_len(n)) {
      v <- rawcol[i]
      if (!is.na(v) && v >= 1 && v <= n) {
        v <- as.integer(v)
        if (!(v %in% used)) {
          used <- c(used, v)
          result_col[i] <- v
        }
      }
    }
    
    remaining_vals <- setdiff(seq_len(n), used)
    na_positions <- which(is.na(result_col))
    if (length(na_positions) > 0) {
      sampled <- sample(remaining_vals, length(na_positions))
      result_col[na_positions] <- sampled
    }
    if (!all(sort(result_col) == seq_len(n))) result_col <- sample(seq_len(n), n)
    
    initial_lhd[, j] <- as.integer(result_col)
  }
  storage.mode(initial_lhd) <- "numeric"
  
  # ------------------------
  # Repair functions
  # ------------------------
  conservative_repair <- function(mat) {
    m <- as.matrix(mat)
    for (j in seq_len(ncol(m))) {
      col <- as.integer(round(m[, j]))
      col[!(col %in% seq_len(n))] <- NA_integer_
      non_na_idx <- which(!is.na(col))
      if (length(non_na_idx) > 0) {
        vals <- col[non_na_idx]; keep_mask <- !duplicated(vals)
        keep_rows <- non_na_idx[keep_mask]; dup_rows <- non_na_idx[!keep_mask]
      } else {
        keep_rows <- integer(0); dup_rows <- integer(0)
      }
      kept_vals <- if (length(keep_rows) > 0) col[keep_rows] else integer(0)
      missing_vals <- setdiff(seq_len(n), unique(kept_vals))
      replace_positions <- c(dup_rows, which(is.na(col)))
      if (length(replace_positions) > 0) {
        if (length(missing_vals) < length(replace_positions)) missing_vals <- rep(missing_vals, length.out = length(replace_positions))
        col[replace_positions] <- missing_vals[seq_along(replace_positions)]
      }
      if (!all(sort(col) == seq_len(n))) {
        remaining <- setdiff(seq_len(n), unique(col))
        for (r in seq_len(n)) {
          if (duplicated(col)[r] || is.na(col[r])) {
            col[r] <- remaining[1]; remaining <- remaining[-1]
          }
        }
      }
      m[, j] <- as.integer(col)
    }
    storage.mode(m) <- "numeric"
    m
  }
  
  random_repair <- function(mat) {
    m <- as.matrix(mat)
    for (j in seq_len(ncol(m))) m[, j] <- sample(seq_len(n), n, replace = FALSE)
    storage.mode(m) <- "numeric"
    m
  }
  
  repair_lhd <- if (repair_mode == "conservative") conservative_repair else random_repair
  
  # ------------------------
  # Fitness helpers (compact base-R key; no serialize/digest)
  # ------------------------
  fitness_cache <- new.env(hash = TRUE)
  
  is_valid_LHD <- function(mat, levels, factors) {
    m <- tryCatch(as.matrix(mat), error = function(e) return(FALSE))
    if (any(is.na(m))) return(FALSE)
    if (nrow(m) != levels || ncol(m) != factors) return(FALSE)
    for (j in 1:ncol(m)) if (!all(sort(as.integer(round(m[, j]))) == seq_len(levels))) return(FALSE)
    TRUE
  }
  
  calculate_fitness <- function(mat, w1, w2, w3) {
    mnum <- as.matrix(mat); storage.mode(mnum) <- "numeric"
    m1 <- MAC(mnum); m2 <- PhipMeasure(mnum, q = 2); m3 <- Maxpro_Measure(mnum)
    w1 * m1 + w2 * m2 + w3 * m3
  }
  
  # compact matrix key (numbers flattened). faster than serialize->paste of bytes.
  get_fitness <- function(mat, w1, w2, w3) {
    # key built from the matrix entries row-major (transpose to preserve column order)
    key <- paste(as.integer(as.vector(t(mat))), collapse = "_")
    if (exists(key, envir = fitness_cache, inherits = FALSE)) return(get(key, envir = fitness_cache))
    val <- tryCatch({
      f <- calculate_fitness(mat, w1, w2, w3)
      assign(key, f, envir = fitness_cache)
      f
    }, error = function(e) Inf)
    val
  }
  
  # ------------------------
  # Population generation
  # ------------------------
  generate_individual_random <- function(initial_lhd) {
    m <- matrix(nrow = nrow(initial_lhd), ncol = ncol(initial_lhd))
    for (j in seq_len(ncol(initial_lhd))) {
      col <- as.integer(initial_lhd[, j]); col <- col[!is.na(col)]
      if (length(col) == 0) col <- seq_len(nrow(initial_lhd))
      m[, j] <- sample(col, nrow(initial_lhd), replace = FALSE)
    }
    storage.mode(m) <- "numeric"; m
  }
  
  generate_individual_preserve <- function(initial_lhd) {
    m <- matrix(nrow = nrow(initial_lhd), ncol = ncol(initial_lhd))
    for (j in seq_len(ncol(initial_lhd))) {
      col <- as.integer(initial_lhd[, j])
      sh <- fast_sample_int(nrow(initial_lhd), 1) - 1
      m[, j] <- circshift(col, sh)
    }
    storage.mode(m) <- "numeric"; m
  }
  
  generate_individual <- if (generate_mode == "preserve") generate_individual_preserve else generate_individual_random
  generate_population <- function(initial_lhd, pop_size) replicate(pop_size, generate_individual(initial_lhd), simplify = FALSE)
  evaluate_population <- function(pop, w1, w2, w3) sapply(pop, function(indiv) {
    if (!is_valid_LHD(indiv, n, k)) return(Inf)
    tryCatch(get_fitness(indiv, w1, w2, w3), error = function(e) Inf)
  }, simplify = TRUE)
  
  # ------------------------
  # Mutation operators
  # ------------------------
  mutate_swap <- function(vec) { L <- length(vec); if (L <= 1) return(vec); idx <- fast_sample_int(L, 2); tmp <- vec[idx[1]]; vec[idx[1]] <- vec[idx[2]]; vec[idx[2]] <- tmp; vec }
  mutate_inversion <- function(vec) { L <- length(vec); if (L <= 2) return(vec); i <- fast_sample_int(L, 2); a <- min(i); b <- max(i); vec[a:b] <- rev(vec[a:b]); vec }
  mutate_scramble <- function(vec) { L <- length(vec); if (L <= 2) return(vec); i <- fast_sample_int(L, 2); a <- min(i); b <- max(i); vec[a:b] <- sample(vec[a:b], length(a:b)); vec }
  mutate_insertion <- function(vec) { L <- length(vec); if (L <= 2) return(vec); i <- fast_sample_int(L, 2); from <- i[1]; to <- i[2]; val <- vec[from]; vec <- vec[-from]; if (to == 1) vec <- c(val, vec) else if (to >= length(vec) + 1) vec <- c(vec, val) else vec <- append(vec, val, after = to - 1); vec }
  
  apply_mut_op <- function(mat, op_name, mut_prob) {
    m <- as.matrix(mat)
    for (j in seq_len(ncol(m))) if (runif(1) < mut_prob) m[, j] <- switch(op_name,
                                                                          swap = mutate_swap(m[, j]),
                                                                          inversion = mutate_inversion(m[, j]),
                                                                          scramble = mutate_scramble(m[, j]),
                                                                          insertion = mutate_insertion(m[, j]),
                                                                          mutate_swap(m[, j]))
    repair_lhd(m)
  }
  
  # ------------------------
  # Crossover operators
  # ------------------------
  uniform_crossover <- function(p1, p2) { child <- p1; for (j in 1:ncol(p1)) if (runif(1) < 0.5) child[, j] <- p2[, j]; repair_lhd(child) }
  single_point_crossover <- function(p1, p2) { nc <- ncol(p1); if (nc <= 1) return(uniform_crossover(p1, p2)); cp <- fast_sample_int(nc - 1, 1); child <- p1; child[, (cp + 1):nc] <- p2[, (cp + 1):nc]; repair_lhd(child) }
  two_point_crossover <- function(p1, p2) { nc <- ncol(p1); if (nc <= 2) return(single_point_crossover(p1, p2)); points <- sort(fast_sample_int(nc - 1, 2)); child <- p1; child[, (points[1] + 1):points[2]] <- p2[, (points[1] + 1):points[2]]; repair_lhd(child) }
  matrix_swap_crossover <- function(p1, p2) { nc <- ncol(p1); block_size <- fast_sample_int(max(1, nc), 1); start_col <- fast_sample_int(max(1, nc - block_size + 1), 1); end_col <- start_col + block_size - 1; child <- p1; child[, start_col:end_col] <- p2[, start_col:end_col]; repair_lhd(child) }
  geometric_crossover <- function(p1, p2) uniform_crossover(p1, p2)
  
  safe_select_crossover <- function(p1, p2, op_name) {
    switch(op_name,
           uniform = uniform_crossover(p1, p2),
           single_point = single_point_crossover(p1, p2),
           two_point = two_point_crossover(p1, p2),
           matrix_swap = matrix_swap_crossover(p1, p2),
           geometric = geometric_crossover(p1, p2),
           uniform_crossover(p1, p2))
  }
  
  # ------------------------
  # Operator probabilities init
  # ------------------------
  mut_op_probs <- setNames(mut_op_probs_init / sum(mut_op_probs_init), mut_ops)
  mut_success_counts <- setNames(rep(1e-6, length(mut_ops)), names(mut_op_probs))
  crossover_op_probs <- setNames(crossover_op_probs_init / sum(crossover_op_probs_init), crossover_ops)
  crossover_success_counts <- setNames(rep(1e-6, length(crossover_ops)), crossover_ops)
  
  # ------------------------
  # MAIN GA (with progress bar)
  # ------------------------
  population <- generate_population(initial_lhd, pop_size)
  fitness_scores <- evaluate_population(population, w1, w2, w3)
  
  best_score <- Inf; best_design <- NULL
  history <- numeric(generations)
  
  tournament_select_idx <- function(fitness_vec, k) {
    L <- length(fitness_vec); if (L <= 1) return(1L)
    k_use <- min(max(1L, as.integer(k)), L)
    ids <- sample(seq_len(L), size = k_use, replace = (k_use > L))
    fin <- fitness_vec; fin[!is.finite(fin)] <- .Machine$double.xmax/4; ids[which.min(fin[ids])]
  }
  
  # create progress bar if verbose
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = generations, style = 3)
  }
  
  for (gen in seq_len(generations)) {
    ord <- order(fitness_scores, na.last = TRUE); population <- population[ord]; fitness_scores <- fitness_scores[ord]
    history[gen] <- fitness_scores[1]
    if (is.finite(fitness_scores[1]) && fitness_scores[1] < best_score - 1e-12) {
      best_score <- fitness_scores[1]; best_design <- population[[1]]}
    
    new_population <- list()
    nkeep <- min(elite_keep, length(population)); if (nkeep > 0) for (i in seq_len(nkeep)) new_population[[i]] <- population[[i]]
    
    while (length(new_population) < pop_size) {
      p1_idx <- tournament_select_idx(fitness_scores, tournament_k)
      p2_idx <- tournament_select_idx(fitness_scores, tournament_k)
      if (p2_idx == p1_idx && length(population) > 1) p2_idx <- sample(setdiff(seq_along(population), p1_idx), 1)
      parent1 <- population[[p1_idx]]; parent2 <- population[[p2_idx]]
      
      crossover_choice <- sample(names(crossover_op_probs), 1, prob = as.numeric(crossover_op_probs))
      child <- safe_select_crossover(parent1, parent2, crossover_choice)
      
      mut_choice <- sample(names(mut_op_probs), 1, prob = as.numeric(mut_op_probs))
      child_mut <- apply_mut_op(child, mut_choice, mut_prob)
      
      child_full <- if (!is_valid_LHD(child_mut, n, k)) repair_lhd(child_mut) else child_mut
      
      if (identical(child_full, parent1)) {
        child_fit <- fitness_scores[p1_idx]
      } else if (identical(child_full, parent2)) {
        child_fit <- fitness_scores[p2_idx]
      } else {
        child_fit <- tryCatch({ if (!is_valid_LHD(child_full, n, k)) Inf else get_fitness(child_full, w1, w2, w3) }, error = function(e) Inf)
      }
      
      parent_baseline <- min(fitness_scores[c(p1_idx, p2_idx)])
      if (is.finite(child_fit) && child_fit < parent_baseline - 1e-12) {
        crossover_success_counts[crossover_choice] <- crossover_success_counts[crossover_choice] + 1
        mut_success_counts[mut_choice] <- mut_success_counts[mut_choice] + 1
      }
      new_population[[length(new_population) + 1]] <- child_full
    }
    
    population <- new_population
    fitness_scores <- evaluate_population(population, w1, w2, w3)
    
    # update progress bar (percentage)
    if (verbose) setTxtProgressBar(pb, gen)
    
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
  }
  
  # close progress bar
  if (verbose) {
    close(pb)
    cat("\n") # newline after the progress bar
  }
  
  # compute and print average CPU time per iteration
  total_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  gen_done <- max(1, generations)          # avoid divide-by-zero
  avg_secs <- total_secs / gen_done
  #print(sprintf("average CPU time per iteration is: %.4f seconds", avg_secs))
  
  if (is.null(best_design)) { ord <- order(fitness_scores, na.last = TRUE); best_design <- population[[ord[1]]] }
  
  final_mat <- repair_lhd(best_design)
  final_fit <- get_fitness(final_mat, w1, w2, w3)
  final_metrics <- c(MAC = MAC(final_mat), Phip = PhipMeasure(final_mat, q = 2), MaxPro = Maxpro_Measure(final_mat))
  if (verbose) {
    #cat("Final fitness:", final_fit, "\n")
    #cat("Metrics: MAC =", final_MAC, "Phip =", final_Phip, "MaxPro =", final_MaxPro, "\n")
    cat("Time elapsed:", round(difftime(Sys.time(), t0, units = "secs"), 2), "secs\n")
  }
  result <- list(design = final_mat,
                 fitness = final_fit,
                 metrics = final_metrics,
                 parameters = list(levels = levels, factors = factors, weights = c(w1 = w1, w2 = w2, w3 = w3), generations_completed = generations),
                 history = history)
  class(result) <- "wtLHD"
  return(result$design)
}
