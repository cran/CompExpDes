#' Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Prime Numbers
#'
#' @param levels Range of levels,L is F<=L<=F^2, where, F is number of factors.
#' @param factors Any number of prime factors, F >=3.
#' @param w1 Weight of maximum absolute correlation. Between 0 to 1. So that w1+w2+w3=1.
#' @param w2 Weight of Phi_p criterion. Between 0 to 1. So that w1+w2+w3=1.
#' @param w3 Weight of Maxpro criterion. Between 0 to 1. So that w1+w2+w3=1.
#' @param pop_size Default population size is 30. 
#' @param generations Default number of generations is 100.
#' @param mut_prob Mutation probability, by default it is 0.05.
#' @returns Generates Latin hypercube designs for a given number of factor-level combinations.
#' @export
#' @examples
#' library(CompExpDes)
#' wtLHDs_prime(9,3,1,0,0)
wtLHDs_prime = function (levels, factors, w1, w2, w3, pop_size = 30,generations = 100, mut_prob = 0.05) {
  t0 <- Sys.time()
  
  if (sum(w1, w2, w3) != 1) {
    return(message("Please enter weights such that sum of weights = 1."))
  }
  
  is.prime <- function(value) {
    if (value <= 1) return(FALSE)
    if (value <= 3) return(TRUE)
    if (value %% 2 == 0 || value %% 3 == 0) return(FALSE)
    i <- 5
    while (i * i <= value) {
      if (value %% i == 0 || value %% (i + 2) == 0) return(FALSE)
      i <- i + 6
    }
    return(TRUE)
  }
  
  if (levels > factors^2 || levels < factors || !is.prime(factors) || factors < 3) {
    return(message("Factors, F (>2) is a prime number and Levels, L should be in the range from F to (F^2)"))
  }
  
  s <- factors
  row_wise_add <- s * (1:(s - 1)) + 1
  column_wise_add <- (s^2 + 1) - row_wise_add
  basic <- matrix(1:s^2, nrow = s, ncol = s, byrow = TRUE)
  
  list_of_first_rows <- lapply(1:(s - 1), function(k) {
    base_row <- 1
    for (l in 1:(s - 1)) {
      base_row <- c(base_row, base_row[length(base_row)] + row_wise_add[k])
    }
    base_row <- base_row %% (s^2)
    base_row[base_row == 0] <- s^2
    base_row
  })
  
  list_of_arrays <- lapply(1:(s - 1), function(m) {
    array <- t(as.matrix(list_of_first_rows[[m]]))
    for (n in 1:(s - 1)) {
      array <- rbind(array, array[nrow(array), ] + column_wise_add[m])
    }
    array
  })
  
  list_of_arrays1 <- lapply(list_of_arrays, function(mat) mat %% s^2)
  list_of_arrays2 <- lapply(list_of_arrays1, function(mat) { mat[mat == 0] <- s^2; mat })
  list_of_arrays3 <- append(list(basic), list_of_arrays2)
  
  new <- outer(0:(s - 1), 1:s, "+") %% s
  new[new == 0] <- s
  
  Final_list <- lapply(1:s, function(i) list_of_arrays3[[i]][, new[, i]])
  
  deleted_values <- tail(1:factors^2, factors^2 - levels)
  LHD_with_NA <- do.call(rbind, Final_list)
  LHD_with_NA[LHD_with_NA %in% deleted_values] <- NA
  LHD_revised <- t(t(apply(LHD_with_NA, 2, function(col) col[!is.na(col)])))
  
  split_matrix_rows <- function(mat) {
    n <- nrow(mat)
    k <- ncol(mat)
    chunk_size <- k
    split_sizes <- rep(chunk_size, n %/% chunk_size)
    if (n %% chunk_size != 0) split_sizes <- c(split_sizes, n %% chunk_size)
    split_list <- list()
    start_idx <- 1
    for (size in split_sizes) {
      end_idx <- start_idx + size - 1
      split_list[[length(split_list) + 1]] <- mat[start_idx:end_idx, , drop = FALSE]
      start_idx <- end_idx + 1
    }
    split_list
  }
  
  revised_list <- split_matrix_rows(LHD_revised)
  one_row <- FALSE
  if (nrow(revised_list[[length(revised_list)]]) == 1) {
    one_row <- TRUE
    the_last_row <- revised_list[[length(revised_list)]]
    revised_list <- revised_list[-length(revised_list)]
  }
  
  fitness_function <- function(mat, w1, w2, w3) {
    m1 <- MAC(mat)
    m2 <- PhipMeasure(mat, q = 2)
    m3 <- Maxpro_Measure(mat)
    return(w1 * m1 + w2 * m2 + w3 * m3)
  }
  
  generate_individual <- function(revised_list) {
    lapply(revised_list, function(submat) {
      apply(submat, 2, sample)
    })
  }
  
  generate_population <- function(revised_list, pop_size) {
    replicate(pop_size, generate_individual(revised_list), simplify = FALSE)
  }
  
  crossover <- function(parent1, parent2) {
    lapply(seq_along(parent1), function(i) {
      if (runif(1) < 0.5) parent1[[i]] else parent2[[i]]
    })
  }
  
  mutate_individual <- function(individual, mut_prob) {
    lapply(individual, function(submat) {
      for (j in 1:ncol(submat)) {
        if (runif(1) < mut_prob) {
          submat[, j] <- sample(submat[, j])
        }
      }
      submat
    })
  }
  
  evaluate_population <- function(population, w1, w2, w3) {
    sapply(population, function(indiv) {
      full_LHD <- do.call(rbind, indiv)
      fitness_function(full_LHD, w1, w2, w3)
    })
  }
  
  run_GA_optimize_LHD <- function(revised_list, w1, w2, w3, pop_size, generations, mut_prob) {
    population <- generate_population(revised_list, pop_size)
    best_score <- Inf
    best_design <- NULL
    for (gen in 1:generations) {
      fitness_scores <- evaluate_population(population, w1, w2, w3)
      ranked <- order(fitness_scores)
      population <- population[ranked]
      if (fitness_scores[ranked[1]] < best_score) {
        best_score <- fitness_scores[ranked[1]]
        best_design <- population[[1]]
        cat("Generation", gen, "Best Score:", best_score, "\n")
      }
      new_population <- list()
      new_population[[1]] <- population[[1]]
      new_population[[2]] <- population[[2]]
      while (length(new_population) < pop_size) {
        parents <- sample(1:5, 2)
        child <- crossover(population[[parents[1]]], population[[parents[2]]])
        child <- mutate_individual(child, mut_prob)
        new_population[[length(new_population) + 1]] <- child
      }
      population <- new_population
    }
    if (one_row) {
      return(rbind(do.call(rbind, best_design), the_last_row))
    } else {
      return(do.call(rbind, best_design))
    }
  }
  
  best_ga_lhd <- run_GA_optimize_LHD(revised_list,
                                     w1 = w1, w2 = w2, w3 = w3,
                                     pop_size = pop_size,
                                     generations = generations,
                                     mut_prob = mut_prob)
  return(best_ga_lhd)
}

