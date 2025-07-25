#' Weighted Criteria-Based Latin Hypercube Designs (LHDs) for Any Numbers of Factors (>=2)
#'
#' @param levels Range of levels,L is F<=L<=choose(F+2,2), where, F is number of factors.
#' @param factors Any number of factors, F >=2.
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
#' wtLHDs_prime(9,3,0.5,0.5,0)
wtLHDs<-function (levels, factors, w1,w2,w3,pop_size=30,generations=100,mut_prob=0.05){
  k = factors
  n = levels
  p = 2 * k + 1
  if (levels > choose(p, 2) || levels < factors || factors < 
      2) {
    return(message("Factors, F (>2) is a prime number and Levels, L should be in the range from F to (F^2)"))
  }
  rotation_mat <- function(matrix) {
    a = as.matrix(matrix)
    ini = c(1:ncol(a))
    final = NULL
    for (i in 1:(ncol(a) - 1)) {
      x = ini - i
      final = rbind(final, x)
    }
    final = rbind(ini, final)
    final = final%%ncol(a)
    final[final == 0] <- ncol(a)
    fmat = NULL
    for (j in 1:ncol(a)) {
      output = NULL
      for (k in 1:ncol(a)) {
        output = cbind(output, (a[, final[j, k]]))
      }
      fmat = rbind(fmat, output)
    }
    return(fmat)
  }
  t0 = Sys.time()
  v = choose(p, 2)
  seq1 = c(1:(p - 1))
  seq2 = c((p - 1):2)
  list <- list()
  j = 1
  for (j in seq1) {
    a = c(j)
    b <- c(j)
    i = 1
    while (i <= length(seq2)) {
      b = c(b, a + seq2[i])
      a <- b[length(b)]
      i = i + 1
    }
    list <- append(list, list(b))
    seq2 = seq2[-length(seq2)]
  }
  Tri1 <- list
  Tri2 <- Tri1[length(Tri1):1]
  array1 <- NULL
  for (i in 1:length(Tri1)) {
    array1 <- rbind(array1, c(Tri1[[i]], Tri2[[i]]))
  }
  final_array = t(array1)
  partition1 <- final_array[, 1:((p - 1)/2)]
  des <- rotation_mat(partition1)
  reg_store <- list()
  for (a in 1:(nrow(des)/p)) {
    reg_store <- append(reg_store, list(des[((a * p) - (p - 
                                                          1)):(a * p), ]))
  }
  full_design<-do.call(rbind,reg_store)
  LHD2<-full_design
  
  ###find >n numbers
  deleted_numbers<-tail(1:nrow(LHD2),nrow(LHD2)-n)
  LHD2[LHD2 %in% c(deleted_numbers)] <- NA
  LHD2_revised<-t(t(apply(LHD2, 2, function(col) col[!is.na(col)])))
  
  ###putting them in a list
  split_matrix_rows <- function(mat) {
    n <- nrow(mat)
    k <- ncol(mat)
    chunk_size=2*k+1
    # Indices to split
    split_sizes <- rep(chunk_size, n %/% chunk_size)
    if (n %% chunk_size != 0) {
      split_sizes <- c(split_sizes, n %% chunk_size)
    }
    
    # Create list of matrices
    split_list <- list()
    start_idx <- 1
    for (size in split_sizes) {
      end_idx <- start_idx + size - 1
      split_list[[length(split_list) + 1]] <- mat[start_idx:end_idx, , drop = FALSE]
      start_idx <- end_idx + 1
    }
    
    return(split_list)
  }
  #######
  s=2*factors+1
  if(n==factors){
    revised_list2<-list(LHD2_revised)
  }
  s_plus<-FALSE
  s_minus<-FALSE
  if(n>s){
    revised_list2<-split_matrix_rows(LHD2_revised) 
    if(nrow(revised_list2[[length(revised_list2)]])==1){
      s_plus<-TRUE
      s_plus_row<-revised_list2[length(revised_list2)]
      revised_list2<-revised_list2[-length(revised_list2)]
    }
  }else if(n>=(factors+1)){
    
    upto_F_part<-as.matrix(LHD2_revised[1:factors,])
    rest_part<-as.matrix(LHD2_revised[(1+factors):n,])
    revised_list2<-list(upto_F_part,rest_part)
    if(nrow(rest_part)==1 || ncol(rest_part)==1){
      s_minus<-TRUE
      s_minus_row<-c(rest_part)
      revised_list2<-list(upto_F_part)
    }
  }
  
  #revised_list2<-list(clean_mat,residual_mat)
  ####^^^^^^^^^^^^^^^^^^^^^ GA
  fitness_function <- function(mat, w1, w2, w3) {
    m1 <- MAC(mat)
    m2 <- PhipMeasure(mat,q=2)
    m3 <- Maxpro_Measure(mat)
    
    return(w1 * m1 + w2 * m2 + w3 * m3)
  }
  generate_individual <- function(revised_list2) {
    lapply(revised_list2, function(submat) {
      apply(submat, 2, sample)  # permutes each column
    })
  }
  
  generate_population <- function(revised_list2, pop_size) {
    replicate(pop_size, generate_individual(revised_list2), simplify = FALSE)
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
    scores <- sapply(population, function(indiv) {
      full_LHD <- do.call(rbind, indiv)
      fitness_function(full_LHD, w1, w2, w3)
    })
    return(scores)
  }
  ###########
  run_GA_optimize_LHD <- function(revised_list2, w1, w2, w3, pop_size, generations, mut_prob) {
    population <- generate_population(revised_list2, pop_size)
    
    best_score <- Inf
    best_design <- NULL
    
    for (gen in 1:generations) {
      # Evaluate fitness
      fitness_scores <- evaluate_population(population, w1, w2, w3)
      ranked <- order(fitness_scores)
      population <- population[ranked]
      
      if (fitness_scores[ranked[1]] < best_score) {
        best_score <- fitness_scores[ranked[1]]
        best_design <- population[[1]]
        cat("Generation", gen, "Best Score:", best_score, "\n")
      }
      
      # Next generation
      new_population <- list()
      
      # Elitism: keep top 2
      new_population[[1]] <- population[[1]]
      new_population[[2]] <- population[[2]]
      
      # Reproduce rest
      while (length(new_population) < pop_size) {
        parents <- sample(1:5, 2)  # sample among top 5
        child <- crossover(population[[parents[1]]], population[[parents[2]]])
        child <- mutate_individual(child, mut_prob)
        new_population[[length(new_population) + 1]] <- child
      }
      
      population <- new_population
    }
    
    return(do.call(rbind, best_design))
    
  }
  
  best_ga_lhd<-run_GA_optimize_LHD(revised_list2 ,w1 = w1,w2 = w2,w3 = w3,pop_size = pop_size,
                                   generations = generations,
                                   mut_prob = mut_prob)
  if(s_plus==TRUE){
    des<-rbind(best_ga_lhd,unlist(s_plus_row))
    rownames(des)<-NULL
    return(des)
    return()
  }else if(s_minus==TRUE){
    des<-rbind(best_ga_lhd,s_minus_row)
    rownames(des)<-NULL
    return(des)
  }else{
    return(best_ga_lhd)
  }
}