SLHDs<- function(slices, factors, levels) {
  t=slices
  k=factors
  m=levels
  t0<-Sys.time()
  ###MAC function
  MAC<-function(matrix){
    matrix=as.matrix(matrix)
    mat<-cor(matrix)
    upper_tri_values<-abs(mat[upper.tri(mat)])
    return(max(upper_tri_values))
  }
  ##
  ####Random LHD generation
  rLHD<-function(m,k){
    rlhd<-NULL
    repeat{
      rlhd<-cbind(rlhd,(sample(1:m)))
      if(ncol(rlhd)==k){
        return(rlhd)
      }
    }
  }
  # Step 1: Generate initial LHDs
  lhds <- lapply(1:t, function(i) rLHD(m, k))
  # Step 2: Optimize full LHD (combine and permute until good MAC value is achieved)
  repeat {
    # Apply sampling to each slice and combine into full design
    dl <- lapply(lhds, function(mat) apply(mat, 2, sample))
    d <- do.call(rbind, dl)
    
    # Break loop if MAC condition is met
    if (MAC(d) < 0.1) break
  }
  
  # Step 3: Create level mapping (efficient replacement of levels)
  set <- sample(1:(t * m))  # Create the mapping set
  sliced_indices <- split(set, rep(1:t, each = m))  # Split set into t slices
  
  # Replace levels in each slice
  lhds_new <- Map(function(slice, levels) {
    apply(slice, 2, function(col) levels[match(col, 1:m)])
  }, lhds, sliced_indices)
  
  # Step 4: Combine updated slices into final SLHD
  slicedLHD <- do.call(rbind, lhds_new)
  slicedLHD1 <- cbind(rep(1:t, each = m), slicedLHD)  # Add slice indices
  # Rename only the first column, remove names from others
  colnames(slicedLHD1) <- c("Slices", paste0("V", 1:(ncol(slicedLHD1)-1)))
  # Step 5: Calculate measures and return results
  t1<-Sys.time()
  lm <- list(
    "SLHD" = slicedLHD1,"Factors"=ncol(slicedLHD1)-1,"Levels"=nrow(slicedLHD1),"Required time"=t1-t0
  )
  
  return(lm)
}