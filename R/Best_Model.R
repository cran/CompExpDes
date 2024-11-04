# Define the function
Best_Model <- function(model, data) {
  polynomial=model
  data=as.data.frame(data)
  # Initialize an empty list to store the best models for each subset order
  final_models <- list()
  ##########
  # if(nrow(data)<length(polynomial)){
  #   len=nrow(data)-1
  # }else{
  #   len=length(polynomial)
  # }
  len=length(polynomial)
  ######
  # Loop over each subset size from 1 to the full length of the polynomial
  for (subset_order in 1:len) {
    # Generate all combinations of the polynomial terms based on current subset size
    combns <- t(combn(length(polynomial), subset_order))
    
    # Initialize variables to store the best model for this subset size
    best_r_squared_sum <- -Inf
    best_model_summary <- NULL
    
    # Iterate over all combinations of terms for the current subset size
    for (i in 1:nrow(combns)) {
      
      # Build the formula dynamically using selected combinations
      formula_str <- paste("data[,ncol(data)] ~", paste(unlist(polynomial[combns[i,]]), collapse = "+"))
      
      # Fit the model
      mdl <- lm(as.formula(formula_str), data = data)
      
      # Extract p-values
      p_values <- summary(mdl)$coefficients[, 4]
      
      # Get the model summary
      model_summary <- summary(mdl)
      
      # Access R-squared and Adjusted R-squared
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      
      # Check if all p-values for the predictors (excluding the intercept) are less than 0.05
      if (all(p_values < 0.05)) {  # Exclude the intercept's p-value
        
        # Calculate the sum of R-squared and Adjusted R-squared
        r_squared_sum <- r_squared + adj_r_squared
        
        # Store the model if it has the highest sum of R-squared values for this subset size
        if (r_squared_sum > best_r_squared_sum) {
          best_r_squared_sum <- r_squared_sum
          best_model_summary <- model_summary
        }
      }
    }
    
    # If a valid model was found for this subset size, add it to the final list
    if (!is.null(best_model_summary)) {
      final_models[[subset_order]] <- best_model_summary
    }
  }
  
  # Check all stored models and select the best one based on the lowest p-value of estimates and highest R-squared + Adjusted R-squared sum
  if (length(final_models) > 0) { 
    return(final_models)
    # adj_r_sq<-unlist(lapply(final_models,function(x)return(x$adj.r.squared)))*100
    # adj_r_sq_adj<-round(adj_r_sq)
    # stop_check<-NULL
    # for(i in 1:(length(adj_r_sq_adj)-1)){
    #   stop_check<-c(stop_check,abs(adj_r_sq_adj[i+1]-adj_r_sq_adj[i]))
    # }
    # ######
    # stop_check_reverse<-rev(stop_check)
    # final_models_reverse<-rev(final_models)
    # for(j in 1:length(stop_check_reverse)){
    #   if(stop_check_reverse[j]>=4){
    #     return(final_models_reverse[[j]])
    #     #return(final_models)
    #   }
    # }
    
  } else {
    return("No valid model found where all p-values are less than 0.05")
  }
}

