# Import libraries
library(fields)
library(car)
library(xtable)
library(ggplot2)
library(gridExtra)
library(magick)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(Matrix)
library(lme4)
library(flexplot)
library(reshape2)
library(numDeriv)


setwd("C:/Thesis/Final Files")

# Goodness of fit
# Setting the seed for reproducibility
set.seed(312)

# Function for creating expected table
calculate_expected_racs <- function(df, 
                                    n_shoes,
                                    n_subareas, 
                                    skip_shoes = NULL) {
  # Create an empty data frame for the output table
  expected_racs_table <- data.frame(shoe = integer(),
                                    sub_area = integer(),
                                    theta_j = numeric(),
                                    prob = numeric(),
                                    expected_racs = numeric(),
                                    actual_racs = numeric(),
                                    contact_surface = numeric())
  
  # Calculate expected RACs for each shoe and sub-area
  for (shoe in 1:n_shoes) {
    # Skip specified shoes
    if (!is.null(skip_shoes) && shoe %in% skip_shoes) next
    
    # Subset data for the current shoe
    shoe_data <- df[df$shoe == shoe, ]
    
    # Calculate total RACs for this shoe
    total_RACs <- sum(shoe_data$RACs_num)
    
    # Calculate denominator for probability calculation
    denominator <- sum(exp(log(shoe_data$theta_j) + log(shoe_data$contact_surface)))
    
    # For each sub-area, calculate probability and expected RACs
    for (sub_area in 1:n_subareas) {
      # Get current row data
      current_row <- shoe_data[shoe_data$sub_area == sub_area, ]
      
      if (nrow(current_row) > 0) {
        # Extract values
        cur_theta <- current_row$theta_j
        cur_surface <- current_row$contact_surface
        cur_actual_RACs <- current_row$RACs_num
        
        # Calculate multinomial probability according to the formula
        prob <- exp(log(cur_theta) + log(cur_surface)) / denominator
        
        # Calculate expected RACs
        expected_RACs <- total_RACs * prob
        
        # Append the results to the output table
        expected_racs_table <- rbind(expected_racs_table, 
                                     data.frame(shoe = shoe, 
                                                sub_area = sub_area, 
                                                theta_j = cur_theta,
                                                prob = prob,
                                                expected_racs = expected_RACs,
                                                actual_racs = cur_actual_RACs,
                                                contact_surface = cur_surface))
      }
    }
  }
  
  # Check sums of expected vs. observed RACs for each shoe
  check_results <- aggregate(cbind(expected_racs, actual_racs) ~ shoe, data = expected_racs_table, sum)
  
  # Identify mismatches
  mismatch_shoes <- check_results[!all.equal(check_results$expected_racs, check_results$actual_racs, tolerance = 1e-6), ]
  
  # If there are mismatches, print warnings
  if (nrow(mismatch_shoes) > 0) {
    warning("Mismatch detected in the following shoes:")
    print(mismatch_shoes)
  } else {
    message("All sums of expected RACs match the observed values.")
  }
  
  return(expected_racs_table)
}

# Random Estimator
Random <- read.csv('Goodness Of Fit/Random_estimator.csv')

# Add sub_area column to Random
Random$sub_area <- seq_len(nrow(Random))

# Creating expected table for original database
sub_area <- read_csv("DB for analysis/chi_original.csv")
df_original <- read_csv("DB for analysis/result_table_original.csv")
n_shoes_original <- 386
n_subareas_original <- 14

# Rename for clarity
lambda <- Random
names(lambda)[which(names(lambda) == "x")] <- "theta_j"

# Sort the data frame by shoe and then by sub_area
df_original <- df_original[order(df_original$shoe, df_original$sub_area), ]
write.csv(df_original, "DB for analysis/sorted_df_original.csv", row.names = FALSE)

# Join the datasets and adjust scales
df_original <- df_original %>%
  left_join(lambda, by = "sub_area")

df_original$theta_j <- df_original$theta_j/10000
df_original$contact_surface <- df_original$contact_surface*10000

write.csv(df_original, "DB for analysis/df_original_for_goodness_of_fit.csv", row.names = FALSE)

expected_racs_table_original <- calculate_expected_racs(df = df_original,
                                                        n_shoes = n_shoes_original,
                                                        n_subareas = n_subareas_original,
                                                        skip_shoes = 127) 



# Creating expected table for new database
sub_area <- read_csv("DB for analysis/chi_new.csv")
df_new <- read_csv("DB for analysis/result_table_new.csv")
n_shoes_new <- 631
n_subareas_new <- 14

# Rename for clarity
lambda <- Random
names(lambda)[which(names(lambda) == "x")] <- "theta_j"

# Sort the data frame by shoe and then by sub_area
df_new <- df_new[order(df_new$shoe, df_new$sub_area), ]
write.csv(df_new, "sorted_df_new.csv", row.names = FALSE)

# Join the datasets and adjust scales
df_new <- df_new %>%
  left_join(lambda, by = "sub_area")

df_new$theta_j <- df_new$theta_j/10000
df_new$contact_surface <- df_new$contact_surface*10000

write.csv(df_new, "DB for analysis/df_new_for_goodness_of_fit.csv", row.names = FALSE)

expected_racs_table_new <- calculate_expected_racs(df = df_new,
                                                   n_shoes = n_shoes_new,
                                                   n_subareas = n_subareas_new,
                                                   skip_shoes = NULL) 
write.csv(expected_racs_table_new, "DB for analysis/expected_racs_table_new.csv", row.names = FALSE)

# Functions for goodness of fit
# Function to calculate test statistic
expected_racs_table_new <- read.csv("DB for analysis/expected_racs_table_new.csv")
calculate_Q <- function(observed, expected) {
  # Ignore cases where expected is 0
  valid_indices <- expected > 0  
  sum((observed[valid_indices] - expected[valid_indices])^2 / expected[valid_indices])
}

# Function to perform permutation test for a single shoe
permutation_test_shoe <- function(shoe_data, n_permutations = 10000, seed = 123) {
  set.seed(seed)
  
  # Extract data for the current shoe
  observed_racs <- shoe_data$actual_racs
  expected_racs <- shoe_data$expected_racs
  total_racs <- sum(observed_racs)
  probabilities <- shoe_data$prob
  
  # Calculate observed test statistic
  observed_stat <- calculate_Q(observed_racs, expected_racs)
  
  # Store permutation results (including observed at the end)
  perm_stats <- numeric(n_permutations + 1)
  
  # Run permutations
  for (i in 1:n_permutations) {
    # Generate random RAC distribution under multinomial model
    perm_racs <- rmultinom(1, total_racs, probabilities)[,1]
    
    # Calculate test statistic for permuted data
    perm_stats[i] <- calculate_Q(perm_racs, expected_racs)
  }
  
  # Add observed statistic to the permutation statistics
  perm_stats[n_permutations + 1] <- observed_stat
  
  cat("min =", min(perm_stats), "\n")
  # Calculate p-value (proportion of permutation statistics >= observed)
  p_value <- mean(abs(perm_stats) >= abs(observed_stat))
  
  # Return results
  return(list(
    observed_stat = observed_stat,
    perm_stats = perm_stats,
    p_value = p_value,
    significant = p_value < 0.05
  ))
}
permutation_test_shoe(shoe_data = expected_racs_table_new[expected_racs_table_new$shoe==234,], n_permutations = 10000, seed = 123)
# Function to apply permutation test to all shoes or a specific shoe
run_permutation_tests <- function(expected_racs_table, shoes = "all", n_permutations = 10000) {
  
  # Determine which shoes to analyze
  if (length(shoes) == 1 && shoes == "all") {
    shoes_to_analyze <- unique(expected_racs_table$shoe)
  } else {
    shoes_to_analyze <- shoes
  }
  
  # Create results dataframe
  results <- data.frame(
    shoe = integer(),
    observed_stat = numeric(),
    p_value = numeric(),
    significant = logical()
  )
  
  # Create a list to store detailed results for each shoe
  detailed_results <- list()
  
  # Process each shoe
  for (shoe in shoes_to_analyze) {
    cat("Processing shoe", shoe, "...\n")
    
    # Get data for current shoe
    shoe_data <- expected_racs_table[expected_racs_table$shoe == shoe, ]
    
    # Check if there's enough data for this shoe
    if (nrow(shoe_data) < 2) {
      cat("Skipping shoe", shoe, "- insufficient data\n")
      next
    }
    
    # Run permutation test
    test_result <- permutation_test_shoe(shoe_data, n_permutations)
    
    # Add to results dataframe
    results <- rbind(results, data.frame(
      shoe = shoe,
      observed_stat = test_result$observed_stat,
      p_value = test_result$p_value,
      significant = test_result$significant
    ))
    
    # Store detailed results
    detailed_results[[as.character(shoe)]] <- test_result
  }
  
  # Return results
  return(list(
    summary = results,
    details = detailed_results
  ))
}

# Function to visualize permutation test results for a single shoe
plot_permutation_results <- function(test_result, shoe_number) {
  # Create dataframe for plotting
  plot_data <- data.frame(
    statistic = test_result$perm_stats
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = statistic)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = test_result$observed_stat, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = test_result$observed_stat, y = 0, 
             label = sprintf("Observed\np-value = %.3f", test_result$p_value),
             hjust = -0.1, vjust = 0) +
    labs(
      title = paste("Permutation Test Results for Shoe", shoe_number),
      subtitle = ifelse(test_result$significant, "Significant (p < 0.05)", "Not Significant (p >= 0.05)"),
      x = "Q Statistic",
      y = "Frequency"
    ) +
    theme_minimal()
  
  return(p)
}

# Run analysis on all shoes - original data
result_all_original <- run_permutation_tests(expected_racs_table_original, shoes = "all", n_permutations = 10000)
write.csv(result_all_original$summary, "Goodness Of Fit/permutation_test_results_original.csv", row.names = FALSE)

# Additional analysis of results
significant_count_original <- sum(result_all_original$summary$significant)
significant_percent_original <- round(significant_count_original / nrow(result_all_original$summary) * 100, 1)

cat("Number of shoes with significant deviation from model:", significant_count_original, 
    "(", significant_percent_original, "% of all shoes)\n")

# Create histogram of p-values
p_values_hist_original <- ggplot(result_all_original$summary, aes(x = p_value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of p-values across all shoes",
    x = "p-value",
    y = "Frequency"
  ) +
  theme_minimal()

print(p_values_hist_original)
ggsave("Goodness Of Fit/p_values_distribution_original.png", p_values_hist_original, width = 8, height = 6)

# Create a table of shoes with significant deviation
if (significant_count_original > 0) {
  significant_shoes_original <- result_all_original$summary[result_all_original$summary$significant, ]
  significant_shoes_original <- significant_shoes_original[order(significant_shoes_original$p_value), ]
  write.csv(significant_shoes_original, "significant_deviations_original.csv", row.names = FALSE)
  
  cat("Top 5 shoes with most significant deviation from model:\n")
  print(head(significant_shoes_original, 5))
}

# Run analysis on all shoes - new data
result_all_new <- run_permutation_tests(expected_racs_table_new, shoes = "all", n_permutations = 10000)
write.csv(result_all_new$summary, "Goodness Of Fit/permutation_test_results_new.csv", row.names = FALSE)

# Additional analysis of results
significant_count_new <- sum(result_all_new$summary$significant)
significant_percent_new <- round(significant_count_new / nrow(result_all_new$summary) * 100, 1)

cat("Number of shoes with significant deviation from model:", significant_count_new, 
    "(", significant_percent_new, "% of all shoes)\n")

# Create histogram of p-values
p_values_hist_new <- ggplot(result_all_new$summary, aes(x = p_value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of p-values across all shoes",
    x = "p-value",
    y = "Frequency"
  ) +
  theme_minimal()

print(p_values_hist_new)
ggsave("Goodness Of Fit/p_values_distribution_new.png", p_values_hist_new, width = 8, height = 6)

# Create a table of shoes with significant deviation
if (significant_count_new > 0) {
  significant_shoes_new <- result_all_new$summary[result_all_new$summary$significant, ]
  significant_shoes_new <- significant_shoes_new[order(significant_shoes_new$p_value), ]
  write.csv(significant_shoes_new, "Goodness Of Fit/significant_deviations_new.csv", row.names = FALSE)
  
  cat("Top 5 shoes with most significant deviation from model:\n")
  print(head(significant_shoes_new, 5))
}

# Combined P-values
p_vals <- read.csv('Goodness Of Fit/permutation_test_results_new.csv')
# Fisher's Method (Chi-square)
combine_fisher <- function(p_values) {
  # Compute the test statistic: -2 * sum(log(p))
  test_stat <- -2 * sum(log(p_values))
  # Under the null, this follows a chi-square distribution with 2k degrees of freedom
  p_combined <- pchisq(test_stat, df = 2 * length(p_values), lower.tail = FALSE)
  return(p_combined)
}

# Stouffer's Method (Z-score method)
combine_stouffer <- function(p_values) {
  # Convert p-values to Z-scores using the inverse normal CDF
  z_scores <- qnorm(1 - p_values)
  # Combine z-scores (standardize by sqrt(k))
  z_combined <- sum(z_scores) / sqrt(length(z_scores))
  # Convert back to a p-value (upper tail)
  p_combined <- pnorm(z_combined, lower.tail = FALSE)
  return(p_combined)
}

# Logistic Method (Tippett-style logit transform)
combine_logistic <- function(p_values) {
  # Apply the logit transform to each p-value
  logit_vals <- log(p_values / (1 - p_values))
  # Sum of transformed values
  stat <- sum(logit_vals)
  # Mean and variance of the logistic distribution under the null
  mean_logit <- 0
  var_logit <- pi^2 / 3
  # Compute a z-score from the summed statistic
  z <- (stat - mean_logit * length(p_values)) / sqrt(var_logit * length(p_values))
  # Convert z-score to a p-value
  p_combined <- pnorm(z, lower.tail = FALSE)
  return(p_combined)
}

# Combine p-values using Fisher's method
combine_fisher(p_vals$p_value)

# Sanitize p-values: replace 0 with a very small value and 1 with just under 1
clean_p <- pmin(pmax(p_vals$p_value, 1e-15), 1 - 1e-15)
# Run Stouffer on sanitized p-values
combine_stouffer(clean_p)
length(p_vals$p_value)
length(clean_p)

# Find where clean_p is not equal to p_vals$p_value
differences <- which(clean_p != p_vals$p_value)

# View the indices where the difference occurs
print(differences)

# Combine p-values using Logistic method
combine_logistic(p_vals$p_value) 
