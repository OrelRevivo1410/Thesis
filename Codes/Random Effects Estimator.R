# Import libraries
library(hglm)
library(Matrix)
library(readr)
library(dplyr)

setwd("C:/Thesis/Final Files")

# Random effects estimator

# Creating random effects estimator
# DF created by "creating DB" for analysis
df <- read_csv("DB for analysis/result_table_original.csv")

# Sorting the data frame first by "shoe" and then by "subarea"
df <- df[order(df$shoe, df$sub_area), ]
write.csv(df, "DB for analysis/sorted_df.csv", row.names = FALSE)

# Creating filtered df where contact surface is not 0
filtered_df <- df[df$contact_surface != 0, ]

# Adding an indicator to df that there is a contact surface
df$ind_cont <- ifelse(df$contact_surface > 0, 1, 0)

# Defining variables for calculations
n_areas = 14
n_shoes = 387
n_subareas <- 14

# Creating X matrix
# Creating an initial identity matrix of size (n_shoes * n_areas) x n_areas
identity_matrix <- Diagonal(n_subareas)

# Stacking the identity matrices vertically to form the base X matrix
X <- kronecker(rep(1, n_shoes), identity_matrix)

# Looping through each row of the table
for (row in seq_len(nrow(df))) {
  
  # Extracting shoe, sub_area, and contact_surface
  shoe <- df$shoe[row]
  sub_area <- df$sub_area[row]
  contact_surface <- df$contact_surface[row]
  
  # Calculating the row index in X
  row_index <- (shoe - 1) * 14 + sub_area
  
  # Checking if contact_surface is zero and update X
  if (contact_surface == 0) {
    X[row_index, sub_area] <- 0
  }
}

# Converting X to a regular matrix
X <- as.matrix(X)

# Checking the dimensions of the matrix
dim(X)  

# Creating Z
total_rows <- nrow(df)

# Initializing Z
Z <- matrix(0, nrow = n_areas * n_shoes, ncol = n_shoes)

# Filling the matrix
for (shoe in 1:n_shoes) {
  row_start <- (shoe - 1) * n_areas + 1
  row_end <- shoe * n_areas
  Z[row_start:row_end, shoe] <- 1
}


# Updating Z based on contact_surface
for (row in seq_len(nrow(df))) {
  # Extracting data
  shoe <- df$shoe[row]
  sub_area <- df$sub_area[row]
  contact_surface <- df$contact_surface[row]
  
  # Calculating row index
  row_index <- (shoe - 1) * n_areas + sub_area
  
  # Updating Z if no contact surface
  if (contact_surface == 0) {
    Z[row_index, shoe] <- 0
  }
}

# Outcome variable (RACs_num) for each shoe-area combination
y <- df$RACs_num

# Identify rows and columns where there is a contact surface (non-zero entries in Z)
non_zero_rows <- which(rowSums(X) > 0)

# Subsetting the Z matrix to include only the rows and columns with contact surfaces
X_filtered <- X[non_zero_rows,]
Z_filtered <- Z[rowSums(Z == 1) > 0, ]
y_filtered <- y[df$ind_cont == 1]

# Checking y, X, and Z
if (any(is.na(y) | is.nan(y) | is.infinite(y))) stop("Invalid values in y")
if (any(is.na(X) | is.nan(X) | is.infinite(X))) stop("Invalid values in X")
if (any(is.na(Z) | is.nan(Z) | is.infinite(Z))) stop("Invalid values in Z")

# Defining the offset paramteter
offset_parameter <- filtered_df$contact_surface

# Fitting the model using HGLM package
fit <- hglm(
  y = y_filtered,
  X = X_filtered,
  Z = Z_filtered,
  family = poisson(link = log),
  rand.family = Gamma(link = log),
  offset=log(offset_parameter),
  vcovmat = TRUE,
  method = "EQL1",
)

# Summary of the model
summary(fit)

# Extracting and displaying only the fixed-effect coefficients
random_lambda <- exp(fit$fixef)
random_lambda_adj <- random_lambda/10000

# Extracting standard errors
random_std_errors <- fit$SeFe

# 95% confidence level
Z <- 1.96
CI_lower_Ramdom_adj <- exp(fit$fixef - Z * random_std_errors)/10000
CI_upper_Ramdom_adj <- exp(fit$fixef + Z * random_std_errors)/10000

# Combining results into a data frame
Random_adj_results <- data.frame(
  Subarea = 1:14,
  Lambda = random_lambda_adj,
  Variance = random_std_errors^2,
  CI_Lower = CI_lower_Ramdom_adj,
  CI_Upper = CI_upper_Ramdom_adj
)

# Save the reults
write.csv(Random_adj_results, "Goodness Of Fit/Random_adj_results.csv", row.names = FALSE)
write.csv(random_lambda, "Goodness Of Fit/Random_estimator.csv", row.names = FALSE)