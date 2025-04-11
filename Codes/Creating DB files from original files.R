# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(magick)

# Setting work directory
setwd("C:/Thesis/Final Files")

# Functions for working on the data
# Function to classify pixels into subareas - RAC's data
classify_pixel <- function(x, y) {
  if (y <= -0.35) {
    if (x < 0) {
      return(1)
    } else {
      return(10)
    }
  } else if (y <= -0.15) {
    if (x < 0) {
      return(2)
    } else {
      return(9)
    }
  } else if (y <= 0) {
    if (x < 0) {
      return(3)
    } else {
      return(8)
    }
  } else if (y <= 0.15) { 
    if (x < 0) {
      if (y >= 12 * (x + 0.02)^2 - 0.07) {
        return(11)
      } else {
        return(4)
      }
    } else {
      if (y >= 16 * (x + 0.02)^2 - 0.07) {
        return(14)
      } else {
        return(7)
      }
    }
  } else {  
    if (x < 0) {
      if (y < -12 * (x + 0.01)^2 + 0.4) {
        return(12)
      } else {
        return(5)
      }
    } else {
      if (y  < -20 * (x + 0.01)^2 + 0.4) {
        return(13)
      } else {
        return(6)
      }
    }
  }
  return(NA) 
}

# Function to process contact data file into row_shoeXcol_shoe matrices
process_contact_data <- function(file_path,row_shoe, col_shoe,num_shoe,byrow=FALSE) {
  # Read the entire file using readLines to get data row by row
  contact_dat <- readLines(file_path)
  
  # Number of shoes (rows in the file)
  num_shoe <- length(contact_dat)
  
  # Create empty list to store matrices
  data_list <- vector("list", num_shoe)
  
  # Process each shoe's data
  for(i in 1:num_shoe) {
    # Split the string into individual characters
    char_vector <- strsplit(contact_dat[i], "")[[1]]
    
    # Convert characters to numeric
    num_vector <- as.numeric(char_vector)
    
    # Reshape
    data_list[[i]] <- matrix(num_vector, 
                             nrow = row_shoe, 
                             ncol = col_shoe, 
                             byrow = byrow)
  }
  
  return(data_list)
}

# Function for creating sub_shaped and dividing the shoes to areas
subshape <- function(j,int) {
  res <- matrix(0,row_shoe,col_shoe)
  sep <- -c(-0.5,-0.35,-0.15,0,0.15,0.5)
  sepPix <- floor(rel_row_shoe*sep+row_shoe/2)
  
  if(j<=5) {
    y <- sepPix[c(j+1,j)]
    x<- c(1,(col_shoe-1)/2)
  }
  else {
    y <- sepPix[c(12-j,11-j)]
    x<- c((col_shoe-1)/2+1,col_shoe)
  }
  
  res[((y[1])+1):(y[2]),(x[1]):(x[2])] <- 1
  
  bound <- function(x,k) {
    ifelse(k==4,return(min(max(floor(-(12*((x-col_shoe/2)/rel_row_shoe+0.02)^2-0.07)*rel_row_shoe+row_shoe/2),1),row_shoe)),
           ifelse(k==5,return(min(max(floor(-(-12*((x-col_shoe/2)/rel_row_shoe+0.01)^2+0.4)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                  ifelse(k==6,return(min(max(floor(-(-20*((x-col_shoe/2)/rel_row_shoe+0.01)^2+0.4)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                         ifelse(k==7,return(min(max(floor(-(16*((x-col_shoe/2)/rel_row_shoe+0.02)^2-0.07)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                                return(1)
                         ))))
  }
  
  if(j==4) {
    if(int==0){
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,4)-1),i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[bound(i,4):row_shoe,i] <- 0
      }
    }
  }
  if(j==5)  {
    if(int ==0){
      for(i in (x[1]):(x[2])) {
        res[bound(i,5):row_shoe,i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,5)-1),i] <- 0
      }
    }
  }
  
  if(j==6){
    if(int ==0){
      for(i in (x[1]):(x[2])) {
        res[bound(i,6):row_shoe,i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,6)-1),i] <- 0
      }
    }
  }
  
  if(j==7) {
    if(int==0){
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,7)-1),i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[bound(i,7):row_shoe,i] <- 0
      }
    }
  }
  
  return(res)
  
}

center_and_flip_contact_data <- function(row_shoe, col_shoe, data, rel_Y_cord, rel_x_cord) {
  num_shoe <- length(data)
  for (i in 1:num_shoe) {
    temp <- matrix(0, row_shoe, col_shoe)
    y_range <- floor(row_shoe / 2 - rel_Y_cord):floor(row_shoe / 2 + rel_Y_cord)
    x_range <- floor(col_shoe / 2 - rel_x_cord):floor(col_shoe / 2 + rel_x_cord)
    
    temp[y_range, x_range] <- data[[i]][y_range, x_range]
    
    data[[i]] <- temp
    data[[i]] <- data[[i]][, ncol(data[[i]]):1]  # Flip horizontally
  }
  return(data)
}

save_matrices_as_pdfs <- function(data, path) {
  for (i in seq_along(data)) {
    img <- image_read(as.raster(data[[i]]))  # Convert matrix to raster and read as an image
    pdf_path <- sprintf("%s/image_%d.pdf", path, i)  # Define PDF filename
    image_write(img, path = pdf_path, format = "pdf")  # Save as PDF
  }
}

subAreaCont_adj <- function(i, j, data) {
  
  data[[i]] * subshapes[[j]]
}

process_RACs_data_txt <- function(file_path) {
  # Read the space-separated file into a data frame
  data <- read.table(file_path, header = FALSE, col.names = c("shoe_num", "accidental_index", "x", "y", "accidental_type"))
  
  return(data)
}

# Create original data bases -JESA 

# RACs data
# Reading the data
data_RAC_original<-read.csv("Original Files/locations_data.CSV",header=TRUE) 

# Relevant numbers
row_shoe_original <-307 #307 is the number of rows in each shoe
col_shoe_original <-395 #395 is the number of columns in each shoe
num_shoe_original <-387  #387 is the number of shoes but 386 is the number of shoes with RACs - shoe 127 has no RACS
rel_col_shoe_original<-150  #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
rel_row_shoe_original<-300  #out of the 395 rows only 300 are relevant (contain non zero pixels in some shoes)
rel_x_cord_original<-0.25 #using coordinates as in the locations_data.CSV file the relevant x coordinates are between -.25 and 0.25
rel_Y_cord_original<-0.5 #the relevant Y coordinates are between -0.5 and 0.5

# Classify the RAC's according to sub_area
for (i in 1:length(data_RAC_original$x)){
  data_RAC_original$sub_area[i] = classify_pixel(data_RAC_original$x[i], data_RAC_original$y[i])
}

write.csv(data_RAC_original, "DB for analysis/locations_data_include_subset_original.CSV",row.names=FALSE)

# contact_data
contact_data_original <-  process_contact_data("Original Files/contacts_data.txt",
                                               row_shoe = row_shoe_original, 
                                               col_shoe = col_shoe_original,
                                               num_shoe = num_shoe_original,
                                               byrow = FALSE)

#flipping shoe 9 -Shoe 9 should be mirrored to correspond to all other shoes
shoe9rev <- contact_data_original[[9]]
contact_data_original[[9]] <- contact_data_original[[9]][,ncol(contact_data_original[[9]]):1]
shoe9 <- contact_data_original[[9]]

# Adapting the result to originals
contact_data_original <- lapply(contact_data_original, t)

dims <-dim(as.matrix(contact_data_original[[1]]))
row_shoe <- dims[1]
col_shoe <- dims[2]
num_shoe <- length(contact_data_original)

rel_row_shoe <- 300
rel_col_shoe <- 150
rel_x_cord <-0.25
rel_Y_cord <- 0.5

for(i in 1:num_shoe) {
  temp <- matrix(0,row_shoe,col_shoe)
  temp[floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),
       floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)] <- contact_data_original[[i]][floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),
                                                                                                                          floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)]
  contact_data_original[[i]]<-temp
  contact_data_original[[i]]<-contact_data_original[[i]][,ncol(contact_data_original[[i]]):1]
}


print(dim(contact_data_original[[1]]))
image(contact_data_original[[1]])
save_matrices_as_pdfs(data = contact_data_original, 
                      path= "original_shoes")


subshapes <- list()
for(i in 1:10) {
  subshapes[[i]] <- subshape(i,0)
}

for(i in 4:7) {
  subshapes[[i+7]] <- subshape(i,1)
}


chi_original<-matrix(0,num_shoe,14)

#For each shoe and each subset, the number of pixels that have a contact surface 
for (i in 1:num_shoe)
{
  for(j in 1:14)
  {
    chi_original[i,j]<-sum(subAreaCont_adj(i,j,contact_data_original))
  }
}
chi_original<-chi_original[-127,] # there are no  RACS in shoe 127 and therefore it is removed
chi_original<-chi_original/10000 # chi is divided by 10000 and so now it is be between zero and 1
write.csv(chi_original, "DB for analysis/chi_original.CSV",row.names=FALSE)


# Create original data bases -FAST
# RACs data
data_RAC_new<-process_RACs_data_txt("Original Files/results_data_naomi_fast.txt")

# Relevant numbers
row_shoe_new <-367 #367 is the number of rows in each shoe
col_shoe_new <-255 #255 is the number of columns in each shoe
num_shoe_new <-631  #631 is the number of shoes
rel_col_shoe_new<-150  #Maybe, needs to be fixed - out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
rel_row_shoe_new<-300  #Maybe, needs to be fixed - out of the 395 rows only 300 are relevant (contain non zero pixels in some shoes)
rel_x_cord_new<-0.25 #using coordinates as in the locations_data.CSV file the relevant x coordinates are between -.25 and 0.25
rel_Y_cord_new<-0.5 #the relevant Y coordinates are between -0.5 and 0.5

# Classify the RAC's according to sub_area
for (i in 1:length(data_RAC_new$x)){
  data_RAC_new$sub_area[i] = classify_pixel(data_RAC_new$x[i], data_RAC_new$y[i])
}

names(data_RAC_new)[names(data_RAC_new) == "shoe_num"] <- "shoe"
write.csv(data_RAC_new, "DB for analysis/locations_data_include_subset_new.CSV",row.names=FALSE)

# Working on contact_data
contact_data_new <-  process_contact_data("Original Files/contacts_data_naomi_fast.txt",
                                          row_shoe = row_shoe_new, 
                                          col_shoe = col_shoe_new,
                                          num_shoe = num_shoe_new,
                                          byrow = TRUE)

dims <-dim(as.matrix(contact_data_new[[1]]))
row_shoe <- dims[1]
col_shoe <- dims[2]
num_shoe <- length(contact_data_new)

rel_row_shoe <- 300
rel_col_shoe <- 150
rel_x_cord <-0.25
rel_Y_cord <- 0.5

for(i in 1:num_shoe) {
  temp <- matrix(0,row_shoe,col_shoe)
  temp[floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),
       floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)] <- contact_data_new[[i]][floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),
                                                                                                                     floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)]
  contact_data_new[[i]]<-temp
  contact_data_new[[i]]<-contact_data_new[[i]][,ncol(contact_data_new[[i]]):1]
}

save_matrices_as_pdfs(data = contact_data_new, 
                      path= "new_shoes")


subshapes <- list()
for(i in 1:10) {
  subshapes[[i]] <- subshape(i,0)
}

for(i in 4:7) {
  subshapes[[i+7]] <- subshape(i,1)
}


chi_new<-matrix(0,num_shoe,14)

#For each shoe and each subset, the number of pixels that have a contact surface 
for (i in 1:num_shoe)
{
  for(j in 1:14)
  {
    chi_new[i,j]<-sum(subAreaCont_adj(i,j,contact_data_new))
  }
}
chi_new<-chi_new/10000 # chi is divided by 10000 and so now it is be between zero and 1
write.csv(chi_new, "DB for analysis/chi_new.CSV",row.names=FALSE)


# Functions for creating DB for estimators
# Function to create the aggregated data table
generate_result_table <- function(sub_area_file, data_file, output_file, 
                                  shoe_range = 1:387, sub_area_range = 1:14, 
                                  remove_shoe = NULL) {
  # Read data
  sub_area <- read_csv(sub_area_file)
  data <- read_csv(data_file)
  
  # Summarize RACs by shoe and subarea
  db <- data %>%
    group_by(across(all_of(c('shoe', 'sub_area')))) %>%
    summarize(RACs_num = n(), .groups = "drop")
  
  # Convert columns to factors
  db[['shoe']] <- factor(db[['shoe']])
  db[['sub_area']] <- factor(db[['sub_area']])
  
  # Adjust shoe numbering to remove the specified shoe(s)
  if (!is.null(remove_shoe)) {
    sub_area <- sub_area %>%
      mutate(
        shoe = case_when(
          row_number() < min(remove_shoe) ~ row_number(),
          row_number() >= min(remove_shoe) ~ row_number() + length(remove_shoe)
        )
      )
  } else {
    sub_area <- sub_area %>% mutate(shoe = row_number())
  }
  
  # Create a table with all combinations of shoes and subareas
  result_table <- expand.grid(
    shoe = factor(shoe_range), 
    sub_area = factor(sub_area_range)
  )
  
  # Merge summarized data with the complete table of combinations
  result_table <- result_table %>%
    left_join(db, by = c("shoe", "sub_area")) %>%
    replace_na(list(RACs_num = 0))  # Replace missing RACs_num with 0
  
  # Reshape sub_area to long format
  sub_area_long <- sub_area %>%
    pivot_longer(
      cols = -shoe,
      names_to = "sub_area",
      values_to = "contact_surface"
    ) %>%
    mutate(
      sub_area = as.numeric(gsub("V", "", sub_area)) # Extract numeric sub-area ID
    ) %>%
    arrange(shoe, sub_area)
  
  # Write sub_area_long to CSV
  write.csv(sub_area_long, file = "sub_area_long.csv", row.names = FALSE)
  
  # Convert columns to factors
  sub_area_long$sub_area <- as.factor(sub_area_long$sub_area)
  sub_area_long$shoe <- as.factor(sub_area_long$shoe)
  
  # Merge contact surface data with result_table
  result_table <- result_table %>%
    left_join(sub_area_long, by = c("shoe", "sub_area")) %>%
    replace_na(list(contact_surface = 0))
  
  # Write final result table to CSV
  write.csv(result_table, file = output_file, row.names = FALSE)
  
  return(result_table)
}

# Create original DB
result_table_original <- generate_result_table(sub_area_file = 'DB for analysis/chi_original.csv', 
                                               data_file = "DB for analysis/locations_data_include_subset_original.csv", 
                                               output_file = "DB for analysis/result_table_original.csv",
                                               shoe_range = 1:387, 
                                               sub_area_range = 1:14, 
                                               remove_shoe = 127)

# Create new DB
result_table_new <- generate_result_table(sub_area_file = 'DB for analysis/chi_new.csv',
                                          data_file = "DB for analysis/locations_data_include_subset_new.csv", 
                                          output_file = "DB for analysis/result_table_new.csv",
                                          shoe_range = 1:631, 
                                          sub_area_range = 1:14, 
                                          remove_shoe = NULL)

