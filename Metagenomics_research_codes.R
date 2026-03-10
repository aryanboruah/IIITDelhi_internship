
#code to calculate correlation values and then create LOOCV(Leave one out cross validation) matrix 

markers <- c("Age", "ALP", "ALT", "AST", "BMI", "Cholesterol", "Creatinine", "Diastolic", "Glucose", "HBA1C", "HDL",
             "HSCRP", "Insulin", "LDL", "Systolic", "Triglyceride")
names_studies <- c("reduced_list_Age","reduced_list_ALP" ,"reduced_list_ALT","reduced_list_AST","reduced_list_BMI",         
"reduced_list_cholesterol","reduced_list_Creatinine","reduced_list_Diastolic","reduced_list_Glucose",     
"reduced_list_HBA1C","reduced_list_HDL","reduced_list_HSCRP","reduced_list_Insulin",     
"reduced_list_LDL","reduced_list_Systolic","reduced_list_triglyceride")

load("C:/Users/aryan/Desktop\rf_tg_env.RData") # set the path accordingly(differ individually)
library(randomForest)
library(caret)
library(dplyr)
rfmodel10 <- list()
num_iterations <- 10


for (study_names in names_studies){
  current_study <- get(study_names)
  print(study_names)
  for(marker in markers){
    print(marker)
  
    matrix_name <- paste("CorrMat_", marker, sep = "")
    x <- matrix(ncol = length(current_study), nrow = num_iterations)
    colnames(x) <- names(current_study)
    assign(matrix_name, x)
    
  for (study_num in 1:length(current_study)){ 
    print(paste("List: ", study_names, sep = ""))
    
    for (i in 1:num_iterations) {
      print(paste("Iteration:", i))
      # Extract data for all other studies except the current study
      leave_out_data <- do.call(rbind, current_study[-study_num])
      # Train Random Forest model
      formula <- as.formula(paste(current_study[[study_num]][i], "~ ."))
      rf_model <- randomForest(formula , data = leave_out_data, ntree = 10, mtry = 10)
      print(paste("rfmodel:", i , "complete"))
      # Store the Random Forest model in the list
      rfmodel5[[names(current_study[study_num])]][[paste("rf_model", i)]] <- rf_model
      print(paste("rfmodel:", i , "for", marker, names(current_study[study_num]) , "stored in the list"))
      
      predictions <- predict(rfmodel5[[names(current_study[study_num])]][[i]], newdata=current_study[[study_num]])
      
      corr_value <- cor(predictions, current_study[[study_num]][[i]])
      print(corr_value)
      
      z <- get(matrix_name)
      z[i, names(current_study[study_num])] <- corr_value
      print(paste("correlation value:", i , "stored in the matrix"))
      
    }
    Dynamic_file_path <- paste("C:/Users/aryan/Desktop/boot_df/list_", marker , ".RData", sep="")
    save.image(file = Dynamic_file_path)
    print(paste("Dynamic file saved for Study",marker))
    
  }
    
}
}  



# code to create out of the box correlation matrix 

CorrMatrix <- matrix(ncol = 14, nrow = 14)
colnames(CorrMatrix) <- names(Triglyceride_datasets)
rownames(CorrMatrix) <- names(Triglyceride_datasets)
rfmodel_same <- list()
rfmodel_diff <- list()
load(...............)


for(i in names(Triglyceride_datasets)){
  for (j in names(Triglyceride_datasets)){
    if(i==j){
      current_data <- Triglyceride_datasets[[i]]
      rf_model <- randomForest(Triglyceride ~ ., data = current_data, ntree = 20, mtry = 10)
      rfmodel_same[["Barton_2018"]] <- rf_model
      print("rfmodel stored")
      
      corr_value <- cor(rfmodel_same[[i]]$y, rfmodel_same[[i]]$predicted)
      print(paste("CorrValue:", i,"&", j,  corr_value, sep = ""))
      CorrMatrix[i,j] <- corr_value
      print("CorrValue added in the matrix")
    }else if(i !=j){
      
      
      train_data <- Triglyceride_datasets[[i]]
      test_data <- Triglyceride_datasets[[j]]
      
      rfmodel <- randomForest(Triglyceride ~ ., data = train_data, ntree = 50, mtry = 10)
      rfmodel_diff[[paste(i, "&", j, sep = "_")]] <- rfmodel
      print("rfmodel")
      predictions <- predict(rfmodel_diff[[paste(i, "&", j, sep = "_")]], newdata = test_data)
      corr_value <- cor(predictions, test_data$Triglyceride)
      print(paste("CorrValue:", i, "&", j, corr_value))
      CorrMatrix[i,j] <- corr_value
      print("CorrValue added in the matrix")
    }
    
    
  }
  if(i=="keohaneD_2020" & j == "keohaneD_2020"){
    file_path <- "........."
    save.image(file = file_path)
    print("file saved")
    
  }
  
}

#code to subset common columns from multiple data frames in a list and create another list with those subset columns

columns_to_subset <- c("Triglyceride",rownames(filtered_BatchSpecies_REM_Triglyceride[filtered_BatchSpecies_REM_Triglyceride$dir %in% c(2, 3, -2, -3), ]))
Triglyceride_datasets_boot <- Triglyceride_datasets
for(i in 1:14)
{   
`Triglyceride datasets`[[i]] <- `Triglyceride datasets`[[i]][, c("Age", "BMI", columns_to_subset)]
}


#code to remove n number columns from n number of data frames in a list 
remove_columns <- function(df) {
  df %>%
    select(-Age, -BMI)
}
Triglyceride_datasets <- lapply(`Triglyceride datasets`, remove_columns)

#Replace all the NAs with 0 from all the data frame in a list 

for (i in 1:length(`Triglyceride datasets`)) {
  for (k in order_list_datasets) {
    for (j in 1:length(datasets[[i]])) {
      na_indices <- is.na(datasets[[i]][[j]][[k]])
      
      if (any(na_indices)) {
        print(paste("NA values detected in", k, "for", j, "in", i))
        print(which(na_indices))
        datasets[[i]][[j]][[k]][na_indices] <- 0
        print(paste("NA values removed in", k, "for", j, "in", i))
      }
    }
  }
}


#code to create summary table of metadata

metadata_cols <- final_table[, 2:17]

# Step 2: Create an empty summary table with unique study names as rows and metadata columns as columns
summary_table <- data.frame(study_name = unique(final_table$study_name), stringsAsFactors = FALSE)

# Step 3: Loop through each metadata column and check for each unique study name
for (col in names(metadata_cols)) {
  for (study in unique(final_table$study_name)) {
    rows_for_study <- final_table$study_name == study
    if (any(!is.na(metadata_cols[rows_for_study, col]))) {
      summary_table[summary_table$study_name == study, col] <- 1
    } else {
      summary_table[summary_table$study_name == study, col] <- 0
    }
  }
}
# Step 4: Add a column to count the number of samples for each study
summary_table$sample_count <- sapply(unique(final_table$study_name), function(study) {
  sum(final_table$study_name == study)
})
# Step 5: Print the summary table
print(summary_table)
write.table(summary_table,"summary_table.txt",sep="\t",row.names = TRUE)

#Code to create a web application to plot data of Penguin(R inbuilt database)

# code to create loocv and oob matrix for all the studies in triglyceride in both the healthy and high class. 
library(randomForest)
library(dplyr)
oob_high_tg_levels <- matrix(ncol = 9, nrow = 9)
colnames(oob_high_tg_levels) <- names(high_tg_list)
rownames(oob_high_tg_levels) <- names(high_tg_list)
loocv_high_tg_levels <- matrix(ncol = 9, nrow = 1)
colnames(loocv_high_tg_levels) <- names(high_tg_list)
oob_healthy_tg_levels <- matrix(ncol = 9, nrow = 9)
colnames(oob_healthy_tg_levels) <- names(healthy_tg_list)
rownames(oob_healthy_tg_levels) <- names(healthy_tg_list)
loocv_healthy_tg_levels <- matrix(ncol = 9, nrow = 1)
colnames(loocv_healthy_tg_levels) <- names(healthy_tg_list)
rfmodels_tg_high <- list()
rfmodels_tg_healthy <- list()
rfmodel_same_high <- list()
rfmodel_diff_high <- list()
rfmodel_same_healthy <- list()
rfmodel_diff_healthy <- list()


loocv_high_subsampling_corrmat <- matrix(nrow = 25, ncol = 9)
colnames(loocv_high_subsampling_corrmat) <- names(high_tg_list)
loocv_healthy_subsampling_corrmat <- matrix(nrow = 25, ncol = 9)
colnames(loocv_healthy_subsampling_corrmat) <- names(healthy_tg_list)

healthy_studiesnames <- names(healthy_tg_list)
high_studynames <- names(high_tg_list)
num_subsampling <- 25

for(i in high_studynames){
  print(i)
  for(j in 1:num_subsampling){
    print(paste("Subsampling: ", j))
    study_leave_out <- i
    leave_out_data <- do.call(rbind, high_tg_list[-which(names(high_tg_list) == study_leave_out)])
    rfmodel <- randomForest(Triglyceride ~ ., data = leave_out_data, mtry = 13, ntree = 1000)
    rfmodels_tg_high[[i]][[j]] <- rfmodel
    print("rfmodel saved")
    predictions  <- predict( rfmodels_tg_high[[i]][[j]], newdata = high_tg_list[[i]])
    corr_value <- cor(predictions, high_tg_list[[i]][[3]])
    print(corr_value)
    loocv_high_subsampling_corrmat[j, i] <- corr_value
    print("Corr value stored")
    
    if(i == "YuJ_2015" & j == 25){
      for(k in names(high_tg_list)){
        for (l in names(high_tg_list)){
          if(k==l){
            current_data <- high_tg_list[[k]]
            rf_model <- randomForest(Triglyceride ~ ., data = current_data, ntree = 1000, mtry = 13)
            rfmodel_same_high[[i]] <- rf_model
            print("rfmodel stored")
            
            corr_value <- cor(rfmodel_same[[k]]$y, rfmodel_same[[k]]$predicted)
            print(paste("CorrValue:", k,"&", l,  corr_value, sep = "_"))
            oob_high_tg_levels[k,l] <- corr_value
            print("CorrValue added in the matrix")
          }else if(k != l){
            
            
            train_data <- high_tg_list[[k]]
            test_data <- high_tg_list[[l]]
            
            rfmodel <- randomForest(Triglyceride ~ ., data = train_data, ntree = 1000, mtry = 13)
            rfmodel_diff_high[[paste(k, "&", l, sep = "_")]] <- rfmodel
            print("rfmodel")
            predictions <- predict(rfmodel_diff[[paste(k, "&", l, sep = "_")]], newdata = test_data)
            corr_value <- cor(predictions, test_data$Triglyceride)
            print(paste("CorrValue:", k, "&", l, corr_value))
            oob_high_tg_levels[k,l] <- corr_value
            print("CorrValue added in the matrix")
          }
          
          
        }
        if(i =="YuJ_2015" & j == "YuJ_2015"){
          file_path <- "/storage/shivangi/aryan_test_code/oob_output/oob_healthy_env.RData"
          save.image(file = file_path)
          print("file saved")
          
        }
        
      }
      
    }
  }
  for(m in healthy_studiesnames){
    print(m)
    for(n in 1:num_subsampling){
      print(paste("Subsampling: ", n))
      study_leave_out <- m
      leave_out_data <- do.call(rbind, healthy_tg_list[-which(names(healthy_tg_list) == study_leave_out)])
      rfmodel <- randomForest(Triglyceride ~ ., data = leave_out_data, mtry = 13, ntree = 1000)
      rfmodels_tg_healthy[[m]][[n]] <- rfmodel
      print("rfmodel saved")
      predictions  <- predict(rfmodels_tg_healthy[[m]][[n]], newdata = healthy_tg_list[[m]])
      corr_value <- cor(predictions, healthy_tg_list[[m]][[3]])
      print(corr_value)
      loocv_healthy_subsampling_corrmat[n, m] <- corr_value
      print("Corr value stored")
      
      if(m =="YuJ_2015" & n == 25){
        for(o in names(healthy_tg_list)){
          for (p in names(healthy_tg_list)) {
            if(o==p){
              current_data <- healthy_tg_list[[o]]
              rf_model <- randomForest(Triglyceride ~ ., data = current_data, ntree = 1000, mtry = 13)
              rfmodel_same_healthy[[o]] <- rf_model
              print("rfmodel stored")
              
              corr_value <- cor(rfmodel_same[[k]]$y, rfmodel_same[[k]]$predicted)
              print(paste("CorrValue:", o,"&", p,  corr_value, sep = "_"))
              oob_healthy_tg_levels[o,p] <- corr_value
              print("CorrValue added in the matrix")
            }else if(o != p){
              train_data <- healthy_tg_list[[o]]
              test_data <- healthy_tg_list[[p]]
              
              rfmodel <- randomForest(Triglyceride ~ ., data = train_data, ntree = 1000, mtry = 13)
              rfmodel_diff_healthy[[paste(o, "&", p, sep = "_")]] <- rfmodel
              print("rfmodel")
              predictions <- predict(rfmodel_diff[[paste(o, "&", p, sep = "_")]], newdata = test_data)
              corr_value <- cor(predictions, test_data$Triglyceride)
              print(paste("CorrValue:", o, "&", p, corr_value))
              oob_healthy_tg_levels[o,p] <- corr_value
              print("CorrValue added in the matrix")
            }
            
          }
        }
      }
    }
  }
}




#code to create summary table 

for (i in rownames(summary_table)) {
  for (j in colnames(summary_table)) {
    if (any(is.na(final_table[final_table$study_name == i, j]))) {
      summary_table[i, j] <- 0
    } else {
      summary_table[i, j] <- 1
    }
  }
}


