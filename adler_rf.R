library(pacman)
p_load(data.table, foreach, doParallel, rpart)

# data_long <- fread("test.csv")
# 
# test.data <- data_long[, .(Ri, ID, continuous_time, Z1, Z2, window, LOS.exacerbed, inc.LOS, event_in_subject), 
#                        c(grep("^X", names(data_long), value = TRUE))]
# ID.vec <- test.data$ID
# cov.mat <- test.data[, !c("Ri", "ID"), with = FALSE]

create_subject_df <- function(subject_id, unique_times) {
  return(data.table(ID = subject_id, continuous_time = unique_times))
}



# ok, rf seems to be doing ok.

calibrated.forest <- function(n.iter = 500, 
                              data = test.data,
                              chosen.seed = NULL,
                              time.name,
                              cov.mat,
                              ID.vec,
                              nsamp = NULL,
                              minbuck = NULL,
                              always.split = NULL,
                              loss.function) {
  
  if(is.null(chosen.seed)){
    chosen.seed = sample(1:10^5, size = 1)
  }
  if(is.null(nsamp)) {
    nsamp <- round(sqrt(ncol(cov.mat)))
  }
  nsub <- length(unique(ID.vec))
  if(is.null(minbuck)) {
    minbuck <- round(.1*nsub)
  }
  if (time.name != "continuous_time") {
    setnames(data, old = time.name, new = "continuous_time")
    always.split = c(setdiff(always.split, time.name), "continuous_time")
  }

  set.seed(chosen.seed)
  out <- foreach(l = 1:n.iter, 
                 .combine = rbind, 
                 .packages = c("data.table", "rpart", "matrixStats", "stats")
  ) %dopar% {
    
    create_subject_df <- function(subject_id, unique_times) {
      return(data.table(ID = subject_id, continuous_time = unique_times))
    }
    
    list_of_dfs <- data[, .(unique_times = list(unique(data[, continuous_time]))), by = ID]
    list_of_dfs[, subject_df := lapply(.I, function(i) create_subject_df(ID[i], unique_times[[i]]))]
    list_of_dfs[, "unique_times" := NULL]
    
    
    bagged <- sample(unique(ID.vec), size = nsub, replace = TRUE)
    n.unique = length(unique(bagged))
    bagged.df <- data.table(ID = bagged)[, .N, by = ID]
    
    single_time_point_per_id <- data[ID %in% bagged.df$ID,
                                     .SD[sample(.N, 1)], by = "ID"]
    
    # Ensure we have the correct replication using the N column from bagged.df
    setkey(single_time_point_per_id, ID)
    setkey(bagged.df, ID)
    
    # Merge and replicate rows based on count in bagged.df
    sampled_data <- single_time_point_per_id[bagged.df, {
      
      # Replicating while ensuring column retention (including ID)
      reps_list <- replicate(N, .SD, simplify = FALSE)
      
      # Bind the list into a data.table
      replicate_data <- rbindlist(reps_list)
      replicate_data[, ID := .BY$ID]  # Reassign ID directly to the results
      
      # Output the replicated data
      replicate_data
      
    }, by = .EACHI]
    
    # Sample covariates
    cov.samp <- c(sample(c(names(cov.mat), "continuous_time"), size = nsamp))
    if(!is.null(always.split)) {
      cov.samp <- unique(c(cov.samp, sample(always.split, 5)))
    }
    
    form <- as.formula(paste0("Ri ~ ", paste(cov.samp, collapse = "+")))
    
    ran.tree <- rpart(form, data = sampled_data, 
                      method = loss.function,
                      control = rpart.control(minbucket = minbuck,
                                              cp = -1))
    
    out.bag <- data[!ID %in% bagged]
    # paths <- path.rpart(ran.tree, nodes = rownames(ran.tree$frame)[ran.tree$frame$var == "<leaf>"])
    # 
    # # Create covariate matrices using data.table
    # create_covar_matrix <- function(data, paths) {
    #   sapply(paths, function(path) {
    #     as.integer(eval(parse(text=paste(setdiff(path, "root"), collapse = " & ")), envir=data))
    #   })
    # }
# 
#     in.bag = data[ID %in% bagged]
# 
#     bagged.covar_matrix <- create_covar_matrix(in.bag, paths)
#     out_bag_covar_matrix <- create_covar_matrix(out.bag, paths)
# 
#     ref.column <- which(colSums2(bagged.covar_matrix) == max(colSums2(bagged.covar_matrix)))
# 
#     init.glm <- glm(Ri ~ . , data = cbind(in.bag[, "Ri"], bagged.covar_matrix[, - ref.column]), model=FALSE)
#     # intercept.init <- init.glm$coefficients["(Intercept)"]
# 
#     intercept.calib_glm <- glm(Ri ~ ., data =cbind(out.bag[, "Ri"], out_bag_covar_matrix[, -ref.column] ))
# 
#     init.glm$coefficients["(Intercept)"] <- intercept.calib_glm$coefficients["(Intercept)"]

    # Predictions for out-of-bag
    if(loss.function == "class"){
      out.bag[, init.rf.pred := predict(ran.tree, newdata = out.bag, type = "vector")-1]
    }else{
      out.bag[, init.rf.pred := predict(ran.tree, newdata = out.bag, type = "vector")]
    }
    # 
    # out.bag[, calibed.preds := predict(init.glm, newdata = as.data.frame(out_bag_covar_matrix) , type = "response")]
    # 
    # 
    
    
    list_of_dfs[, subject_df := lapply(subject_df, function(df) {
      
      if (!anyNA(out.bag[["init.rf.pred"]]) #& !anyNA(out.bag[["calibed.preds"]]) 
          ) { # If no NAs are present, proceed with the merge
        setkey(df, "ID", "continuous_time")
        setkey(out.bag, "ID", "continuous_time")
        df <- out.bag[df, .SD, .SDcols = c("ID", "continuous_time", "init.rf.pred") #, "calibed.preds")
                      , nomatch = 0]
        df = df[, num.tree.in:= l]
      }
      # Return df unchanged if there's any NA, ensure it's a data.table
      setDT(df)
      return(df)
    })]
    
    filtered_subject_df_list <- Filter(function(dt) nrow(dt) > 0, list_of_dfs$subject_df)
    out_one = rbindlist(filtered_subject_df_list)
    # sink("log.txt", append=TRUE)
    # print(paste("Iteration:", i, "Output Rows:", nrow(out_one)))
    # sink()
    out_one
    
    
    
  }
  # Aggregate predictions
  out_avg = out[, .(
    avg.rf.pred = mean(init.rf.pred)#,
    # avg.calib.pred = mean(calibed.preds)
    ), by = .(ID, continuous_time)]
  
  return(out_avg)
  
}

