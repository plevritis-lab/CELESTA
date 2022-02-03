#############################################################################################
#' Create CELESTA object
#' @export
Celesta <- setClass("Celesta",
                    slots = c(marker_exp_matrix = "matrix", # normalized expression from segmentation
                              original_exp ="matrix", # original expression from segmentation
                              prior_info = "data.frame", # store information from cell-type prior knowledge
                              cell_ID = "numeric", 
                              lineage_info = "data.frame",
                              coords = "matrix", #coordinates
                              marker_exp_prob = "matrix", # marker expression probability
                              cell_prob = "matrix", # cell type probability for each cell
                              final_cell_type_assignment = "matrix", 
                              project_name = "character",
                              nb_list = "matrix", # N-nearest neighbor list
                              total_rounds = "numeric",
                              cell_nb_in_bandwidth = "ANY", #Cells locates within a bandwidth to cell i
                              cell_nb_dist = "ANY", #The distance of each cell to cell i within a bandwidth
                              current_scoring_matrix = "matrix", #scoring function
                              initial_pri_matrix = "matrix", #user defined cell-type marker matrix
                              current_pri_matrix = "matrix", #updated cell-type marker matrix
                              current_cell_prob = "matrix", #cell probability for updates
                              current_cell_type_assignment = "matrix",
                              anchor_cell_type_assignment = "matrix",
                              starting_cell_type_assignment = "matrix",
                              current_beta = "matrix",
                              dist_from_nearest_assigned_cell = "matrix",
                              nb_cell_type = "ANY", #cell types of the neighboring cells for index cells
                              all_cell_nb_cell_type = "ANY", # cell types of the neighboring cells for all cells
                              unassigned_cells = "numeric", # store cells to check in each round and each iteration
                              assigned_cells = "numeric")) # cells already assigned cell type

#############################################################################################
#############################################################################################
#' Get protein marker expression
#' @export
GetMarkerExpMatrix <- function(CelestaObj,prior_marker_info,imaging_data_file,
                               cofactor,transform_type){
  markers_to_use <- colnames(prior_marker_info)[3:dim(prior_marker_info)[2]]
  matching_markers <- match(markers_to_use,colnames(imaging_data_file))
  if(length(which(is.na(matching_markers)==TRUE))>0){
    print("Please double check the protein markers in the cell-type marker matrix and 
          imaging input file")
  }else{
    if(transform_type==0){#no transform
      marker_exp_matrix <- data.matrix(imaging_data_file[,matching_markers])
      CelestaObj@marker_exp_matrix <- marker_exp_matrix
    }else if(transform_type==1){#arcsinh
      marker_exp_matrix <- data.matrix(imaging_data_file[,matching_markers])
      marker_exp_transformed <- asinh(marker_exp_matrix/cofactor)
      CelestaObj@marker_exp_matrix <- marker_exp_transformed
    }
    CelestaObj@original_exp <- data.matrix(imaging_data_file[,matching_markers])
    cellIDs <- seq(1,dim(marker_exp_matrix)[1],by=1)
    CelestaObj@cell_ID <- cellIDs
    return (CelestaObj)
  }
}
############################################################################################
#' Get prior knowledge on cell types
#' @export
GetPirorInfo <- function(CelestaObj,prior_marker_info){
  CelestaObj@prior_info <- prior_marker_info
  lineage_column <- prior_marker_info[,2]
  if(grepl("_", lineage_column[1], fixed = TRUE)){
    round <- integer()
    previous_cell_type <- integer()
    cell_type_number <- integer()
    for(i in 1:dim(prior_marker_info)[1]){
      info <- strtoi(unlist(strsplit(prior_marker_info[i,2],"_")))
      round[i] <- info[1]
      previous_cell_type[i] <- info[2]
      cell_type_number[i] <- info[3]
    }
    CelestaObj@lineage_info <- data.frame(Round=round,
                                          Previous_cell_type=previous_cell_type,
                                          Cell_type_number=cell_type_number)
    total_rounds <- max(CelestaObj@lineage_info$Round)
    CelestaObj@total_rounds <- total_rounds
  }else{
    print("Warning:the lineage information column has formatting errors")
  }
  return(CelestaObj)
}
#############################################################################################
#' Get coordinates
#' @export
GetCoords <- function(CelestaObj,imaging_data_file){
  Coords <- cbind(imaging_data_file$X,
                  imaging_data_file$Y)
  colnames(Coords) <- c("X","Y")
  CelestaObj@coords <- Coords
  return(CelestaObj)
}
############################################################################################
#' Gaussian mixture model for each marker
#' @export
GMM_fitting <- function(marker_exp,marker_name,figure=FALSE){
  print("Marker:")
  print(marker_name)
  GMM_marker_param <- matrix(nrow=3,ncol=2)
  set.seed(1)
  zero_indices <- which(marker_exp==0)
  zero_percentage <- length(zero_indices)/length(marker_exp)
  if(zero_percentage > 0.1 & zero_percentage<0.2){
    print("Warning: The marker expression potentially has too many zeros for fitting.
          GMM fitting will use input expression data with reduced sparsity")
    num_of_indices_to_remove <- floor(length(marker_exp)*(zero_percentage))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp,2,
                         models=mixmodGaussianModel(family="general",
                                                    listModels = "Gaussian_p_Lk_Ck",
                                                    free.proportions = FALSE,equal.proportions = TRUE))
    ### Check the models information for the Gaussian models, which shows which parameters are constrained. 
    ### Want equal proportions of the two Gaussians
    GMM_marker_param[1,] <- xxx@results[[1]]@parameters@proportions
    GMM_marker_param[2,] <- xxx@results[[1]]@parameters@mean[,1]
    GMM_marker_param[3,1] <- xxx@results[[1]]@parameters@variance[[1]][,1]
    GMM_marker_param[3,2] <- xxx@results[[1]]@parameters@variance[[2]][,1]
  }else if(zero_percentage >= 0.2 & zero_percentage<0.5){
    print("Warning: The marker expression potentially has too many zeros for fitting.
          GMM fitting will use input expression data with reduced sparsity")
    num_of_indices_to_remove <- floor(length(marker_exp)*(zero_percentage - 0.05))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp,2,
                         models=mixmodGaussianModel(family="general",
                                                    listModels = "Gaussian_p_Lk_Ck",
                                                    free.proportions = FALSE,equal.proportions = TRUE))
    ### Check the models information for the Gaussian models, which shows which parameters are constrained. 
    ### Want equal proportions of the two Gaussians
    GMM_marker_param[1,] <- xxx@results[[1]]@parameters@proportions
    GMM_marker_param[2,] <- xxx@results[[1]]@parameters@mean[,1]
    GMM_marker_param[3,1] <- xxx@results[[1]]@parameters@variance[[1]][,1]
    GMM_marker_param[3,2] <- xxx@results[[1]]@parameters@variance[[2]][,1]
  }else if(zero_percentage>=0.5 & zero_percentage<=0.9){
    print("Warning: The marker expression potentially has too many zeros for fitting.
          GMM fitting will use input expression data with reduced sparsity")
    num_of_indices_to_remove <- ceiling(length(marker_exp)*(zero_percentage-0.02))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp,2,
                         models=mixmodGaussianModel(family="general",
                                                    listModels = "Gaussian_p_Lk_Ck",
                                                    free.proportions = FALSE,equal.proportions = TRUE))
    ### Check the models information for the Gaussian models, which shows which parameters are constrained. 
    ### Want equal proportions of the two Gaussians
    GMM_marker_param[1,] <- xxx@results[[1]]@parameters@proportions
    GMM_marker_param[2,] <- xxx@results[[1]]@parameters@mean[,1]
    GMM_marker_param[3,1] <- xxx@results[[1]]@parameters@variance[[1]][,1]
    GMM_marker_param[3,2] <- xxx@results[[1]]@parameters@variance[[2]][,1]
  }else if(zero_percentage>=0.9){
    print("Warning: The marker expression potentially has too many zeros for fitting.
          GMM fitting will use input expression data with reduced sparsity")
    marker_exp <- marker_exp[-zero_indices]
    xxx <- mixmodCluster(marker_exp,2,
                         models=mixmodGaussianModel(family="general",
                                                    listModels = "Gaussian_p_Lk_Ck",
                                                    free.proportions = FALSE,equal.proportions = TRUE))
    ### Check the models information for the Gaussian models, which shows which parameters are constrained. 
    ### Want equal proportions of the two Gaussians
    GMM_marker_param[1,] <- xxx@results[[1]]@parameters@proportions
    GMM_marker_param[2,] <- xxx@results[[1]]@parameters@mean[,1]
    GMM_marker_param[3,1] <- xxx@results[[1]]@parameters@variance[[1]][,1]
    GMM_marker_param[3,2] <- xxx@results[[1]]@parameters@variance[[2]][,1]
  }else{
    xxx <- mixmodCluster(marker_exp,2,
                         models=mixmodGaussianModel(family="general",
                                                    listModels = "Gaussian_p_Lk_Ck",
                                                    free.proportions = FALSE,equal.proportions = TRUE))
    ### Check the models information for the Gaussian models, which shows which parameters are constrained. 
    ### Want equal proportions of the two Gaussians
    GMM_marker_param[1,] <- xxx@results[[1]]@parameters@proportions
    GMM_marker_param[2,] <- xxx@results[[1]]@parameters@mean[,1]
    GMM_marker_param[3,1] <- xxx@results[[1]]@parameters@variance[[1]][,1]
    GMM_marker_param[3,2] <- xxx@results[[1]]@parameters@variance[[2]][,1]
  }
  #print(GMM_marker_param)
  if(figure == TRUE){
    bin_size <- 20
    filename <- paste0(marker_name,"_GMM.png")
    png(filename,width = 5.5, height = 6.5,units = 'in',res = 300)
    h<-hist(marker_exp,breaks=bin_size,xlab="Marker expression",main=paste0("Histogram for ",marker_name))
    highestCount <- max(h$counts)
    multiplier <- h$counts/h$density
    xfit <- seq(min(marker_exp),max(marker_exp),length=length(h$breaks))
    yfit1 <- dnorm(xfit,mean=GMM_marker_param[2,1],sd=sqrt(GMM_marker_param[3,1]))*multiplier[1]
    lines(xfit, yfit1, col="blue", lwd=2)
    yfit2 <- dnorm(xfit,mean=GMM_marker_param[2,2],sd=sqrt(GMM_marker_param[3,2]))*multiplier[1]
    lines(xfit, yfit2, col="red", lwd=2)
    dev.off()
  }
  return(GMM_marker_param)
}
#############################################################################################
#' Build sigmoid function for calculation of expression probability
#' @export
build_sigmoid_function <- function(marker_exp_matrix,figure=FALSE){
  sigmoid_function_parameter <- matrix(nrow=2,ncol=dim(marker_exp_matrix)[2])
  ### For each marker, fit GMM
  for(i in 1:dim(marker_exp_matrix)[2]){
    marker_exp <- marker_exp_matrix[,i]
    marker_name <- colnames(marker_exp_matrix)[i]
    if(typeof(marker_name) != "character"){
      print("Protein marker name in the marker expression matrix has potential problem.")
    }else{
      marker_GMM_model <- GMM_fitting(marker_exp,marker_name,figure)
      weight <- marker_GMM_model[1,]
      mus <- marker_GMM_model[2,]
      sigmas <- marker_GMM_model[3,]
      
      if(mus[1] > mus[2]){ # first Gaussian model is for marker expressed, second is for marker not expressed
        a <- (-0.5 / sigmas[2] + 0.5 /sigmas[1])
        b <- mus[2] / sigmas[2] - mus[1] / sigmas[1]
        c <- 0.5 * (-mus[2]^2 / sigmas[2] + mus[1]^2 / sigmas[1]) + log(weight[2] / weight[1]) + 0.5 * log(sigmas[1] / sigmas[2])
        xroot <- (-b - sqrt(b^2 - 4.0 * a * c) ) / (2.0 * a)
        #slope <- 0.5 * (xroot - mus[2]) / sigmas[2] - 0.5 * (xroot - mus[1]) / sigmas[1]
        slope <- 1
      }else{# second Gaussian model is for marker expressed, first is for marker not expressed
        a <- (-0.5 / sigmas[1] + 0.5 /sigmas[2])
        b <- mus[1] / sigmas[1] - mus[2] / sigmas[2]
        c <- 0.5 * (-mus[1]^2 / sigmas[1] + mus[2]^2 / sigmas[2]) + log(weight[1] / weight[2]) + 0.5 * log(sigmas[2] / sigmas[1])
        xroot <- (-b - sqrt(b^2 - 4.0 * a * c) ) / (2.0 * a)
        #slope <- 0.5 * (xroot - mus[1]) / sigmas[1] - 0.5 * (xroot - mus[2]) / sigmas[2]
        slope <- 1
      }
      if(figure==TRUE){
        filename <- paste0(marker_name,"_sigmoid.png")
        ### plot sigmoid function
        exp_term <- exp(slope*(marker_exp-xroot))
        yyy <- exp_term/(1+exp_term)
        yyy <- (yyy-min(yyy))/(max(yyy)-min(yyy))
        png(filename,width = 4.5, height = 4.5,units = 'in',res = 300)
        plot(marker_exp, yyy, col = "darkblue",
             xlab = "", ylab = "", main = paste0(marker_name," sigmoid function"))
        grid()
        dev.off()
      }
      sigmoid_function_parameter[1,i] <- xroot
      sigmoid_function_parameter[2,i] <- slope
    }
  }
  return(sigmoid_function_parameter)
}
#############################################################################################
#' Calculate expression probability for each marker in the prior matrix
#' @export
marker_exp_probability <- function(CelestaObj,figure=FALSE){
  ### Fit GMM model and get parameters for the activation probabilities
  marker_exp_matrix <- CelestaObj@marker_exp_matrix
  sigmoid_function_parameter <- build_sigmoid_function(marker_exp_matrix,figure)
  ### Marker activation probability matrix
  marker_exp_prob <- matrix(nrow=dim(marker_exp_matrix)[1],ncol=dim(marker_exp_matrix)[2])
  colnames(marker_exp_prob) <- colnames(marker_exp_matrix)
  
  for(i in 1:dim(marker_exp_matrix)[2]){
    exp_term <- exp(sigmoid_function_parameter[2,i]*(marker_exp_matrix[,i]-sigmoid_function_parameter[1,i]))
    y = exp_term/(1+exp_term)
    marker_exp_prob[,i] <- (y-min(y))/(max(y)-min(y))
  }
  CelestaObj@marker_exp_prob <- marker_exp_prob
  return(CelestaObj)
}
#############################################################################################
#' Get neighborhood informtion
#' @export
GetNeighborInfo <- function(CelestaObj,number_of_neighbors=5,bandwidth=100){
  coords <- CelestaObj@coords
  print("Get nearest neighbors.")
  xxx <- knearneigh(coords,k=number_of_neighbors)
  nb_list <- xxx$nn
  colnames(nb_list) <- paste0("neighbor",seq(1,number_of_neighbors,by=1))
  ### Identify N-nearest neighbors for each cell
  CelestaObj@nb_list <- nb_list
  ### Identify cells within a circle bandwidth
  print("Identify neighboring cells within a defined bandwidth.")
  all_cell_nb_in_bandwidth <- dnearneigh(coords, 0, bandwidth, longlat = NULL)
  CelestaObj@cell_nb_in_bandwidth <- all_cell_nb_in_bandwidth
  ### Identify distances for all the cells within the circle
  CelestaObj@cell_nb_dist <- nbdists(all_cell_nb_in_bandwidth, coords)
  return(CelestaObj)
}
##############################################################################
#' Initialize the celesta object
#' @export
initialize_object <- function(CelestaObj){
  total_rounds <- max(CelestaObj@lineage_info$Round)
  current_cell_type_assignment <- matrix(0L,nrow =dim(CelestaObj@marker_exp_matrix)[1],
                                         ncol=total_rounds)
  CelestaObj@current_cell_type_assignment <- current_cell_type_assignment
  CelestaObj@anchor_cell_type_assignment <- current_cell_type_assignment
  CelestaObj@starting_cell_type_assignment <- current_cell_type_assignment
  
  current_scoring_matrix <- matrix(nrow=dim(CelestaObj@marker_exp_matrix)[1],
                                   ncol = dim(CelestaObj@prior_info)[1])
  colnames(current_scoring_matrix) <- CelestaObj@prior_info[,1]
  CelestaObj@current_scoring_matrix <- current_scoring_matrix
  
  current_cell_prob <- matrix(nrow=dim(CelestaObj@marker_exp_matrix)[1],
                              ncol = dim(CelestaObj@prior_info)[1])
  colnames(current_cell_prob) <- CelestaObj@prior_info[,1]
  CelestaObj@current_cell_prob <- current_cell_prob
  return(CelestaObj)
}
#############################################################################################
#' Create CELESTA object
#' @export
CreateCELESTAobj <- function(project_title="Project",prior_marker_info,imaging_data_file,
                             cofactor=10,transform_type=1,
                             number_of_neighbors=5,bandwidth=100){
  CelestaObj <- Celesta(project_name = project_title) 
  ### Get protein marker expressions and cell IDs
  CelestaObj <- GetMarkerExpMatrix(CelestaObj,prior_marker_info,imaging_data_file,cofactor=10,
                                   transform_type = transform_type)
  ### Get user-defined prior knowledge matrix and cell lineage information
  CelestaObj <- GetPirorInfo(CelestaObj,prior_marker_info)
  ### Get coordinates
  CelestaObj <- GetCoords(CelestaObj,imaging_data_file)
  ### Convert marker expressions to marker activation probability
  CelestaObj <- marker_exp_probability(CelestaObj)
  ### Get neighboring cell information
  #CelestaObj <- GetNeighborInfo(CelestaObj,number_of_neighbors=5,bandwidth=100)
  CelestaObj <- GetNeighborInfo(CelestaObj)
  #Initialize the matrices for scoring function and prob matrix
  CelestaObj <- initialize_object(CelestaObj)
  return(CelestaObj)
}
#############################################################################################
#############################################################################################
#' Filter cells that could potentially be artifacts
#' @export
cell_filtering <- function(high_marker_threshold=0.9, low_marker_threshold=0.4,
                           CelestaObj){
  ### Filter out cells that have marker expressions all high or all low
  total_rounds <- CelestaObj@total_rounds
  number_of_marker <- dim(CelestaObj@initial_pri_matrix)[2]
  for(i in 1:dim(CelestaObj@marker_exp_matrix)[1]){
    cell_activation_prob <- CelestaObj@marker_exp_prob[i,]
    if(MarkQuestionableCells(cell_activation_prob,high_marker_threshold,low_marker_threshold)){
      CelestaObj@current_cell_type_assignment[i,1:total_rounds] <- rep(NA,total_rounds)
    }else{
    }
  }
  CelestaObj@starting_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  return(CelestaObj)
}
################################################################################################
#' Mark questionable cells
#' @export
MarkQuestionableCells <- function(cell_activation_prob,high_marker_threshold,low_marker_threshold){
  number_of_marker <- length(cell_activation_prob)
  number_of_low_markers <- length(which(cell_activation_prob<low_marker_threshold))
  number_of_high_markers <- length(which(cell_activation_prob>high_marker_threshold))
  if(number_of_low_markers==number_of_marker | number_of_high_markers==number_of_marker){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
#############################################################################################
#' For each round, need to get the prior knowledge on the cell types
#' @export
get_initial_prior_matrix <- function(CelestaObj,round){
  lineage_info <- CelestaObj@lineage_info
  cell_type_num <- lineage_info$Cell_type_number[which(lineage_info$Round==round)]
  initial_pri_matrix <- data.matrix(prior_marker_info[which(lineage_info$Round==round),
                                                      3:dim(prior_marker_info)[2]])
  CelestaObj@initial_pri_matrix <- initial_pri_matrix
  CelestaObj@current_pri_matrix <- initial_pri_matrix
  return(CelestaObj)
}
################################################################################################
################################################################################################
#' Find cells to check
#' @export
find_unassigned_cells <- function(CelestaObj,round){
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  lineage_info <- CelestaObj@lineage_info
  cell_ID <- CelestaObj@cell_ID
  if(round == 1){
    unassigned_cells <- cell_ID[which(current_cell_type_assignment[,round] == 0)]
  }else{
    previous_level_type <- unique(lineage_info$Previous_cell_type[which(lineage_info$Round==round)])
    previous_level_round <- lineage_info$Round[which(lineage_info$Cell_type_number==previous_level_type)]
    unassigned_cells <- cell_ID[which(current_cell_type_assignment[,round] == 0 & 
                                      (current_cell_type_assignment[,previous_level_round]==previous_level_type))]
  }
  return(unassigned_cells)
}
################################################################################################
#' Find cells with ID assigned
#' @export
find_assigned_cells <- function(CelestaObj,round){
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  lineage_info <- CelestaObj@lineage_info
  cell_ID <- CelestaObj@cell_ID
  if(round == 1){
    assigned_cells <- cell_ID[which(current_cell_type_assignment[,round] != 0 & 
                                     is.na(current_cell_type_assignment[,round])==FALSE)]
  }else{
    previous_level_type <- unique(lineage_info$Previous_cell_type[which(lineage_info$Round==round)])
    assigned_cells <- cell_ID[which(current_cell_type_assignment[,round] != 0 & 
                                     is.na(current_cell_type_assignment[,round])==FALSE  & 
                                     (current_cell_type_assignment[,(round-1)]==previous_level_type))]
  }
  return(assigned_cells)
}
################################################################################################
#' Calculate scores using MSE
#' @export
get_score <- function(activation_prob_to_use,prior_info,non_NA_index){
  score <- apply(activation_prob_to_use[,non_NA_index],1,function(x) (1-sum((x-prior_info)^2)/length(x)))
  return(score)
}
#############################################################################################
#' Function for calculating scoring function
#' @export
scoring_function <- function(CelestaObj,round,unassigned_cells,cell_type_num){
  marker_exp_prob <- CelestaObj@marker_exp_prob
  current_pri_matrix <- CelestaObj@current_pri_matrix
  current_scoring_matrix <- CelestaObj@current_scoring_matrix
  print("Start calculating the scoring function.")
  activation_prob_to_use <- marker_exp_prob[unassigned_cells,]
  for(i in 1:length(cell_type_num)){
    non_NA_index <- which(!is.na(current_pri_matrix[i,]))
    prior_info <- current_pri_matrix[i,non_NA_index]
    current_scoring_matrix[unassigned_cells,cell_type_num[i]] <- get_score(activation_prob_to_use,prior_info,non_NA_index)
  }
  current_scoring_matrix[unassigned_cells,cell_type_num]<-t(apply(current_scoring_matrix[unassigned_cells,cell_type_num],
                                                                1,function(x) x/sum(x)))
  CelestaObj@current_scoring_matrix <- current_scoring_matrix
  return(CelestaObj)
}
################################################################################################
################################################################################################
#' Calculate probability differences
#' @export
find_min_prob_diff <- function(max.prob,max.prob_index,cell_prob_list,unassigned_cells){
  ### max.prob, max.prob_index are calculated only on unassigned_cells
  ### but cell_prob_list has all the cells
  min_prob_diff <- numeric(length=length(unassigned_cells))
  for(i in 1:length(unassigned_cells)){
    min_prob_diff[i] <- min(max.prob[i]-cell_prob_list[unassigned_cells[i],-max.prob_index[i]])
  }
  return(min_prob_diff)
}
################################################################################################
#' Find the cell types based on the scores (anchor cell) or probabilities (index cell)
#' @export
cell_type <- function(CelestaObj,cell_type_num,unassigned_cells,round,
                      min_difference=0,min_prob=0,
                      high_marker_threshold,low_marker_threshold){
  all_cell_prob <- CelestaObj@current_cell_prob
  initial_pri_matrix <- CelestaObj@initial_pri_matrix
  cell_prob_list <- all_cell_prob[,cell_type_num]
  cell_type_assignment <- CelestaObj@current_cell_type_assignment[,round]
  marker_exp_prob <- CelestaObj@marker_exp_prob
  max.prob_index <- apply(cell_prob_list[unassigned_cells,],1,which.max)
  max.prob <- apply(cell_prob_list[unassigned_cells,],1,max)
  min_prob_diff <- find_min_prob_diff(max.prob,max.prob_index,cell_prob_list,unassigned_cells)
  ### Find cells with cell type max probability > threshold and cell type probability difference > threshold
  ########################################
  ### Indexing on unassigned_cells!!!!!!!!
  threshold_cells <- unassigned_cells[which(min_prob_diff > min_difference & max.prob > min_prob)]
  max.prob_index_thresholded <- max.prob_index[which(min_prob_diff > min_difference & max.prob > min_prob)]
  ########################################
  for(i in 1:length(threshold_cells)){
    cell_ID_to_check <- threshold_cells[i]
    high_marker_index <- which(initial_pri_matrix[max.prob_index_thresholded[i],]==1)
    low_marker_index <- which(initial_pri_matrix[max.prob_index_thresholded[i],]==0)
    threshold_index <- cell_type_num[max.prob_index_thresholded[i]]
    if(length(which(marker_exp_prob[cell_ID_to_check,high_marker_index]>=high_marker_threshold[threshold_index]))==length(high_marker_index) &
       length(which(marker_exp_prob[cell_ID_to_check,low_marker_index]<=low_marker_threshold[threshold_index]))==length(low_marker_index)){
      cell_type_assignment[cell_ID_to_check] <- cell_type_num[max.prob_index_thresholded[i]]
    }else{
      #cell_type_assignment[cell_ID_to_check] <- 0
    }
  }
  CelestaObj@current_cell_type_assignment[,round] <- cell_type_assignment
  return(CelestaObj)
}
################################################################################################
#' Cell type count
#' @export
count_cell_type <- function(CelestaObj,cell_type_num,round){
  cell_type_count <- matrix(nrow=(length(cell_type_num)),ncol=3)
  colnames(cell_type_count) <- c("cell_type_number","count","proportion")
  prior_marker_info <- CelestaObj@prior_info
  row.names(cell_type_count) <- prior_marker_info[cell_type_num,1]
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  cell_type_count[,1] <- cell_type_num
  total_cell_number <- dim(current_cell_type_assignment)[1]
  for(i in 1:length(cell_type_num)){
    cell_type_count[i,2] <- length(which(current_cell_type_assignment[,round]==cell_type_num[i]))
    cell_type <- prior_marker_info[cell_type_num[i],1]
    if(cell_type_count[i,2]<1){
      print(paste0("Too few cells identified for: ",cell_type))
      print("Please consider relax threshold.")
    }
  }
  cell_type_count[,3] <- cell_type_count[,2]/total_cell_number
  return(cell_type_count)
}
################################################################################################
# plot_cells_iteration <- function(CelestaObj,cell_number_to_use,round,
#                        cell_type_colors,point_size=0.1,iteration,figure = FALSE){
#   if(figure==TRUE){
#     coords <- CelestaObj@coords
#     current_cell_type_assignment <- CelestaObj@current_cell_type_assignment[,round]
#     project_name <- CelestaObj@project_name
#     prior_marker_info <- CelestaObj@prior_info
#     cell_types <- prior_marker_info[,1]
#     x_min <- min(coords[,1])
#     x_max <- max(coords[,1])
#     y_min <- min(coords[,2])
#     y_max <- max(coords[,2])
#     range <- c(min(x_min,y_min),max(x_max,y_max))
#     
#     filename <- paste0(project_name,paste0(paste0("Round_",round),
#                                            paste0("_Iteration_",paste0(iteration,".png"))))
#     cell_index <- integer()
#     cell_anno <- character()
#     count <- 0
#     for(i in 1:length(cell_number_to_use)){
#       unassigned_cells <- which(current_cell_type_assignment == cell_number_to_use[i])
#       cell_index[(count+1):(count+length(unassigned_cells))] <- unassigned_cells
#       cell_anno[(count+1):(count+length(unassigned_cells))] <- cell_types[cell_number_to_use[i]]
#       count <- count + length(unassigned_cells)
#     }
#     df_plot <- data.frame(x=coords[cell_index,1],
#                           y=coords[cell_index,2],
#                           cell_anno=cell_anno)
#     df_plot$cell_anno <- factor(df_plot$cell_anno,levels = c(cell_types[cell_number_to_use]))
#     color_plot <- cell_type_colors[cell_number_to_use]
#     
#     g<- ggplot(df_plot,aes(x=x,y=y,group=cell_anno))+geom_point(aes(color=cell_anno),size=point_size)+
#       scale_color_manual(values=color_plot)+
#       xlim(range[1],range[2])+ylim(range[1],range[2])+
#       labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(), 
#                           panel.grid.minor = element_blank(),
#                           legend.title = element_blank(),
#                           legend.text=element_text(size=12,face = "bold"),
#                           panel.background = element_rect(fill = 'black'), 
#                           axis.line = element_line(colour = "black"),
#                           axis.title.x=element_blank(),
#                           axis.title.y=element_blank())+
#       guides(colour = guide_legend(override.aes = list(size=10)))
#     ggsave(filename,plot=g,width = 16.5, height = 16,units = 'in',dpi = 300)
#   }
# }
################################################################################################
#' Find the cell types of the neighbors for unassigned_cells
#' @export
neighbor_cell_type <- function(CelestaObj,cell_type_num,round,unassigned_cells){
  ### Only has information for cells to check
  nb_list <- CelestaObj@nb_list
  cell_type_assignment <- CelestaObj@current_cell_type_assignment[,round]
  same_type_nb <- matrix(rep(list(),length(cell_type_num)*length(unassigned_cells)),
                         nrow=length(unassigned_cells),ncol=length(cell_type_num))
  row.names(same_type_nb) <- unassigned_cells
  colnames(same_type_nb) <- cell_type_num
  for(j in 1:length(unassigned_cells)){
    current_cell_ID <- unassigned_cells[j]
    neighbors <- nb_list[current_cell_ID,]
    neighbor_types <- cell_type_assignment[neighbors]
    for(i in 1:length(cell_type_num)){
      same_type_nb[j,i][[1]] <- neighbors[which(neighbor_types == cell_type_num[i])]
    }
  }
  CelestaObj@nb_cell_type <- same_type_nb
  return(CelestaObj)
}
################################################################################################
################################################################################################
#' Get distance from nearest assigned cells
#' @export
get_dist_from_nearest_assigned_cells <- function(CelestaObj,cell_type_num,unassigned_cells,
                                                 assigned_cells,round){
  print("Get distance from nearest assigned cells.")
  all_cell_nb_in_circle <- CelestaObj@cell_nb_in_bandwidth
  all_cell_nb_circle_dist <- CelestaObj@cell_nb_dist
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  dist_nearest_assigned_cell <- matrix(nrow = length(unassigned_cells),
                                       ncol = length(cell_type_num))
  colnames(dist_nearest_assigned_cell) <- cell_type_num
  for(i in 1:dim(dist_nearest_assigned_cell)[1]){
    cell_to_check <- unassigned_cells[i]
    matching <- match(all_cell_nb_in_circle[[cell_to_check]],assigned_cells)
    index <- matching[which(is.na(matching)==FALSE)]
    if(length(index)==0){
      
    }else{
      nb_cell_with_ID <- assigned_cells[index]
      nb_cell_type <- current_cell_type_assignment[nb_cell_with_ID,round]
      unique_nb_cell_type <- unique(nb_cell_type)
      nb_cell_dist <- all_cell_nb_circle_dist[[cell_to_check]][which(is.na(matching)==FALSE)]
      for(j in 1:length(unique_nb_cell_type)){
        type_j <- which(nb_cell_type == unique_nb_cell_type[j])
        dist_nearest_assigned_cell[i,which(cell_type_num==unique_nb_cell_type[j])] <- min(nb_cell_dist[type_j])
      }
    }
  }
  CelestaObj@dist_from_nearest_assigned_cell <- dist_nearest_assigned_cell
  return(CelestaObj)
}
#############################################################################################
#' Function to calcualte beta
#' @export
calculate_beta <- function(CelestaObj,scale_factor=5,bandwidth=100){
  dist_from_nearest_assigned_cell <- CelestaObj@dist_from_nearest_assigned_cell
  beta <- scale_factor*(1-dist_from_nearest_assigned_cell/bandwidth)
  beta[is.na(beta)] <- 0
  CelestaObj@current_beta <- beta
  return(CelestaObj)
}
################################################################################################
#' Function to calculate probability for index cells
#' @export
cell_prob <- function(CelestaObj,cell_type_num,unassigned_cells,round){
  # This function uses mean field estimation to calculate probability
  # For each cell, a probability is calculated for each cell type to check
  current_cell_prob_list <- CelestaObj@current_cell_prob[,cell_type_num] #all cells*cell_type_num
  u <- current_cell_prob_list #all cells
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment[,round] # all cells
  current_beta <- CelestaObj@current_beta #cells to check
  nb_cell_type <- CelestaObj@nb_cell_type #cells to check
  current_scoring_matrix <- CelestaObj@current_scoring_matrix # all cells*all cell types
  for(i in 1:length(unassigned_cells)){
    cell_ID_to_check <- unassigned_cells[i]
    u_i <- numeric(length=length(cell_type_num))
    number_of_nb <- lengths(nb_cell_type[i,])
    for(j in 1:length(number_of_nb)){
      current_same_type_nb <- unlist(nb_cell_type[i,j][[1]])
      u_i[j] <- exp(current_scoring_matrix[cell_ID_to_check,cell_type_num[j]])*
        exp(current_beta[i,j]*
              sum(current_cell_prob_list[current_same_type_nb,j]))
    }
    u[cell_ID_to_check,] <- u_i/sum(u_i)
  }
  print("Cell probability updating done.")
  CelestaObj@current_cell_prob[,cell_type_num] <- u
  return(CelestaObj)
}
################################################################################################
#' Function to update prior knowledge matrix of the cell type signatures
#' @export
update_prior_matrix <- function(CelestaObj,round,cell_type_num){
  updated_prior_matrix <- CelestaObj@current_pri_matrix
  initial_pri_matrix_data <- CelestaObj@initial_pri_matrix
  current_pri_matrix_data <- CelestaObj@current_pri_matrix
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  all_marker_pro_matrix <- CelestaObj@marker_exp_prob
  for(i in 1:length(cell_type_num)){ ### for each cell type
    cell_type_to_check <- cell_type_num[i]
    for(j in 1:dim(current_pri_matrix_data)[2]){
      if(is.na(initial_pri_matrix_data[i,j])==TRUE){
      }else{
        cells_of_current_cell_type <- which(current_cell_type_assignment[,round] == cell_type_to_check)
        updated_prior_matrix[i,j] <- (mean(all_marker_pro_matrix[cells_of_current_cell_type,j])+
                                        initial_pri_matrix_data[i,j])/2
      }
    }
  }
  CelestaObj@current_pri_matrix <- updated_prior_matrix
  return(CelestaObj)
}
##############################################################################
### For different rounds
# plot_marker_exp <- function(CelestaObj,cell_type_colors=c(palette()[2:7],"white"),
#                             cell_type_num,round){
#     sample_name <- CelestaObj@project_name
#     current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
#     marker_exp_matrix <- CelestaObj@marker_exp_matrix
#     plot_matrix <- matrix(0L,nrow=length(cell_type_num),ncol=dim(marker_exp_matrix)[2])
#     row.names(plot_matrix) <- as.character(CelestaObj@prior_info[cell_type_num,1])
#     colnames(plot_matrix) <- colnames(marker_exp_matrix)
#     for(i in 1:length(cell_type_num)){
#       plot_matrix[i,] <- colMeans(marker_exp_matrix[which(current_cell_type_assignment[,round] == cell_type_num[i]),])
#     }
#     df <- as.data.frame(cbind(row.names(plot_matrix),plot_matrix))
#     colnames(df) <- c("cell_types",colnames(marker_exp_matrix))
#     df.m <- melt(df, id.var = "cell_types")
#     
#     df.m$value <- as.numeric(df.m$value)
#     df.m$cell_types <- factor(df.m$cell_types,levels = row.names(plot_matrix))
#     
#     filename <- paste0(sample_name,"_")
#     filename1 <- paste0(filename,round)
#     filename2 <- paste0(filename1,"_ave_marker_exp.png")
#     
#     g<-ggplot(df.m,aes(x=variable,y=value,group=cell_types,color=cell_types)) + geom_point() + geom_line() +
#       scale_color_manual(values=cell_type_colors[cell_type_num])+xlab("Marker")+
#       ylab("Expression")+theme_bw()+
#       theme(legend.title = element_blank())+
#       theme(axis.text.x = element_text(angle = 80, hjust = 1,size=12,face="bold"),
#             legend.text=element_text(size=12,face = "bold"),
#             panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#     ggsave(filename2,plot=g,width=13.5,height=9,units = 'in',dpi = 300)
# }
################################################################################################
################################################################################################
#' Get final results
#' @export
get_final_inferred_cell_types <- function(total_rounds,CelestaObj,imaging_data){
  current_cell_type_assignment <- CelestaObj@current_cell_type_assignment
  anchor_cell_assignment <- CelestaObj@anchor_cell_type_assignment
  cell_type_name_assigned <- matrix(nrow=dim(current_cell_type_assignment),ncol=total_rounds)
  anchor_cell_type_name_assigned <- matrix(nrow=dim(current_cell_type_assignment),ncol=total_rounds)
  prior_marker_info <- CelestaObj@prior_info
  lineage_info <- CelestaObj@lineage_info
  final_cell_type_assignment <- rep(0,length=dim(current_cell_type_assignment)[1])
  for(i in 1:total_rounds){
    current_pri_matrix_num <- i
    cell_type_name_assigned[,i] <- prior_marker_info[match(current_cell_type_assignment[,i],
                                                           lineage_info$Cell_type_number),1]
    cell_type_name_assigned[which(current_cell_type_assignment[,i]==0),i] <- "Unknown"
    anchor_cell_type_name_assigned[,i] <- prior_marker_info[match(anchor_cell_assignment[,i],
                                                                  lineage_info$Cell_type_number),1]
    anchor_cell_type_name_assigned[which(anchor_cell_assignment[,i]==0),i] <- "Unknown"
    if(current_pri_matrix_num == 1){
      final_cell_type_assignment <- current_cell_type_assignment[,current_pri_matrix_num]
    }else{
      previous_level_type <- unique(lineage_info$Previous_cell_type[which(lineage_info$Round==current_pri_matrix_num)])
      assignment <- current_cell_type_assignment[which(final_cell_type_assignment==previous_level_type &
                                                         current_cell_type_assignment[,current_pri_matrix_num]!=0),current_pri_matrix_num]
      final_cell_type_assignment[which(final_cell_type_assignment==previous_level_type &
                                         current_cell_type_assignment[,current_pri_matrix_num]!=0)] <- assignment
    }
  }
  final_cell_names <- character(length=dim(current_cell_type_assignment)[1])
  final_cell_names <- prior_marker_info[match(final_cell_type_assignment,lineage_info$Cell_type_number),1]
  final_cell_names[which(final_cell_type_assignment==0)] <- "Unknown"
  final_result <- cbind(cell_type_name_assigned,final_cell_type_assignment,final_cell_names)
  round_name <- paste("Round",seq(1,total_rounds,by=1))
  colnames(final_result) <- c(round_name,"Cell type number","Final cell type")
  filename <- paste0(CelestaObj@project_name,"_final_cell_type_assignment.csv")
  write.csv(cbind(final_result,imaging_data),file=filename,row.names = FALSE)
  filename <- paste0(CelestaObj@project_name,"_anchor_cell_assignment.csv")
  write.csv(anchor_cell_type_name_assigned,file=filename)
  CelestaObj@final_cell_type_assignment <- final_result
  return(CelestaObj)
}
#############################################################################################
#############################################################################################
#' Plot the cells using XY coordinates
#' @export
plot_cells_any_combination <- function(cell_type_assignment_to_plot,CelestaObj,
                                       cell_number_to_use,
                                       cell_type_colors=c(palette()[2:7],"white"),
                                       test_size=1){
  ### Cannot plot more than 7 cell types
  current_cell_type_assignment <- cell_type_assignment_to_plot
  coords <- CelestaObj@coords
  cell_types <- CelestaObj@prior_info[cell_number_to_use,1]
  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))
  
  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for(i in 1:length(cell_number_to_use)){
    unassigned_cells <- which(current_cell_type_assignment == cell_number_to_use[i])
    cell_index[(count+1):(count+length(unassigned_cells))] <- unassigned_cells
    cell_anno[(count+1):(count+length(unassigned_cells))] <- cell_types[i]
    count <- count + length(unassigned_cells)
  }
  df_plot <- data.frame(x=coords[cell_index,1],
                        y=coords[cell_index,2],
                        cell_anno=cell_anno)
  df_plot$cell_anno <- factor(df_plot$cell_anno,levels = c(cell_types))
  color_plot <- cell_type_colors[1:length(cell_number_to_use)]
  
  g<- ggplot(df_plot,aes(x=x,y=y,group=cell_anno))+geom_point(aes(color=cell_anno),size=test_size)+
    scale_color_manual(values=color_plot)+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        panel.background = element_rect(fill = 'black'), 
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.text = element_text(size=12,face="bold") )+
    guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave(filename="plot_cell_assignment.png",plot=g,width = 12, height = 12,units = 'in',dpi = 300)
}
#############################################################################################
#' Plot the expression probabilities of cells in the tissue
#' @export
plot_exp_prob <- function(CelestaObj,size_to_use=1,width_to_use=5,height_to_use=4){
  coords <- CelestaObj@coords
  marker_exp_prob <- CelestaObj@marker_exp_prob
  prior_marker_info <- CelestaObj@prior_info
  palette <- colorRampPalette(colors=c("white", "blue4"))
  cols <- palette(6)
  #plot(1:6, col=cols, pch=16, cex=3)
  
  markers_to_check <- as.character(colnames(prior_marker_info)[3:dim(prior_marker_info)[2]])
  for(i in 1:length(markers_to_check)){
    marker_to_use <- markers_to_check[i]
    marker_exp_prob_to_use <- marker_exp_prob[,which(colnames(marker_exp_prob)==marker_to_use)]
    cols_anno <- character(length=length(marker_exp_prob_to_use))
    cols_anno[which(marker_exp_prob_to_use>0.9)] <- ">0.9"
    cols_anno[which(marker_exp_prob_to_use>0.8 & marker_exp_prob_to_use<=0.9)] <- ">0.8"
    cols_anno[which(marker_exp_prob_to_use>0.7 & marker_exp_prob_to_use<=0.8)] <- ">0.7"
    cols_anno[which(marker_exp_prob_to_use>0.5 & marker_exp_prob_to_use<=0.7)] <- ">0.5"
    cols_anno[which(marker_exp_prob_to_use<=0.5)] <- "<=0.5"
    
    mca <- data.frame(Coords_1 = round(coords[,1],digits = 2),
                      Coords_2 = round(coords[,2],digits = 2),
                      Exp_quantile = round(marker_exp_prob_to_use,digits = 2),
                      Col_anno=cols_anno)
    row.names(mca) <- NULL
    colnames(mca) <- c("X","Y","Expression","Color_anno")
    mca$Color_anno <- factor(mca$Color_anno,levels=c("<=0.5",">0.5",">0.7",">0.8",">0.9"))
    
    x_min <- min(coords[,1])
    x_max <- max(coords[,1])
    y_min <- min(coords[,2])
    y_max <- max(coords[,2])
    range <- c(min(x_min,y_min),max(x_max,y_max))
    
    filename <- paste0(marker_to_use,"_exp_prob.png")
    g <- ggplot(mca,aes(x=X,y=Y,color=Color_anno)) + 
      xlim(range[1],range[2])+ylim(range[1],range[2])+
      geom_point(shape=20,size=size_to_use) +
      ggtitle(marker_to_use)+theme_bw()+
      scale_colour_manual(values=c(cols[1],cols[2],cols[3],cols[4],cols[6]))+
      #scale_colour_manual(values=cols)+
      theme(legend.title = element_blank(),
            legend.text = element_text(size=14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
      guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename,plot=g,width=width_to_use,height=height_to_use,units = 'in',dpi = 300)
  }
}
#############################################################################################
#' Function to assign cell types through iterations
#' @export
assign_cell_main <- function(CelestaObj,max_iteration=10,cell_change_threshold=0.01,
                             min_diff=0,min_probability=0,
                             high_marker_threshold_anchor=rep(0.7,length=50),
                             low_marker_threshold_anchor=rep(0.9,length=50),
                             high_marker_threshold_iteration=rep(0.5,length=50),
                             low_marker_threshold_iteration=rep(1,length=50)){
  total_rounds <- CelestaObj@total_rounds
  ### This loop is the main part for cell type assignment
  ### Cell type assignment function (normally should finish within 10min for ~100k cells)
  ### It runs pretty fast for below 50k cells
  for(i in 1:total_rounds){
    round <- i
    CelestaObj@current_cell_type_assignment[,round] <- CelestaObj@starting_cell_type_assignment[,round]
    
    current_number_of_cells_changed <- numeric()
    loglikelihood <- numeric()
    lineage_info <- CelestaObj@lineage_info
    cell_type_num <- lineage_info$Cell_type_number[which(lineage_info$Round==round)]
    CelestaObj <- get_initial_prior_matrix(CelestaObj,round)
    unassigned_cells <- find_unassigned_cells(CelestaObj,round)
    number_of_cells_to_find_identity <- length(unassigned_cells)
    print(number_of_cells_to_find_identity)
    ### Get scoring function 
    CelestaObj <- scoring_function(CelestaObj,round,unassigned_cells,cell_type_num)
    ### Initialize the cell probability with initial scores
    CelestaObj@current_cell_prob <- CelestaObj@current_scoring_matrix
    ### Assign anchor cells
    old_cell_assignment <- CelestaObj@current_cell_type_assignment[,round]
    CelestaObj <- cell_type(CelestaObj,cell_type_num,unassigned_cells,round,
                            min_difference=min_diff,
                            min_prob=min_probability,
                            high_marker_threshold=high_marker_threshold_anchor,
                            low_marker_threshold=low_marker_threshold_anchor)
    # 
    iteration <- 1
    CelestaObj@anchor_cell_type_assignment[,round] <- CelestaObj@current_cell_type_assignment[,round]
    print(cell_type_count <- count_cell_type(CelestaObj,cell_type_num,round))
    if(length(which(cell_type_count[,2]<1))==length(cell_type_num)){
      print("Too few cells identified for certain cell type, please consider relaxing threshold.")
      return(CelestaObj)
      break
    }
    current_number_of_cells_changed[iteration] <- 1
    #############
    ### Find cells to check
    unassigned_cells <- find_unassigned_cells(CelestaObj,round)
    assigned_cells <- find_assigned_cells(CelestaObj,round)
    #############
    ### Calculate beta
    CelestaObj <- neighbor_cell_type(CelestaObj,cell_type_num,round,unassigned_cells)
    CelestaObj <- get_dist_from_nearest_assigned_cells(CelestaObj,cell_type_num,
                                                       unassigned_cells,assigned_cells,round)
    CelestaObj <- calculate_beta(CelestaObj,scale_factor = 5,bandwidth = 100)
    ### Iterative cell type assignment
    while(iteration < max_iteration & current_number_of_cells_changed[iteration] > cell_change_threshold){
      iteration <- iteration + 1
      ### Calculate cell type probabilities
      CelestaObj <- cell_prob(CelestaObj, cell_type_num,unassigned_cells,round)
      old_cell_assignment <- CelestaObj@current_cell_type_assignment[,round]
      ### Update cell type assignment
      CelestaObj <- cell_type(CelestaObj,cell_type_num,unassigned_cells,round,
                              min_difference=min_diff,
                              min_prob=min_probability,
                              high_marker_threshold=high_marker_threshold_iteration,
                              low_marker_threshold=low_marker_threshold_iteration)
      print(cell_type_count <- count_cell_type(CelestaObj,cell_type_num,round))
      current_number_of_cells_changed[iteration] <- length(which((old_cell_assignment-CelestaObj@current_cell_type_assignment[,round])!=0))/number_of_cells_to_find_identity
      print(current_number_of_cells_changed[iteration])
      if(current_number_of_cells_changed[iteration] < cell_change_threshold){
        break
      }
      #############
      ### Find cells to check
      unassigned_cells <- find_unassigned_cells(CelestaObj,round)
      assigned_cells <- find_assigned_cells(CelestaObj,round)
      if(length(unassigned_cells)==0){
        break
      }
      #############
      ### Calculate beta
      CelestaObj <- neighbor_cell_type(CelestaObj,cell_type_num,round,unassigned_cells)
      CelestaObj <- get_dist_from_nearest_assigned_cells(CelestaObj,cell_type_num,
                                                         unassigned_cells,assigned_cells,round)
      CelestaObj <- calculate_beta(CelestaObj,scale_factor = 5,bandwidth = 100)
      ############
      ### Update prior cell-type marker matrix
      CelestaObj <- update_prior_matrix(CelestaObj,round,cell_type_num)
      ### Update scoring function
      CelestaObj <- scoring_function(CelestaObj,round,unassigned_cells,cell_type_num)
    }
  }
  CelestaObj <- get_final_inferred_cell_types(total_rounds,CelestaObj,imaging_data)
  return(CelestaObj)
}
#############################################################################################
#############################################################################################