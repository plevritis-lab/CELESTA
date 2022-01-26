################################################################################
################################################################################
#' Celesta
#'
#' @description Celesta object definition
#'
#' @slot project_name name of the project (used in file names)
#' @slot prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @slot marker_exp_matrix transformed protein marker expression (or original
#' segmentation protein marker expression if transformation is not specified)
#' @slot original_exp original protein marker expression (containing only the
#' protein markers specified in `prior_info`)
#' @slot cell_ID the IDs of the cells (from 1 to the total number of cells)
#' @slot lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @slot coords the x, y coordinates of each cell
#' @slot cell_prob cell type probability for each cell
#' @slot final_cell_type_assignment the final cell type assignments
#' @slot nb_list the list of N-nearest neighbors
#' @slot total_rounds the maximum round value
#' @slot cell_nb_in_bandwidth the cells located within a bandwidth to cell *c*
#' @slot cell_nb_dist the distance of each cell to cell *c* within a bandwidth
#' @slot initial_pri_matrix user defined cell-type marker matrix for a specific
#' round
#' @slot anchor_cell_type_assignment the anchor cell type assignments
#' @slot dist_from_nearest_assigned_cell the distance from the nearest assigned
#' cell
#' @slot nb_cell_type cell types of the neighboring cells for index cells
#'
#' @slot marker_exp_prob the marker expression probability for each cell
#' @slot current_scoring_matrix the current scoring matrix
#' (number_cells x number_cell_type)
#' @slot current_pri_matrix the updated cell-type marker matrix
#' @slot current_cell_prob the current cell probability
#' (number_cells x number_cell_type)
#' @slot current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @slot starting_cell_type_assignment the initial cell type assignments
#' (number_cells x total_rounds)
#' @slot current_beta the current beta values
#' @slot unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @slot assigned_cells cells with an assigned cell type
#'
#' @export
#' @md
Celesta <- setClass("Celesta",
  slots = c(
    # STATIC FIELDS: remain untouched after initialization
    project_name = "character",
    prior_info = "data.frame",
    marker_exp_matrix = "matrix",
    original_exp = "matrix",
    cell_ID = "numeric",
    lineage_info = "data.frame",
    coords = "matrix",
    cell_prob = "matrix",
    final_cell_type_assignment = "matrix",
    nb_list = "matrix",
    total_rounds = "numeric",
    cell_nb_in_bandwidth = "ANY",
    cell_nb_dist = "ANY",
    initial_pri_matrix = "matrix",
    anchor_cell_type_assignment = "matrix",
    dist_from_nearest_assigned_cell = "matrix",
    nb_cell_type = "ANY",
    # NON-STATIC FIELDS: are updated after initialization
    marker_exp_prob = "matrix",
    current_scoring_matrix = "matrix",
    current_pri_matrix = "matrix",
    current_cell_prob = "matrix",
    current_cell_type_assignment = "matrix",
    starting_cell_type_assignment = "matrix",
    current_beta = "matrix",
    unassigned_cells = "numeric",
    assigned_cells = "numeric"
  )
)
################################################################################
################################################################################
#
# PUBLIC FUNCTIONS
#
################################################################################
################################################################################
#' CreateCelestaObject
#'
#' @description Initializes the following fields of the Celesta object:
#' * `cell_ID`
#' * `original_exp`
#' * `marker_exp_matrix`
#' * `prior_info`
#' * `lineage_info`
#' * `total_rounds`
#' * `coords`
#' * `marker_exp_prob`
#' * `nb_list`
#' * `cell_nb_in_bandwidth`
#' * `cell_nb_dist`
#' * `current_cell_type_assignment`
#' * `anchor_cell_type_assignment`
#' * `starting_cell_type_assignment`
#' * `current_scoring_matrix`
#' * `current_cell_prob`
#'
#' @param project_title *required* name of the project (used in file names)
#' @param prior_marker_info *required* user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param imaging_data_file *required* segmented imaging data.
#' The first column must contain the cell types to be inferred. The second
#' column must contain the lineage information with the following format
#' (without spaces): # _ # _ #.
#'
#' * The first number indicates round. Cell types with the same lineage level
#' are inferred at the same round. An higher number indicates higher cell-type
#' resolution. For example, immune cells -> CD3+ T cells -> CD4+ T cells.
#'
#' * The second number indicates the previous lineage cell type number for the
#' current cell type. For example, the second number for CD3+ T cell is 5
#' because it is a subtype of immune cells which have cell type number 5.
#'
#' * The third number is a number assigned to the cell type
#' (i.e. cell type number).
#'
#' The third column and beyond are columns for protein markers.
#'
#' * If a protein marker is known to be expressed for that cell type, then it
#' is denoted by a "1".
#' * If a protein marker is known to not express for a cell type, then it is
#' denoted by a "0".
#' * If the protein marker is irrelevant or uncertain to express for a cell
#' type, then it is denoted by "NA".
#'
#' @param cofactor value used to calculate the arcsinh transform on the protein
#' marker expressions
#' @param transform_type indicates a transform type for the protein marker
#' expressions (0 = no transform, 1 = arcsinh transform)
#' @param number_of_neighbors the number of cells in a single neighborhood
#' @param bandwidth the upper distance bound used when calculating
#' neighborhoods by distance
#' @return an initialized Celesta object
#' @export
#' @md
CreateCelestaObject <- function(project_title,
                                prior_marker_info,
                                imaging_data_file,
                                cofactor = 10,
                                transform_type = 1,
                                number_of_neighbors = 5,
                                bandwidth = 100) {
  celesta_obj <- Celesta(project_name = project_title)
  # Get protein marker expressions and cell IDs
  c(cell_ids, original_exp, marker_exp_matrix) %<-% GetMarkerExpMatrix(
    prior_marker_info,
    imaging_data_file,
    cofactor = 10,
    transform_type = transform_type
  )
  celesta_obj@cell_ID <- cell_ids
  celesta_obj@original_exp <- original_exp
  celesta_obj@marker_exp_matrix <- marker_exp_matrix

  # Get user-defined prior knowledge matrix and cell lineage information
  c(lineage_info, total_rounds) %<-% GetPriorInfo(prior_marker_info)
  celesta_obj@prior_info <- prior_marker_info
  celesta_obj@lineage_info <- lineage_info
  celesta_obj@total_rounds <- total_rounds

  # Get coordinates
  celesta_obj@coords <- GetCoords(imaging_data_file)

  # Convert marker expressions to marker activation probability
  celesta_obj@marker_exp_prob <- CalcMarkerActivationProbability(
    celesta_obj@marker_exp_matrix
  )

  # Get neighboring cell information
  c(nb_list, all_cell_nb_in_bandwidth, cell_nb_dist) %<-% GetNeighborInfo(
    celesta_obj@coords,
    number_of_neighbors
  )
  celesta_obj@nb_list <- nb_list
  celesta_obj@cell_nb_in_bandwidth <- all_cell_nb_in_bandwidth
  celesta_obj@cell_nb_dist <- cell_nb_dist

  # Initialize the matrices for scoring function and probability matrix
  c(current_cell_type_assignment, current_scoring_matrix, current_cell_prob) %<-%
    InitializeCellAndScoringMatrices(
      celesta_obj@lineage_info,
      celesta_obj@marker_exp_matrix, celesta_obj@prior_info
    )
  celesta_obj@current_cell_type_assignment <- current_cell_type_assignment
  celesta_obj@anchor_cell_type_assignment <- current_cell_type_assignment
  celesta_obj@starting_cell_type_assignment <- current_cell_type_assignment
  celesta_obj@current_scoring_matrix <- current_scoring_matrix
  celesta_obj@current_cell_prob <- current_cell_prob

  return(celesta_obj)
}
################################################################################
################################################################################
#' FilterCells
#'
#' @description Filters out artifact cells from the cell type assignments
#'
#' @param celesta_obj an initialized Celesta object (provided by
#' `CreateCelestaObject`)
#' @param high_marker_threshold upper bound used to filter out questionable
#' cells
#' @param low_marker_threshold lower bound used to filter out questionable
#' cells
#' @return a Celesta object with questionable cells marked with NA
#' @export
FilterCells <- function(celesta_obj,
                        high_marker_threshold = 0.9,
                        low_marker_threshold = 0.4) {
  current_cell_type_assignment <- FilterArtifactCells(
    celesta_obj@total_rounds,
    celesta_obj@marker_exp_matrix,
    celesta_obj@marker_exp_prob,
    celesta_obj@current_cell_type_assignment,
    high_marker_threshold,
    low_marker_threshold
  )
  celesta_obj@starting_cell_type_assignment <- current_cell_type_assignment
  celesta_obj@current_cell_type_assignment <- current_cell_type_assignment
  return(celesta_obj)
}
################################################################################
################################################################################
#' AssignCells
#'
#' @description Iteratively assigns cells based on spatial and protein
#' expression information.
#'
#' @param celesta_obj an initialized and filtered Celesta object (provided by
#' `FilterCells`)
#' @param max_iteration the maximum number of iterations
#' @param cell_change_threshold user defined threshold on when the iterative
#' cell-type assignment stops. The default value is 0.01, which means that if
#' the percentage of additional assigned cells is smaller than 1% of the
#' unassigned cells, then cell-type assignment will stop. The recommended range
#' is 0.01 - 0.05. Note that the higher the cell change threshold, the more
#' cells are left unassigned.
#' @param min_diff user defined threshold on how much the largest cell-type
#' probability needs to be higher than the second largest cell-type probability.
#' The default value is 0. It is recommended to not change this value.
#' @param min_probability user defined threshold on the maximum probability
#' (i.e. a cell-type probability needs to be higher than this threshold for a
#' cell to be assigned to that cell type). The default value is 0. It is
#' recommended to not set this value higher than 0.5.
#' @param high_expression_threshold_anchor the upper threshold for each cell type
#' @param low_expression_threshold_anchor the lower threshold for each cell type
#' @param high_expression_threshold_index user defined marker expression
#' probability threshold for high expression for non-anchor cells
#' @param low_expression_threshold_index user defined marker expression
#' probability threshold for low expression for non-anchor cells
#' @param progress progress object used for the Shiny app. Do not specify
#' manually.
#' @return a fully initialized Celesta object
#' @export
AssignCells <- function(celesta_obj,
                        max_iteration = 10,
                        cell_change_threshold = 0.01,
                        min_diff = 0,
                        min_probability = 0,
                        high_expression_threshold_anchor = rep(0.7,
                          length = 50
                        ),
                        low_expression_threshold_anchor = rep(0.9,
                          length = 50
                        ),
                        high_expression_threshold_index = rep(0.5,
                          length = 50
                        ),
                        low_expression_threshold_index = rep(1,
                          length = 50
                        ),
                        progress = NULL,
                        save_result = T) {
  # Cell type assignment should normally finish within 10 minutes for ~100k
  # cells and runs pretty fast for <50k cells
  for (round in 1:celesta_obj@total_rounds) {
    celesta_obj@current_cell_type_assignment[, round] <-
      celesta_obj@starting_cell_type_assignment[, round]

    current_number_of_cells_changed <- numeric()
    lineage_info <- celesta_obj@lineage_info

    initial_pri_matrix <- GetInitialPriorMatrix(
      lineage_info,
      celesta_obj@prior_info,
      round
    )
    celesta_obj@initial_pri_matrix <- initial_pri_matrix
    celesta_obj@current_pri_matrix <- initial_pri_matrix

    unassigned_cells <- FindCellsToCheck(
      celesta_obj@current_cell_type_assignment,
      celesta_obj@lineage_info,
      celesta_obj@cell_ID,
      round
    )
    number_of_cells_to_find_identity <- length(unassigned_cells)
    print(number_of_cells_to_find_identity)

    # Calculate scores using scoring function
    cell_type_num <-
      lineage_info$Cell_type_number[which(lineage_info$Round == round)]
    celesta_obj@current_scoring_matrix <- CalculateScores(
      celesta_obj@marker_exp_prob,
      celesta_obj@current_pri_matrix,
      celesta_obj@current_scoring_matrix,
      round,
      unassigned_cells,
      cell_type_num
    )

    # Initialize the cell probability with initial scores
    celesta_obj@current_cell_prob <- celesta_obj@current_scoring_matrix

    # Assign anchor cells
    old_cell_assignment <- celesta_obj@current_cell_type_assignment[, round]
    celesta_obj@current_cell_type_assignment[, round] <- AssignCellTypes(
      celesta_obj@initial_pri_matrix,
      celesta_obj@current_cell_prob,
      celesta_obj@current_cell_type_assignment,
      celesta_obj@marker_exp_prob,
      cell_type_num,
      unassigned_cells,
      round,
      high_marker_threshold = high_expression_threshold_anchor,
      low_marker_threshold = low_expression_threshold_anchor,
      min_difference = min_diff,
      min_prob = min_probability
    )

    celesta_obj@anchor_cell_type_assignment[, round] <-
      celesta_obj@current_cell_type_assignment[, round]
    cell_type_count <- CountCellType(
      celesta_obj@prior_info,
      celesta_obj@current_cell_type_assignment,
      cell_type_num,
      round
    )
    print(cell_type_count)

    if (length(which(cell_type_count[, 2] < 1)) == length(cell_type_num)) {
      print("Too few cells identified for certain cell type,
            please consider relaxing threshold.")
      return(celesta_obj)
      break
    }

    # Find cells to check
    unassigned_cells <- FindCellsToCheck(
      celesta_obj@current_cell_type_assignment,
      celesta_obj@lineage_info,
      celesta_obj@cell_ID,
      round
    )
    assigned_cells <- FindCellsWithId(
      celesta_obj@current_cell_type_assignment,
      celesta_obj@lineage_info,
      celesta_obj@cell_ID,
      round
    )

    # Calculate beta
    celesta_obj@nb_cell_type <- NeighborCellType(
      celesta_obj@nb_list,
      celesta_obj@current_cell_type_assignment,
      cell_type_num,
      round,
      unassigned_cells
    )
    celesta_obj@dist_from_nearest_assigned_cell <-
      GetDistFromNearestAssignedCells(
        celesta_obj@cell_nb_in_bandwidth,
        celesta_obj@cell_nb_dist,
        celesta_obj@current_cell_type_assignment,
        cell_type_num,
        unassigned_cells,
        assigned_cells,
        round
      )
    celesta_obj@current_beta <- CalculateBeta(
      celesta_obj@dist_from_nearest_assigned_cell,
      scale_factor = 5,
      bandwidth = 100
    )

    iteration <- 1
    current_number_of_cells_changed[iteration] <- 1

    # Iterative cell type assignment
    while (iteration < max_iteration &
      current_number_of_cells_changed[iteration] > cell_change_threshold) {
      iteration <- iteration + 1

      # Calculate cell type probabilities
      celesta_obj@current_cell_prob[, cell_type_num] <- CalculateIndexCellProb(
        celesta_obj@current_cell_prob,
        celesta_obj@current_cell_type_assignment,
        celesta_obj@current_beta,
        celesta_obj@nb_cell_type,
        celesta_obj@current_scoring_matrix,
        cell_type_num,
        unassigned_cells,
        round
      )

      # Update cell type assignment
      old_cell_assignment <- celesta_obj@current_cell_type_assignment[, round]
      celesta_obj@current_cell_type_assignment[, round] <- AssignCellTypes(
        celesta_obj@initial_pri_matrix,
        celesta_obj@current_cell_prob,
        celesta_obj@current_cell_type_assignment,
        celesta_obj@marker_exp_prob,
        cell_type_num,
        unassigned_cells,
        round,
        high_marker_threshold = high_expression_threshold_index,
        low_marker_threshold = low_expression_threshold_index,
        min_difference = min_diff,
        min_prob = min_probability
      )
      cell_type_count <- CountCellType(
        celesta_obj@prior_info,
        celesta_obj@current_cell_type_assignment,
        cell_type_num,
        round
      )
      print(cell_type_count)

      current_number_of_cells_changed[iteration] <-
        length(which((old_cell_assignment -
          celesta_obj@current_cell_type_assignment[, round]) != 0)) /
          number_of_cells_to_find_identity
      print(current_number_of_cells_changed[iteration])
      if (current_number_of_cells_changed[iteration] < cell_change_threshold) {
        break
      }

      # Find cells to check
      unassigned_cells <- FindCellsToCheck(
        celesta_obj@current_cell_type_assignment,
        celesta_obj@lineage_info,
        celesta_obj@cell_ID,
        round
      )
      assigned_cells <- FindCellsWithId(
        celesta_obj@current_cell_type_assignment,
        celesta_obj@lineage_info,
        celesta_obj@cell_ID,
        round
      )
      if (length(unassigned_cells) == 0) {
        break
      }

      # Calculate beta
      celesta_obj@nb_cell_type <- NeighborCellType(
        celesta_obj@nb_list,
        celesta_obj@current_cell_type_assignment,
        cell_type_num,
        round,
        unassigned_cells
      )
      celesta_obj@dist_from_nearest_assigned_cell <-
        GetDistFromNearestAssignedCells(
          celesta_obj@cell_nb_in_bandwidth,
          celesta_obj@cell_nb_dist,
          celesta_obj@current_cell_type_assignment,
          cell_type_num,
          unassigned_cells,
          assigned_cells,
          round
        )
      celesta_obj@current_beta <- CalculateBeta(
        celesta_obj@dist_from_nearest_assigned_cell,
        scale_factor = 5,
        bandwidth = 100
      )

      # Update prior cell-type marker matrix
      celesta_obj@current_pri_matrix <- UpdatePriorMatrix(
        celesta_obj@current_pri_matrix,
        celesta_obj@initial_pri_matrix,
        celesta_obj@current_cell_type_assignment,
        celesta_obj@marker_exp_prob,
        round,
        cell_type_num
      )
      # Update scoring function
      celesta_obj@current_scoring_matrix <- CalculateScores(
        celesta_obj@marker_exp_prob,
        celesta_obj@current_pri_matrix,
        celesta_obj@current_scoring_matrix,
        round,
        unassigned_cells,
        cell_type_num
      )
    }

    if (!is.null(progress)) {
      currValue <- progress$getValue()
      value <- 100 / celesta_obj@total_rounds
      detail <- ifelse(round == celesta_obj@total_rounds, "Assignment complete",
        paste0(
          "Round ",
          round + 1,
          "/",
          celesta_obj@total_rounds
        )
      )
      progress$set(
        value = currValue + value,
        detail = detail
      )
    }
  }
  celesta_obj@final_cell_type_assignment <- GetFinalInferredCellTypes(
    celesta_obj@project_name,
    celesta_obj@total_rounds,
    celesta_obj@current_cell_type_assignment,
    celesta_obj@anchor_cell_type_assignment,
    celesta_obj@prior_info,
    celesta_obj@lineage_info,
    celesta_obj@coords, 
    celesta_obj@original_exp,
    save_result = save_result
  )

  if (dim(celesta_obj@final_cell_type_assignment)[2] == 0) {
    print("No cells were assigned. Please adjust the CELESTA paramters.")
  }
  return(celesta_obj)
}
################################################################################
################################################################################
#' PlotCellsAnyCombination
#'
#' @description Plots the cells using x, y coordinates with their assigned cell
#' types
#'
#' @param cell_type_assignment_to_plot the final cell type assignment for each
#' cell
#' @param coords the x, y coordinates of each cell
#' @param prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param cell_number_to_use the row number of the cell types to plot from
#' `prior_info`
#' @param cell_type_colors the colors for the cell types
#' @param test_size the size of the points in the plot
#' @param save_plot whether to save the plot
#' @param output_dir the path to the directory to where the plot will be
#' outputted. This defaults to the directory containing CELESTA_functions.R.
#' Note that the directory must exist.
#' @return writes the final cell type assignment plot
#' @export
#' @md
PlotCellsAnyCombination <- function(cell_type_assignment_to_plot,
                                    coords,
                                    prior_info,
                                    cell_number_to_use,
                                    cell_type_colors = c(
                                      palette()[2:7],
                                      "white"
                                    ),
                                    test_size = 1,
                                    save_plot = TRUE,
                                    output_dir = ".") {
  # Cannot plot more than 7 cell types
  cell_types <- c("Unknown", prior_info[cell_number_to_use, 1])
  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  y_min <- min(coords[, 2])
  y_max <- max(coords[, 2])
  range <- c(min(x_min, y_min), max(x_max, y_max))

  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for (i in 1:length(cell_number_to_use)) {
    unassigned_cells <- which(
      cell_type_assignment_to_plot == cell_number_to_use[i]
    )
    cell_index[(count + 1):(count + length(unassigned_cells))] <-
      unassigned_cells
    cell_anno[(count + 1):(count + length(unassigned_cells))] <- cell_types[i]
    count <- count + length(unassigned_cells)
  }
  df_plot <- data.frame(
    x = coords[cell_index, 1],
    y = coords[cell_index, 2],
    cell_anno = cell_anno
  )
  df_plot$cell_anno <- factor(df_plot$cell_anno, levels = c(cell_types))
  color_plot <- cell_type_colors[1:length(cell_number_to_use)]

  g <- ggplot(df_plot, aes(x = x, y = y, group = cell_anno)) +
    geom_point(aes(color = cell_anno), size = test_size) +
    scale_color_manual(values = color_plot) +
    xlim(range[1], range[2]) +
    ylim(range[1], range[2]) +
    labs(main = "") +
    theme(
      aspect.ratio = 1, panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "black"),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.text = element_text(size = 12, face = "bold")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5)))

  if (save_plot) {
    ggsave(
      path = output_dir,
      filename = "plot_cell_assignment.png",
      plot = g,
      width = 12,
      height = 12,
      units = "in",
      dpi = 300
    )
  }
  return(g)
}
################################################################################
################################################################################
#' PlotExpProb
#'
#' @description Plots the expression probabilities of cells in the tissue
#'
#' @param coords the x, y coordinates of each cell
#' @param marker_exp_prob the marker expression probability for each cell
#' @param prior_marker_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param size_to_use the size of the points in the plot
#' @param width_to_use the width of the plot
#' @param height_to_use the height of the plot
#' @param save_plot whether to save the plot
#' @param output_dir the path to the directory to where the plot will be
#' outputted. This defaults to the directory containing CELESTA_functions.R.
#' Note that the directory must exist.
#' @return writes a plot of the expression probabilities for each marker
#' @export
PlotExpProb <- function(coords,
                        marker_exp_prob,
                        prior_marker_info,
                        size_to_use = 1,
                        width_to_use = 5,
                        height_to_use = 4,
                        save_plot = TRUE,
                        output_dir = ".") {
  palette <- colorRampPalette(colors = c("white", "blue4"))
  cols <- palette(6)

  markers_to_check <- as.character(
    colnames(prior_marker_info)[3:dim(prior_marker_info)[2]]
  )
  for (i in 1:length(markers_to_check)) {
    g <- PlotSingleExpProb(
      coords,
      marker_exp_prob,
      cols,
      markers_to_check[i],
      size_to_use,
      width_to_use,
      height_to_use,
      save_plot,
      output_dir
    )
  }
  return(g)
}
################################################################################
################################################################################
#
# PRIVATE FUNCTIONS
#
################################################################################
################################################################################
#' PlotSingleExpProb
#'
#' @description Plots the expression probabilities of cells in the tissue. This
#' is use soley for the Shiny app.
#'
#' @param coords the x, y coordinates of each cell
#' @param marker_exp_prob the marker expression probability for each cell
#' @param cols the color palette for the plot
#' @param marker_to_use marker to plot
#' @param size_to_use the size of the points in the plot
#' @param width_to_use the width of the plot
#' @param height_to_use the height of the plot
#' @param save_plot whether to save the plot
#' @param output_dir the path to the directory to where the plot will be
#' outputted. This defaults to the directory containing CELESTA_functions.R.
#' Note that the directory must exist.
#' @return generates a plot of the expression probabilities for a specified
#' marker
#' @export
PlotSingleExpProb <- function(coords,
                              marker_exp_prob,
                              cols = NULL,
                              marker_to_use,
                              size_to_use = 1,
                              width_to_use = 5,
                              height_to_use = 4,
                              save_plot = TRUE,
                              output_dir = ".") {
  if (is.null(cols)) {
    palette <- colorRampPalette(colors = c("white", "blue4"))
    cols <- palette(6)
  }

  marker_exp_prob_to_use <- marker_exp_prob[, which(
    colnames(marker_exp_prob) == marker_to_use
  )]
  cols_anno <- character(length = length(marker_exp_prob_to_use))
  cols_anno[which(marker_exp_prob_to_use > 0.9)] <- ">0.9"
  cols_anno[which(marker_exp_prob_to_use > 0.8 &
    marker_exp_prob_to_use <= 0.9)] <- ">0.8"
  cols_anno[which(marker_exp_prob_to_use > 0.7 &
    marker_exp_prob_to_use <= 0.8)] <- ">0.7"
  cols_anno[which(marker_exp_prob_to_use > 0.5 &
    marker_exp_prob_to_use <= 0.7)] <- ">0.5"
  cols_anno[which(marker_exp_prob_to_use <= 0.5)] <- "<=0.5"

  mca <- data.frame(
    Coords_1 = round(coords[, 1], digits = 2),
    Coords_2 = round(coords[, 2], digits = 2),
    Exp_quantile = round(marker_exp_prob_to_use, digits = 2),
    Col_anno = cols_anno
  )
  row.names(mca) <- NULL
  colnames(mca) <- c("X", "Y", "Expression", "Color_anno")
  mca$Color_anno <- factor(mca$Color_anno,
    levels = c("<=0.5", ">0.5", ">0.7", ">0.8", ">0.9")
  )

  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  y_min <- min(coords[, 2])
  y_max <- max(coords[, 2])
  range <- c(min(x_min, y_min), max(x_max, y_max))

  filename <- paste0(marker_to_use, "_exp_prob.png")
  g <- ggplot(mca, aes(x = X, y = Y, color = Color_anno)) +
    xlim(range[1], range[2]) +
    ylim(range[1], range[2]) +
    geom_point(shape = 20, size = size_to_use) +
    ggtitle(marker_to_use) +
    theme_bw() +
    scale_colour_manual(values = c(
      cols[1], cols[2], cols[3], cols[4],
      cols[6]
    )) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 10)))

  if (save_plot) {
    ggsave(
      path = output_dir,
      filename,
      plot = g,
      width = width_to_use,
      height = height_to_use,
      units = "in",
      dpi = 300
    )
  }
  return(g)
}
################################################################################
################################################################################
#' GetMarkerExpMatrix
#'
#' @description Gets the protein marker expressions and assigns each cell a
#' cell ID.
#'
#' Only protein markers specified in `prior_marker_info` are extracted from the
#' `imaging_data_file`. Cells are assigned IDs from 1 to the total number of
#' cells. If `transform_type = 1`, then an arcsinh transform is applied to the
#' protein marker expressions.
#'
#' @param prior_marker_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param imaging_data_file segmented imaging data.
#' The first column must contain the cell types to be inferred. The second
#' column must contain the lineage information with the following format
#' (without spaces): # _ # _ #.
#'
#' * The first number indicates round. Cell types with the same lineage level
#' are inferred at the same round. An higher number indicates higher cell-type
#' resolution. For example, immune cells -> CD3+ T cells -> CD4+ T cells.
#'
#' * The second number indicates the previous lineage cell type number for the
#' current cell type. For example, the second number for CD3+ T cell is 5
#' because it is a subtype of immune cells which have cell type number 5.
#'
#' * The third number is a number assigned to the cell type
#' (i.e. cell type number).
#'
#' The third column and beyond are columns for protein markers.
#'
#' * If a protein marker is known to be expressed for that cell type, then it
#' is denoted by a "1".
#' * If a protein marker is known to not express for a cell type, then it is
#' denoted by a "0".
#' * If the protein marker is irrelevant or uncertain to express for a cell
#' type, then it is denoted by "NA".
#'
#' @param cofactor used to calculate the arcsinh transform on the protein marker
#' expressions
#' @param transform_type indicates a transform type for the protein marker
#' expressions (0 = no transform, 1 = arcsinh transform)
#' @return a list with the following information:
#' \describe{
#'   \item{`cell_ids`}{the IDs of the cells}
#'   \item{`original_exp`}{the original expression matrix (containing only the
#'   protein markers specified by `prior_marker_info`)}
#'   \item{`marker_exp_matrix` or `original_exp`}{the transformed expression
#'   matrix (or original expression matrix if a transform is not specified)}
#' }
#' @export
#' @md
GetMarkerExpMatrix <- function(prior_marker_info,
                               imaging_data_file,
                               cofactor,
                               transform_type) {
  markers_to_use <- colnames(prior_marker_info)[3:dim(prior_marker_info)[2]]
  matching_markers <- match(markers_to_use, colnames(imaging_data_file))
  if (length(which(is.na(matching_markers) == TRUE)) > 0) {
    stop("Please double check that the protein markers in the user-defined
         cell-type signature matrix (prior_marker_info) are included in the
         protein markers in the segmented imaging input file
         (imaging_data_file).")
  }

  original_exp <- data.matrix(imaging_data_file[, matching_markers])
  cell_ids <- seq(1, dim(original_exp)[1], by = 1)

  if (transform_type == 1) { # arcsinh transformation
    marker_exp_matrix <- asinh(original_exp / cofactor)
    return(list(cell_ids, original_exp, marker_exp_matrix))
  }

  return(list(cell_ids, original_exp, original_exp))
}
################################################################################
################################################################################
#' GetPriorInfo
#'
#' @description Extracts the lineage information from  the `prior_marker_info`
#' and determines the total rounds
#'
#' @param prior_marker_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @return a list with the following information:
#' \describe{
#'   \item{`lineage_info`}{the lineage information parsed into round, previous
#'   cell type, and cell type number columns}
#'   \item{`total_rounds`}{the maximum round value}
#' }
#' @export
#' @md
GetPriorInfo <- function(prior_marker_info) {
  if (FALSE %in% grepl("_", prior_marker_info[, 2], fixed = TRUE)) {
    stop("Warning: the lineage information column has formatting errors")
  }

  c(round, previous_cell_type, cell_type_number) %<-%
    list(integer(), integer(), integer())
  for (i in 1:dim(prior_marker_info)[1]) {
    info <- strtoi(unlist(strsplit(prior_marker_info[i, 2], "_")))
    c(round[i], previous_cell_type[i], cell_type_number[i]) %<-%
      list(info[1], info[2], info[3])
  }
  lineage_info <- data.frame(
    Round = round,
    Previous_cell_type = previous_cell_type,
    Cell_type_number = cell_type_number
  )
  total_rounds <- max(lineage_info$Round)
  return(list(lineage_info, total_rounds))
}
################################################################################
################################################################################
#' GetCoords
#'
#' @description Gets the x, y coordinates of each cell
#'
#' @param imaging_data_file segmented imaging data.
#' The first column must contain the cell types to be inferred. The second
#' column must contain the lineage information with the following format
#' (without spaces): # _ # _ #.
#'
#' * The first number indicates round. Cell types with the same lineage level
#' are inferred at the same round. An higher number indicates higher cell-type
#' resolution. For example, immune cells -> CD3+ T cells -> CD4+ T cells.
#'
#' * The second number indicates the previous lineage cell type number for the
#' current cell type. For example, the second number for CD3+ T cell is 5
#' because it is a subtype of immune cells which have cell type number 5.
#'
#' * The third number is a number assigned to the cell type
#' (i.e. cell type number).
#'
#' The third column and beyond are columns for protein markers.
#'
#' * If a protein marker is known to be expressed for that cell type, then it
#' is denoted by a "1".
#' * If a protein marker is known to not express for a cell type, then it is
#' denoted by a "0".
#' * If the protein marker is irrelevant or uncertain to express for a cell
#' type, then it is denoted by "NA".
#'
#' @returns the x, y coordinates of each cell
#' @export
#' @md
GetCoords <- function(imaging_data_file) {
  coords <- cbind(
    imaging_data_file$X,
    imaging_data_file$Y
  )
  colnames(coords) <- c("X", "Y")
  return(coords)
}
################################################################################
################################################################################
#' FitGmmModel
#'
#' @description Fits a Gaussian mixture model for each marker
#'
#' @param marker_exp the expression of the marker for each cell
#' @param marker_name the name of the marker
#' @param figure whether a figure should be generated or not
#' @return the Gaussian mixture model parameters for the marker
#' @export
#' @md
FitGmmModel <- function(marker_exp, marker_name, figure = FALSE) {
  cat("Marker: ", marker_name, "\n")

  gmm_marker_param <- matrix(nrow = 3, ncol = 2)

  set.seed(1)
  zero_indices <- which(marker_exp == 0)
  zero_percentage <- length(zero_indices) / length(marker_exp)

  if (zero_percentage > 0.1 & zero_percentage < 0.2) {
    print("Warning: The marker expression potentially has too many zeros for
          fitting. GMM fitting will use input expression data with reduced
          sparsity")
    num_of_indices_to_remove <- floor(length(marker_exp) * (zero_percentage))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp, 2,
      models = mixmodGaussianModel(
        family = "general",
        listModels = "Gaussian_p_Lk_Ck",
        free.proportions = FALSE, equal.proportions = TRUE
      )
    )
  } else if (zero_percentage >= 0.2 & zero_percentage < 0.5) {
    print("Warning: The marker expression potentially has too many zeros for
          fitting. GMM fitting will use input expression data with reduced
          sparsity")
    num_of_indices_to_remove <- floor(length(marker_exp) * (zero_percentage -
      0.05))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp, 2,
      models = mixmodGaussianModel(
        family = "general",
        listModels = "Gaussian_p_Lk_Ck",
        free.proportions = FALSE, equal.proportions = TRUE
      )
    )
  } else if (zero_percentage >= 0.5) {
    print("Warning: The marker expression potentially has too many zeros for
          fitting. GMM fitting will use input expression data with reduced
          sparsity")
    num_of_indices_to_remove <- ceiling(length(marker_exp) * (zero_percentage -
      0.02))
    marker_exp <- marker_exp[-zero_indices[1:num_of_indices_to_remove]]
    xxx <- mixmodCluster(marker_exp, 2,
      models = mixmodGaussianModel(
        family = "general",
        listModels = "Gaussian_p_Lk_Ck",
        free.proportions = FALSE, equal.proportions = TRUE
      )
    )
  } else {
    xxx <- mixmodCluster(marker_exp, 2,
      models = mixmodGaussianModel(
        family = "general",
        listModels = "Gaussian_p_Lk_Ck",
        free.proportions = FALSE, equal.proportions = TRUE
      )
    )
  }

  # Check the models information for the Gaussian models, which shows which
  # parameters are constrained
  # Want equal proportions of the two Gaussians
  gmm_marker_param[1, ] <- xxx@results[[1]]@parameters@proportions
  gmm_marker_param[2, ] <- xxx@results[[1]]@parameters@mean[, 1]
  gmm_marker_param[3, 1] <- xxx@results[[1]]@parameters@variance[[1]][, 1]
  gmm_marker_param[3, 2] <- xxx@results[[1]]@parameters@variance[[2]][, 1]


  if (figure == TRUE) {
    bin_size <- 20
    filename <- paste0(marker_name, "_GMM.png")
    png(filename, width = 5.5, height = 6.5, units = "in", res = 300)
    h <- hist(marker_exp,
      breaks = bin_size, xlab = "Marker expression", main =
        paste0("Histogram for ", marker_name)
    )
    highestCount <- max(h$counts)
    multiplier <- h$counts / h$density

    xfit <- seq(min(marker_exp), max(marker_exp), length = length(h$breaks))
    yfit1 <- dnorm(xfit,
      mean = gmm_marker_param[2, 1],
      sd = sqrt(gmm_marker_param[3, 1])
    ) * multiplier[1]
    lines(xfit, yfit1, col = "blue", lwd = 2)
    yfit2 <- dnorm(xfit,
      mean = gmm_marker_param[2, 2],
      sd = sqrt(gmm_marker_param[3, 2])
    ) * multiplier[1]
    lines(xfit, yfit2, col = "red", lwd = 2)
    dev.off()
  }
  return(gmm_marker_param)
}
################################################################################
################################################################################
#' BuildSigmoidFunction
#'
#' @description Builds the sigmoid function for the calculation of the
#' expression probability
#'
#' @param marker_exp_matrix transformed protein marker expression (or original
#' segmentation protein marker expression if transformation is not specified)
#' @param figure whether a figure should be generated or not
#' @return the sigmoid function parameter, containing the \eqn{x_root} and slope
#' @export
#' @md
BuildSigmoidFunction <- function(marker_exp_matrix, figure = FALSE) {
  sigmoid_function_parameter <- matrix(
    nrow = 2,
    ncol = dim(marker_exp_matrix)[2]
  )
  # For each marker, fit GMM
  for (i in 1:dim(marker_exp_matrix)[2]) {
    marker_exp <- marker_exp_matrix[, i]
    marker_name <- colnames(marker_exp_matrix)[i]
    if (typeof(marker_name) != "character") {
      stop("Protein marker name in the marker expression matrix has potential
           problem.")
    }

    marker_GMM_model <- FitGmmModel(marker_exp, marker_name, figure)
    c(weight, mus, sigmas) %<-% list(
      marker_GMM_model[1, ],
      marker_GMM_model[2, ],
      marker_GMM_model[3, ]
    )

    if (mus[1] > mus[2]) {
      # The first Gaussian model is for marker expressed,
      # second is for marker not expressed
      a <- (-0.5 / sigmas[2] + 0.5 / sigmas[1])
      b <- mus[2] / sigmas[2] - mus[1] / sigmas[1]
      c <- 0.5 * (-mus[2]^2 / sigmas[2] + mus[1]^2 / sigmas[1]) +
        log(weight[2] / weight[1]) + 0.5 * log(sigmas[1] / sigmas[2])
    } else {
      # The second Gaussian model is for marker expressed,
      # first is for marker not expressed
      a <- (-0.5 / sigmas[1] + 0.5 / sigmas[2])
      b <- mus[1] / sigmas[1] - mus[2] / sigmas[2]
      c <- 0.5 * (-mus[1]^2 / sigmas[1] + mus[2]^2 / sigmas[2]) +
        log(weight[1] / weight[2]) + 0.5 * log(sigmas[2] / sigmas[1])
    }
    xroot <- (-b - sqrt(b^2 - 4.0 * a * c)) / (2.0 * a)
    slope <- 1

    if (figure == TRUE) {
      filename <- paste0(marker_name, "_sigmoid.png")

      exp_term <- exp(slope * (marker_exp - xroot))
      yyy <- exp_term / (1 + exp_term)
      yyy <- (yyy - min(yyy)) / (max(yyy) - min(yyy))

      # Plot sigmoid function
      png(filename, width = 4.5, height = 4.5, units = "in", res = 300)
      plot(marker_exp, yyy,
        col = "darkblue",
        xlab = "", ylab = "", main = paste0(marker_name, " sigmoid function")
      )
      grid()
      dev.off()
    }
    sigmoid_function_parameter[1, i] <- xroot
    sigmoid_function_parameter[2, i] <- slope
  }

  return(sigmoid_function_parameter)
}
################################################################################
################################################################################
#' CalcMarkerActivationProbability
#'
#' @description Calculates the activation probability for each marker in the
#' prior matrix
#'
#' @param marker_exp_matrix transformed protein marker expression (or original
#' segmentation
#' protein marker expression if transformation is not specified)
#' @return the protein marker activation probability
#' @export
CalcMarkerActivationProbability <- function(marker_exp_matrix, figure = FALSE) {
  # Fit GMM model and get parameters for the activation probabilities
  sigmoid_function_parameter <- BuildSigmoidFunction(marker_exp_matrix, figure)

  # Marker activation probability matrix
  marker_exp_prob <- matrix(
    nrow = dim(marker_exp_matrix)[1],
    ncol = dim(marker_exp_matrix)[2]
  )
  colnames(marker_exp_prob) <- colnames(marker_exp_matrix)

  for (i in 1:dim(marker_exp_matrix)[2]) {
    exp_term <- exp(sigmoid_function_parameter[2, i] * (marker_exp_matrix[, i] -
      sigmoid_function_parameter[1, i]))
    y <- exp_term / (1 + exp_term)
    marker_exp_prob[, i] <- (y - min(y)) / (max(y) - min(y))
  }
  return(marker_exp_prob)
}
################################################################################
################################################################################
#' GetNeighborInfo
#'
#' @description Gets the neighborhood information, including neighborhoods by
#' number and distance.
#'
#' @param coords the x, y coordinates of each cell
#' @param number_of_neighbors the number of cells in a single neighborhood
#' @param bandwidth the upper distance bound used when calculating neighborhoods
#' by distance
#' @returns a list of the following information
#' \describe{
#'   \item{`nb_list`}{the list of N-nearest neighbors}
#'   \item{`all_cell_nb_in_bandwidth`}{the cells located within a bandwidth to
#'   cell *c*}
#'   \item{`cell_nb_dist`}{the distance of each cell to cell *c* within a
#'   bandwidth}
#' }
#' @export
GetNeighborInfo <- function(coords, number_of_neighbors = 5, bandwidth = 100) {
  print("Getting the nearest neighbors")
  nb_list <- knearneigh(coords, k = number_of_neighbors)$nn
  colnames(nb_list) <- paste0("neighbor", seq(1, number_of_neighbors, by = 1))

  print("Identifying neighboring cells within a defined circle bandwidth")
  all_cell_nb_in_bandwidth <- dnearneigh(coords, 0, bandwidth, longlat = NULL)

  print("Identify distances for all the cells within the circle bandwidth")
  cell_nb_dist <- nbdists(all_cell_nb_in_bandwidth, coords)
  return(list(nb_list, all_cell_nb_in_bandwidth, cell_nb_dist))
}
################################################################################
################################################################################
#' InitializeCellAndScoringMatrices
#'
#' @description Initialize the cell type assignments, cell probabilities, and
#' scoring matrices
#'
#' @param lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @param marker_exp_matrix transformed protein marker expression (or original
#' segmentation protein marker expression if transformation is not specified)
#' @param prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @return a list with the following information
#' \describe{
#'   \item{`current_cell_type_assignment`}{a zero matrix with dimension
#'   (number_cells x total_rounds)}
#'   \item{`current_scoring_matrix`}{a NA matrix with dimension
#'   (number_cells x number_cell_type)}
#'   \item{`current_cell_prob`}{a NA matrix with dimension
#'   (number_cells x number_cell_type)}
#' }
#' @export
InitializeCellAndScoringMatrices <- function(lineage_info,
                                             marker_exp_matrix,
                                             prior_info) {
  total_rounds <- max(lineage_info$Round)
  current_cell_type_assignment <- matrix(
    0L,
    nrow = dim(marker_exp_matrix)[1],
    ncol = total_rounds
  )

  current_scoring_matrix <- matrix(
    nrow = dim(marker_exp_matrix)[1],
    ncol = dim(prior_info)[1]
  )
  colnames(current_scoring_matrix) <- prior_info[, 1]

  current_cell_prob <- matrix(
    nrow = dim(marker_exp_matrix)[1],
    ncol = dim(prior_info)[1]
  )
  colnames(current_cell_prob) <- prior_info[, 1]
  return(list(
    current_cell_type_assignment,
    current_scoring_matrix,
    current_cell_prob
  ))
}
################################################################################
################################################################################
#' FilterArtifactCells
#'
#' @description Filter out cells that could potentially be artifacts
#'
#' @param total_rounds the maximum round value
#' @param marker_exp_matrix transformed protein marker expression (or original
#' segmentation protein marker expression if transformation is not specified)
#' @param marker_exp_prob the marker expression probability for each cell
#' @param current_cell_type_assignment the cell type assignments for each round
#' for each cell
#' @param high_marker_threshold upper bound used to filter out questionable
#' cells
#' @param low_marker_threshold lower bound used to filter out questionable
#' cells
#' @return current cell type assignment, where a questionable cells are marked
#' with a row of NAs.
#' @export
FilterArtifactCells <- function(total_rounds,
                                marker_exp_matrix,
                                marker_exp_prob,
                                current_cell_type_assignment,
                                high_marker_threshold = 0.9,
                                low_marker_threshold = 0.4) {
  # Filter out cells that are questionable
  for (i in 1:dim(marker_exp_matrix)[1]) {
    cell_activation_prob <- marker_exp_prob[i, ]
    if (MarkQuestionableCells(
      cell_activation_prob,
      high_marker_threshold,
      low_marker_threshold
    )) {
      current_cell_type_assignment[i, 1:total_rounds] <- rep(NA, total_rounds)
    }
  }
  return(current_cell_type_assignment)
}
################################################################################
################################################################################
#' MarkQuestionableCells
#'
#' @description Determine if a cell is questionable.
#'
#' A cell is questionable if *all* of its protein marker expressions are below
#' the `lower_marker_threshold` or above the `high_marker_threshold`.
#'
#' @param cell_activation_prob the protein marker expressions for a single cell
#' @param high_marker_threshold upper bound used to filter out questionable
#' cells
#' @param low_marker_threshold lower bound used to filter out questionable
#' cells
#' @return whether a cell is questionable or not
#' @export
MarkQuestionableCells <- function(cell_activation_prob,
                                  high_marker_threshold,
                                  low_marker_threshold) {
  number_of_marker <- length(cell_activation_prob)
  number_of_low_markers <- length(which(cell_activation_prob <
    low_marker_threshold))
  number_of_high_markers <- length(which(cell_activation_prob >
    high_marker_threshold))
  if (number_of_low_markers == number_of_marker | number_of_high_markers ==
    number_of_marker) {
    return(TRUE)
  }
  return(FALSE)
}
################################################################################
################################################################################
#' GetInitialPriorMatrix
#'
#' @description Gets the prior knowledge of the cell types with the specified
#' round.
#'
#' @param lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @param round the current round
#' @return the prior knowledge of the cells types with the specified round.
#' @export
GetInitialPriorMatrix <- function(lineage_info, prior_marker_info, round) {
  initial_pri_matrix <- data.matrix(prior_marker_info[
    which(lineage_info$Round == round),
    3:dim(prior_marker_info)[2]
  ])
  return(initial_pri_matrix)
}
################################################################################
################################################################################
#' FindCellsToCheck
#'
#' @description Find unassigned cells
#'
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @param cell_ID the IDs of the cells (from 1 to the total number of cells)
#' @return the IDs of unassigned cells
#' @export
FindCellsToCheck <- function(current_cell_type_assignment,
                             lineage_info,
                             cell_ID,
                             round) {
  if (round == 1) {
    unassigned_cells <- cell_ID[which(
      current_cell_type_assignment[, round] == 0
    )]
  } else {
    previous_level_type <- unique(lineage_info$Previous_cell_type[which(
      lineage_info$Round == round
    )])
    previous_level_round <- lineage_info$Round[which(
      lineage_info$Cell_type_number == previous_level_type
    )]
    unassigned_cells <- cell_ID[which(
      current_cell_type_assignment[, round] == 0 &
        (current_cell_type_assignment[, previous_level_round] ==
          previous_level_type)
    )]
  }

  return(unassigned_cells)
}
################################################################################
################################################################################
#' FindCellsWithId
#'
#' @description Find cells assigned with a cell type
#'
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @param cell_ID the IDs of the cells (from 1 to the total number of cells)
#' @param round the current round
#' @return cells that have been assigned a cell type
#' @export
FindCellsWithId <- function(current_cell_type_assignment,
                            lineage_info,
                            cell_ID,
                            round) {
  if (round == 1) {
    assigned_cells <- cell_ID[which(current_cell_type_assignment[, round] != 0 &
      is.na(current_cell_type_assignment[, round]) == FALSE)]
  } else {
    previous_level_type <- unique(lineage_info$Previous_cell_type[which(
      lineage_info$Round == round
    )])
    assigned_cells <- cell_ID[which(current_cell_type_assignment[, round] != 0 &
      is.na(current_cell_type_assignment[, round]) == FALSE &
      (current_cell_type_assignment[, (round - 1)] == previous_level_type))]
  }

  return(assigned_cells)
}
################################################################################
################################################################################
#' GetScore
#'
#' @description Calculate scores using MSE
#'
#' @param activation_prob_to_use the marker expression probabilities of the
#' unassigned cells
#' @param prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param non_NA_index the index of the columns in `current_pri_matrix` that do
#' not contain NA for a particular cell
#' @return the score of the cell
#' @export
GetScore <- function(activation_prob_to_use, prior_info, non_NA_index) {
  score <- apply(
    activation_prob_to_use[, non_NA_index], 1,
    function(x) (1 - sum((x - prior_info)^2) / length(x))
  )
  return(score)
}
################################################################################
################################################################################
#' CalculateScores
#'
#' @description Calculate the scores based on the scoring function
#'
#' @param marker_exp_prob the marker expression probability for each cell
#' @param current_pri_matrix the updated cell-type marker matrix
#' @param current_scoring_matrix the current scoring matrix
#' (number_cells x number_cell_type)
#' @param round the current round
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @param cell_type_num the cell types associated with the current round
#' @return the current scoring matrix containing the scores for each cell type
#' associated with the current round for each unassigned cell
#' @export
CalculateScores <- function(marker_exp_prob,
                            current_pri_matrix,
                            current_scoring_matrix,
                            round,
                            unassigned_cells,
                            cell_type_num) {
  print("Start calculating the scoring function.")
  activation_prob_to_use <- marker_exp_prob[unassigned_cells, ]

  for (i in 1:length(cell_type_num)) {
    non_NA_index <- which(!is.na(current_pri_matrix[i, ]))
    prior_info <- current_pri_matrix[i, non_NA_index]
    current_scoring_matrix[unassigned_cells, cell_type_num[i]] <- GetScore(
      activation_prob_to_use, prior_info, non_NA_index
    )
  }

  current_scoring_matrix[unassigned_cells, cell_type_num] <- t(apply(
    current_scoring_matrix[unassigned_cells, cell_type_num],
    1, function(x) x / sum(x)
  ))

  return(current_scoring_matrix)
}
################################################################################
################################################################################
#' CalculateProbabilityDifference
#'
#' @description Calculate the probability differences
#'
#' @param max.prob the maximum marker probability for each cell
#' @param max.prob_index the index of the maximum marker probability for each
#' cell
#' @param cell_prob_list the probabilities of the cells are are not assigned a
#' cell type
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @return the minimum of the difference in probability between the maximum
#' marker probability and other marker probabilities
#' @export
CalculateProbabilityDifference <- function(max.prob,
                                           max.prob_index,
                                           cell_prob_list,
                                           unassigned_cells) {
  # max.prob, max.prob_index are calculated only on unassigned_cells
  # but cell_prob_list has all the cells
  min_prob_diff <- numeric(length = length(unassigned_cells))
  for (i in 1:length(unassigned_cells)) {
    min_prob_diff[i] <- min(max.prob[i] - cell_prob_list[
      unassigned_cells[i],
      -max.prob_index[i]
    ])
  }
  return(min_prob_diff)
}
################################################################################
################################################################################
#' AssignCellTypes
#'
#' @description Find the cell types based on the scores (anchor cell) or
#' probabilities (index cell)
#'
#' @param initial_pri_matrix user defined cell-type marker matrix for a
#' specific round
#' @param current_cell_prob the current cell probability
#' (number_cells x number_cell_type)
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param marker_exp_prob the marker expression probability for each cell
#' @param cell_type_num the cell types associated with the current round
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @param round the current round
#' @param high_marker_threshold the upper threshold for each cell type
#' @param low_marker_threshold the lower threshold for each cell type
#' @param min_difference lower bound used to determine cells that meet the
#' threshold
#' @param min_prob lower bound used to determine cells that meet the threshold
#' @return an updated current cell type assignment (number_cells x total_rounds)
#' with more cells assigned for the current round
#' @export
AssignCellTypes <- function(initial_pri_matrix,
                            current_cell_prob,
                            current_cell_type_assignment,
                            marker_exp_prob,
                            cell_type_num,
                            unassigned_cells,
                            round,
                            high_marker_threshold,
                            low_marker_threshold,
                            min_difference = 0,
                            min_prob = 0) {
  cell_prob_list <- current_cell_prob[, cell_type_num]
  cell_type_assignment <- current_cell_type_assignment[, round]

  max.prob_index <- apply(cell_prob_list[unassigned_cells, ], 1, which.max)
  max.prob <- apply(cell_prob_list[unassigned_cells, ], 1, max)
  min_prob_diff <- CalculateProbabilityDifference(
    max.prob,
    max.prob_index,
    cell_prob_list,
    unassigned_cells
  )

  # Find cells with cell type max probability > threshold
  # and cell type probability difference > threshold
  # Indexing on unassigned_cells
  threshold_cells <- unassigned_cells[which(min_prob_diff > min_difference &
    max.prob > min_prob)]
  max.prob_index_thresholded <- max.prob_index[which(
    min_prob_diff > min_difference &
      max.prob > min_prob
  )]

  for (i in 1:length(threshold_cells)) {
    cell_ID_to_check <- threshold_cells[i]
    high_marker_index <- which(
      initial_pri_matrix[max.prob_index_thresholded[i], ] == 1
    )
    low_marker_index <- which(
      initial_pri_matrix[max.prob_index_thresholded[i], ] == 0
    )
    threshold_index <- cell_type_num[max.prob_index_thresholded[i]]

    if (
      length(which(marker_exp_prob[cell_ID_to_check, high_marker_index] >=
        high_marker_threshold[threshold_index])) ==
        length(high_marker_index) &
        length(which(marker_exp_prob[cell_ID_to_check, low_marker_index] <=
          low_marker_threshold[threshold_index])) ==
          length(low_marker_index)) {
      cell_type_assignment[cell_ID_to_check] <-
        cell_type_num[max.prob_index_thresholded[i]]
    }
  }
  return(cell_type_assignment)
}
################################################################################
################################################################################
#' CountCellType
#'
#' @description Counts the cell type
#'
#' @param prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param cell_type_num the cell types associated with the current round
#' @param round the current round
#' @return the count and proportion for each cell type based on the current cell
#' type assignments
#' @export
CountCellType <- function(prior_info,
                          current_cell_type_assignment,
                          cell_type_num,
                          round) {
  cell_type_count <- matrix(nrow = (length(cell_type_num)), ncol = 3)
  colnames(cell_type_count) <- c("cell_type_number", "count", "proportion")
  row.names(cell_type_count) <- prior_info[cell_type_num, 1]
  cell_type_count[, 1] <- cell_type_num
  total_cell_number <- dim(current_cell_type_assignment)[1]

  for (i in 1:length(cell_type_num)) {
    cell_type_count[i, 2] <- length(which(
      current_cell_type_assignment[, round] == cell_type_num[i]
    ))
    cell_type <- prior_info[cell_type_num[i], 1]
    if (cell_type_count[i, 2] < 1) {
      print(paste0("Too few cells identified for: ", cell_type))
      print("Please consider relaxing the threshold.")
    }
  }

  cell_type_count[, 3] <- cell_type_count[, 2] / total_cell_number
  return(cell_type_count)
}
################################################################################
################################################################################
#' NeighborCellType
#'
#' @description Find the cell types of the neighbors of unassigned cells
#'
#' @param nb_list the list of N-nearest neighbors
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param cell_type_num the cell types associated with the current round
#' @param round the current round
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @return the cell types of the neighbors of unassigned cells
#' @export
NeighborCellType <- function(nb_list,
                             current_cell_type_assignment,
                             cell_type_num,
                             round,
                             unassigned_cells) {
  cell_type_assignment <- current_cell_type_assignment[, round]
  same_type_nb <- matrix(rep(
    list(),
    length(cell_type_num) * length(unassigned_cells)
  ),
  nrow = length(unassigned_cells), ncol = length(cell_type_num)
  )
  row.names(same_type_nb) <- unassigned_cells
  colnames(same_type_nb) <- cell_type_num
  for (j in 1:length(unassigned_cells)) {
    current_cell_ID <- unassigned_cells[j]
    neighbors <- nb_list[current_cell_ID, ]
    neighbor_types <- cell_type_assignment[neighbors]
    for (i in 1:length(cell_type_num)) {
      same_type_nb[j, i][[1]] <-
        neighbors[which(neighbor_types == cell_type_num[i])]
    }
  }
  return(same_type_nb)
}
################################################################################
################################################################################
#' GetDistFromNearestAssignedCells
#'
#' @description Get distance from nearest assigned cells
#'
#' @param cell_nb_in_bandwidth the cells located within a bandwidth to cell *c*
#' @param cell_nb_dist the distance of each cell to cell *c* within a bandwidth
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param cell_type_num the cell types associated with the current round
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @param assigned_cells cells with an assigned cell type
#' @param round the current round
#' @return the distance to the nearest assigned cells
#' @export
GetDistFromNearestAssignedCells <- function(cell_nb_in_bandwidth,
                                            cell_nb_dist,
                                            current_cell_type_assignment,
                                            cell_type_num,
                                            unassigned_cells,
                                            assigned_cells,
                                            round) {
  print("Get distance from nearest assigned cells.")
  dist_nearest_assigned_cell <- matrix(
    nrow = length(unassigned_cells),
    ncol = length(cell_type_num)
  )
  colnames(dist_nearest_assigned_cell) <- cell_type_num

  for (i in 1:dim(dist_nearest_assigned_cell)[1]) {
    cell_to_check <- unassigned_cells[i]
    matching <- match(cell_nb_in_bandwidth[[cell_to_check]], assigned_cells)
    index <- matching[which(is.na(matching) == FALSE)]

    if (length(index) == 0) {
      next
    }

    nb_cell_with_ID <- assigned_cells[index]
    nb_cell_type <- current_cell_type_assignment[nb_cell_with_ID, round]
    unique_nb_cell_type <- unique(nb_cell_type)
    nb_cell_dist <- cell_nb_dist[[cell_to_check]][which(is.na(matching) ==
      FALSE)]

    for (j in 1:length(unique_nb_cell_type)) {
      type_j <- which(nb_cell_type == unique_nb_cell_type[j])
      dist_nearest_assigned_cell[i, which(cell_type_num ==
        unique_nb_cell_type[j])] <- min(nb_cell_dist[type_j])
    }
  }

  return(dist_nearest_assigned_cell)
}
################################################################################
################################################################################
#' CalculateBeta
#'
#' @description Calculates beta
#'
#' @param dist_from_nearest_assigned_cell the distance from the nearest assigned
#' cell
#' @param scale_factor the scale factor
#' @param bandwidth the bandwidth
#' @return the beta value
#' @export
#' @md
CalculateBeta <- function(dist_from_nearest_assigned_cell,
                          scale_factor = 5,
                          bandwidth = 100) {
  beta <- scale_factor * (1 - dist_from_nearest_assigned_cell / bandwidth)
  beta[is.na(beta)] <- 0
  return(beta)
}
################################################################################
################################################################################
#' CalculateIndexCellProb
#'
#' @description Calculates the probability for index cells
#'
#' @param current_cell_prob the current cell probability
#' (number_cells x number_cell_type)
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param current_beta the current beta values
#' @param nb_cell_type cell types of the neighboring cells for index cells
#' @param current_scoring_matrix the current scoring matrix
#' (number_cells x number_cell_type)
#' @param cell_type_num the cell types associated with the current round
#' @param unassigned_cells cells not assigned a cell type for each round and
#' iteration
#' @param round the current round
#' @return calculates the probability for each cell type for unassigned cells
#' @export
#' @md
CalculateIndexCellProb <- function(current_cell_prob,
                                   current_cell_type_assignment,
                                   current_beta,
                                   nb_cell_type,
                                   current_scoring_matrix,
                                   cell_type_num,
                                   unassigned_cells,
                                   round) {
  # This function uses mean field estimation to calculate the probability
  # For each cell, a probability is calculated for each cell type to check
  current_cell_prob_list <- current_cell_prob[, cell_type_num]

  # all cells * cell_type_num
  u <- current_cell_prob_list
  # all cells
  current_cell_type_assignment <- current_cell_type_assignment[, round]

  for (i in 1:length(unassigned_cells)) {
    cell_ID_to_check <- unassigned_cells[i]
    u_i <- numeric(length = length(cell_type_num))
    number_of_nb <- lengths(nb_cell_type[i, ])

    for (j in 1:length(number_of_nb)) {
      current_same_type_nb <- unlist(nb_cell_type[i, j][[1]])
      u_i[j] <- exp(
        current_scoring_matrix[cell_ID_to_check, cell_type_num[j]]
      ) *
        exp(current_beta[i, j] *
          sum(current_cell_prob_list[current_same_type_nb, j]))
    }
    u[cell_ID_to_check, ] <- u_i / sum(u_i)
  }
  print("Updating cell probability done.")
  return(u)
}
################################################################################
################################################################################
#' UpdatePriorMatrix
#'
#' @description Updates prior knowledge matrix of the cell type signatures
#'
#' @param current_pri_matrix the updated cell-type marker matrix
#' @param initial_pri_matrix user defined cell-type marker matrix for a specific
#' round
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param marker_exp_prob the marker expression probability for each cell
#' @param round the current round
#' @param cell_type_num the cell types associated with the current round
#' @return updates the prior knowledge matrix with information from cells
#' assigned to each cell type
#' @export
UpdatePriorMatrix <- function(current_pri_matrix,
                              initial_pri_matrix,
                              current_cell_type_assignment,
                              marker_exp_prob,
                              round,
                              cell_type_num) {
  updated_prior_matrix <- current_pri_matrix
  for (i in 1:length(cell_type_num)) {
    cell_type_to_check <- cell_type_num[i]
    for (j in 1:dim(current_pri_matrix)[2]) {
      cells_of_current_cell_type <- which(
        current_cell_type_assignment[, round] == cell_type_to_check
      )
      updated_prior_matrix[i, j] <- (mean(
        marker_exp_prob[cells_of_current_cell_type, j]
      ) +
        initial_pri_matrix[i, j]) / 2
    }
  }

  return(updated_prior_matrix)
}
################################################################################
################################################################################
#' GetFinalInferredCellTypes
#'
#' @description Get final cell types and writes two files: the final cell type
#' assignments and the anchor cell type assignments.
#'
#' @param total_rounds the maximum round
#' @param current_cell_type_assignment the current cell type assignments
#' (number_cells x total_rounds)
#' @param anchor_cell_type_assignment the anchor cell type assignments
#' @param prior_info user-defined cell-type signature matrix.
#'
#' The data should contain two columns (name X and Y) for the x, y coordinates
#' and a column for each protein marker. Each row represents the data for a
#' single cell, including its x, y coordinates and expression for each protein
#' marker.
#'
#' @param lineage_info the lineage information from `prior_info` parsed into
#' round, previous cell type, and cell type number columns
#' @param coords the x, y coordinates of each cell
#' @param original_exp original protein marker expression (containing only the
#' protein markers specified in `prior_info`)
#' @param save_data whether or not to save the final cell type assignment
#' and anchor cell assignment results
#' @return the final cell type assignments
#' @export
#' @md
GetFinalInferredCellTypes <- function(project_name,
                                      total_rounds,
                                      current_cell_type_assignment,
                                      anchor_cell_type_assignment,
                                      prior_info,
                                      lineage_info,
                                      coords,
                                      original_exp,# TODO: coords & original exp matrix
                                      save_result = T) {
  cell_type_name_assigned <- matrix(
    nrow = dim(current_cell_type_assignment),
    ncol = total_rounds
  )
  anchor_cell_type_name_assigned <- matrix(
    nrow = dim(current_cell_type_assignment),
    ncol = total_rounds
  )
  final_cell_type_assignment <- rep(0,
    length = dim(current_cell_type_assignment)[1]
  )
  for (i in 1:total_rounds) {
    current_pri_matrix_num <- i
    cell_type_name_assigned[, i] <- prior_info[match(
      current_cell_type_assignment[, i],
      lineage_info$Cell_type_number
    ), 1]
    cell_type_name_assigned[which(current_cell_type_assignment[, i] ==
      0), i] <- "Unknown"
    anchor_cell_type_name_assigned[, i] <- prior_info[match(
      anchor_cell_type_assignment[, i],
      lineage_info$Cell_type_number
    ), 1]
    anchor_cell_type_name_assigned[which(anchor_cell_type_assignment[, i] ==
      0), i] <- "Unknown"
    if (current_pri_matrix_num == 1) {
      final_cell_type_assignment <-
        current_cell_type_assignment[, current_pri_matrix_num]
    } else {
      previous_level_type <- unique(
        lineage_info$Previous_cell_type[which(lineage_info$Round ==
          current_pri_matrix_num)]
      )
      assignment <- current_cell_type_assignment[
        which(
          final_cell_type_assignment == previous_level_type &
            current_cell_type_assignment[, current_pri_matrix_num] != 0
        ),
        current_pri_matrix_num
      ]
      final_cell_type_assignment[which(final_cell_type_assignment ==
        previous_level_type &
        current_cell_type_assignment[, current_pri_matrix_num] != 0)] <-
        assignment
    }
  }
  final_cell_names <- character(length = dim(current_cell_type_assignment)[1])
  final_cell_names <- prior_info[match(
    final_cell_type_assignment,
    lineage_info$Cell_type_number
  ), 1]
  final_cell_names[which(final_cell_type_assignment == 0)] <- "Unknown"
  final_result <- cbind(
    cell_type_name_assigned,
    final_cell_type_assignment,
    final_cell_names
  )
  round_name <- paste("Round", seq(1, total_rounds, by = 1))
  colnames(final_result) <- c(round_name, "Cell type number", "Final cell type")

  if (save_result) {
    filename <- paste0(project_name, "_final_cell_type_assignment.csv")
    write.csv(cbind(final_result, coords, original_exp),
      file = filename,
      row.names = FALSE
    )
  }
  return(final_result)
}
################################################################################
################################################################################