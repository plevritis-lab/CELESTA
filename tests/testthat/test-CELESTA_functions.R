library(testthat)
library(dplyr)
library(devtools)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)

# Load the original CELESTA functions
# Note that some of the slot names were changed in the original to match
# the new version for comparison.
source("../CELESTA_functions_orig.R")

# Load data
load(file = "../../data/prior_marker_info.rda")
load(file = "../../data/imaging_data.rda")

compareCelesta <- function(actual, expected) {
  sapply(slotNames(actual), function(x) {
    tryCatch(
      {
        a <- slot(actual, x)
        e <- slot(expected, x)

        if (typeof(a) == "double") {
          expect_equal(e, a)
        } else {
          expect_identical(e, a)
        }
      },
      error = function(e) {
        # If there is not a matching slot (in the case where it was deleted),
        # assume that it is vacuously true that they are equal
        return(TRUE)
      }
    )
  })
}

test_that("CreateCelestaObject", {
  actual <- CELESTA::CreateCelestaObject(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )
  expected <- CreateCELESTAobj(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )
  compareCelesta(actual, expected)
})

test_that("FilterCells", {
  celesta_obj <- CELESTA::CreateCelestaObject(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )
  CelestaObj <- CreateCELESTAobj(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )

  actual <- CELESTA::FilterCells(celesta_obj,
    high_marker_threshold = 0.9,
    low_marker_threshold = 0.5
  )
  expected <- cell_filtering(
    high_marker_threshold = 0.9, low_marker_threshold = 0.5,
    CelestaObj
  )
  compareCelesta(actual, expected)
})

test_that("AssignCells", {
  celesta_obj <- CELESTA::CreateCelestaObject(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )
  CelestaObj <- CreateCELESTAobj(
    project_title = "project_title",
    prior_marker_info,
    imaging_data
  )

  celesta_obj <- CELESTA::FilterCells(celesta_obj,
    high_marker_threshold = 0.9,
    low_marker_threshold = 0.5
  )
  CelestaObj <- cell_filtering(
    high_marker_threshold = 0.9, low_marker_threshold = 0.5,
    CelestaObj
  )

  actual <- CELESTA::AssignCells(celesta_obj,
    max_iteration = 10, cell_change_threshold = 0.01,
    high_expression_threshold_anchor = high_marker_threshold_anchor,
    low_expression_threshold_anchor = low_marker_threshold_anchor,
    high_expression_threshold_index = high_marker_threshold_iteration,
    low_expression_threshold_index = low_marker_threshold_iteration
  )
  expected <- assign_cell_main(CelestaObj,
    max_iteration = 10, cell_change_threshold = 0.01,
    high_marker_threshold_anchor = high_marker_threshold_anchor,
    low_marker_threshold_anchor = low_marker_threshold_anchor,
    high_marker_threshold_iteration = high_marker_threshold_iteration,
    low_marker_threshold_iteration = low_marker_threshold_iteration
  )
  compareCelesta(actual, expected)
})
