# MISTy runner
# Copyright (c) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]


#' @importFrom rlang !! := .data
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mistyR is able to run computationally intensive functions
  in parallel. Please consider specifying a future::plan(). For example by running
  future::plan(future::multisession) before calling mistyR functions.")
}

#' @importFrom dplyr %>%
#' @inherit run_misty examples
#' @export
dplyr::`%>%`

#' Train MISTy models
#'
#' Trains multi-view models for all target markers, estimates the performance,
#' the contributions of the view specific models and the importance of predictor
#' markers for each target marker.
#'
#' If \code{bypass.intra} is set to \code{TRUE} all variable in the intraview
#' the intraview data will be treated as targets only. The baseline intraview
#' model in this case is a trivial model that predicts the average of each
#' target "(basically just an intercept with no slope"). 
#' If the intraview has only one variable this switch is automatically
#' set to \code{TRUE}. Which makes sense, since there would not be any predictors
#' available for the intraview
#'
#' Default values passed to \code{\link[ranger]{ranger}()} for training the
#' view-specific models: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' @param views view composition.
#' @param results.folder path to the top level folder to store raw results.
#' @param seed seed used for random sampling to ensure reproducibility.
#' @param target.subset subset of targets to train models for. If \code{NULL},
#'     models will be trained for markers in the intraview.
#' @param bypass.intra a \code{logical} indicating whether to train a baseline
#'     model using the intraview data (see Details).
#' @param cv.folds number of cross-validation folds to consider for estimating
#'     the performance of the multi-view models.
#' @param cached a \code{logical} indicating whether to cache the trained models
#'     and to reuse previously cached ones if they already exist for this sample.
#' @param append a \code{logical} indicating whether to append the performance
#'     and coefficient files in the \code{results.folder}. Consider setting to
#'     \code{TRUE} when rerunning a workflow with different \code{target.subset}
#'     parameters.
#' @param ... all additional parameters are passed to
#'     \code{\link[ranger]{ranger}()} for training the view-specific models
#'     (see Details for defaults).
#'
#' @return Path to the results folder that can be passed to
#'     \code{\link{collect_results}()}.
#'
#' @seealso \code{\link{create_initial_view}()} for
#'     starting a view composition.
#'
#' @examples
#' # Create a view composition of an intraview and a paraview with radius 10 then
#' # run MISTy for a single sample.
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#' # get the coordinates for each cell
#' pos <- synthetic[[1]] %>% select(row, col)
#'
#' # compose
#' misty.views <- create_initial_view(expr) %>% add_paraview(pos, l = 10)
#'
#' # run with default parameters
#' run_misty(misty.views)
#' @export
run_misty <- function(views, results.folder = "results", seed = 42,
                      target.subset = NULL, bypass.intra = FALSE, cv.folds = 10,
                      cached = FALSE, append = FALSE, ...) {
  normalized.results.folder <- R.utils::getAbsolutePath(results.folder)

  if (!dir.exists(normalized.results.folder)) {
    dir.create(normalized.results.folder, recursive = TRUE)
  }

  on.exit(sweep_cache())

  # Generate a vector to get all view abbreviations
  view.abbrev <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])

  # Generate header for what the coefficient file, containing the
  # coefficient from the linear meta model
  header <- stringr::str_glue("target intercept {views} p.intercept {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  #  Extract the expression values from the intraview (aka baseline model)
  expr <- views[["intraview"]][["data"]]

  # Check whether we actually have more samples (aka spots/cells) than cv folds
  assertthat::assert_that(nrow(expr) >= cv.folds,
    msg = "The data has less rows than the requested number of cv folds."
  )

  # Bypass intra was implemented if we do not want to include the intramodel,
  # which could be the case if we only have a single target in a cell
  # (thus there would be no predictors left).
  # but would this have to be ncol???
  if(nrow(expr) == 1) bypass.intra <- TRUE
  
  # Get the standard deviation for each target
  target.var <- apply(expr, 2, stats::sd)

  # If any targets have no variance, output a warning (in this case the modeling
  # would be kinda unnecessary, simply mean=intercept, slope=0 for linear model)
  if (any(target.var == 0)) {
    warning.message <- paste(
      "Targets",
      paste(names(which(target.var == 0)),
        collapse = ", "
      ),
      "have zero variance."
    )
    warning(warning.message)
  }

  # Now we create the filenames for the txt file that contain the meta model
  # coefficients and the filename for the lock.file
  coef.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt"
  )
  coef.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt.lock"
  )
  on.exit(file.remove(coef.lock), add = TRUE)

  # If append is FALSE, simply write/overwrite file with header
  if (!append) {
    current.lock <- filelock::lock(coef.lock)
    write(header, file = coef.file)
    filelock::unlock(current.lock)
  }

  # Now we will do the same thing for the file that contains there performance measures
  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"

  perf.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt"
  )
  perf.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt.lock"
  )
  on.exit(file.remove(perf.lock), add = TRUE)

  if (!append) {
    current.lock <- filelock::lock(perf.lock)
    write(header, file = perf.file)
    filelock::unlock(current.lock)
  }

  # Check which targets should be predicted
  targets <- switch(class(target.subset),
    # if numbers are supplied use the indexed targets
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    # if a character vector is supplied use these targets
    "character" = target.subset,
    # if no target.subset was specified we simply take all available targets
    "NULL" = colnames(expr),
    NULL
  )

  message("\nTraining models")
  
  # Building the model for each target (see models.R for details)
  targets %>% furrr::future_map_chr(function(target, ...) {
    target.model <- build_model(
      views, target, bypass.intra,
      seed, cv.folds, cached, ...
    )

    # After building and training the model, we basically need to extract
    # all the interesting information/data
    
    # 1. The build, trained, and evaluated meta model
    combined.views <- target.model[["meta.model"]]
    model.summary <- summary(combined.views)

    # 2. The coefficients of that meta model
    # coefficient values and p-values
    # WARNING: hardcoded column index
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])

    # Write to file as specified above, use filelock system to prevent
    # simultaneous access from different processes (only process which 
    # has the lock can proceed with writing, other processes have to wait)
    current.lock <- filelock::lock(coef.lock)
    write(paste(target, paste(coeff, collapse = " ")),
      file = coef.file, append = TRUE
    )
    filelock::unlock(current.lock)

    # 3. The importances of each predictor in each model
    # extract the importances as returned by the ranger random forest implementation
    # TODO: What are those importances actually?
    target.model[["model.views"]] %>% purrr::walk2(
      view.abbrev,
      function(model.view, abbrev) {
        model.view.imps <- model.view$variable.importance
        targets <- names(model.view.imps)

        imps <- tibble::tibble(
          target = targets,
          imp = model.view.imps
        )
        
        # Write tibble to csv
        readr::write_csv(
          imps,
          paste0(
            normalized.results.folder, .Platform$file.sep,
            "importances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )

    # 4. The performance measures of the baseline linear model (based on intraview)
    # and the meta linear model (based on all views)
    # Output warning if negative performance measures detected
    if (sum(target.model[["performance.estimate"]] < 0) > 0) {
      warning.message <-
        paste(
          "Negative performance detected and replaced with 0 for target",
          target
        )
      warning(warning.message)
    }

    # 4. a) Replace all negative values by 0
    performance.estimate <- target.model[["performance.estimate"]] %>%
      dplyr::mutate_if(~ sum(. < 0) > 0, ~ pmax(., 0))
    
    # 4. b) Perform ttest for the RMSE from baseline model and multiview model
    performance.summary <- c(
      performance.estimate %>% colMeans(),
      tryCatch(stats::t.test(performance.estimate %>%
        dplyr::pull(.data$intra.RMSE),
      performance.estimate %>%
        dplyr::pull(.data$multi.RMSE),
      # if multiview model explains more variability of the data, then
      # RMSE should be smaller, thus our alternative is greater
      alternative = "greater"
      # and extract p.value
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of RMSE performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      }),
      # 4. b) Perform ttest for the R2 from baseline model and multiview model
      tryCatch(stats::t.test(performance.estimate %>%
        dplyr::pull(.data$intra.R2),
      performance.estimate %>%
        dplyr::pull(.data$multi.R2),
      # if multiview model explains more variability of the data, then
      # R2 should be greater, thus our alternative is less
      alternative = "less"
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of R2 performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      })
    )

    current.lock <- filelock::lock(perf.lock)
    write(paste(target, paste(performance.summary, collapse = " ")),
      file = perf.file, append = TRUE
    )
    filelock::unlock(current.lock)

    return(target)
  }, ..., .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

  return(normalized.results.folder)
}
