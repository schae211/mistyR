# mistyR model training functions
# Copyright (c) 2020 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Train a multi-view model for a single target
#'
#' Trains individual Random Forest models for each view, a linear meta-model
#' using the out-of-bag predictions of the view specific models and estimates
#' the overall performance by cross-validation.
#'
#' Default values passed to \code{\link[ranger]{ranger}()} for training the
#' view-specific models: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' @inheritParams run_misty
#'
#' @param target name of the target marker to train models for.
#'
#' @return A list containing the trained meta-model, a list of trained
#' view-specific models and performance estimates.
#'
#' @noRd
build_model <- function(views, target, bypass.intra = FALSE, seed = 42,
                        cv.folds = 10, cached = FALSE, ...) {
  
  cache.location <- R.utils::getAbsolutePath(paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  ))

  if (cached && !dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  expr <- views[["intraview"]][["data"]]
  
  # We will basically build a model for each target, using all the other
  # coefficients (meaning genes or whatever else)
  target.vector <- expr %>% dplyr::pull(target)

  # merge ellipsis with default algorithm arguments
  # See ranger::ranger documentation for these arguments
  algo.arguments <- list(
    num.trees = 100, importance = "impurity",
    verbose = FALSE, num.threads = 1, seed = seed,
    dependent.variable.name = target
  )

  ellipsis.args <- list(...)
  ellipsis.args.text <- paste(names(ellipsis.args), ellipsis.args,
    sep = ".", collapse = "."
  )
  
  # So if there are ellipsis ("additional") arguments they are megered here with
  # the algo arguments
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments <- rlist::list.merge(algo.arguments, ellipsis.args)
  }

  # returns a list of models
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
      # For each view we will do the following
      
      # 1. First let's create a file name where the model will be saved
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", view[["abbrev"]], ".", target,
          ".par", ellipsis.args.text, ".rds"
        )

      if (file.exists(model.view.cache.file) & cached) {
        # if cached is true and a file with the right name exists,
        # the model view will simply be read from storage
        model.view <- readr::read_rds(model.view.cache.file)
        
      } else {
        # If we want to bypass intra and currently have the intra view
        # transformed data will simply be the target vector with .novar = 0
        if ((view[["abbrev"]] == "intra") & bypass.intra) {
          transformed.view.data <-
          # Why do we assign the zero variance here?
          tibble::tibble(!!target := target.vector, ".novar" = 0)
          # decomposing the above command
          # We want to create a new tibble where the first column is named according to the target
          # If we would write target = target.vector, the column would literally be named target
          # Therefore, we have to use ":=" that expect an expression (and not a string) on the left side
          # and !!target which forces the evaluation of target (get the value from the variable)
          # so !! basically creates and expression
        } else {
          # We get the data from the view and add the target vector (basically our dependent variable)
          transformed.view.data <- view[["data"]] %>%
            dplyr::mutate(!!target := target.vector)
        }
        
        # Now we use do.call to evaluate the following expression
        model.view <- do.call(
          # Ranger is a random forest implementation
          ranger::ranger,
          c(
            list(data = transformed.view.data),
            algo.arguments
          )
        )
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }

      return(model.view)
    })

  # Now for each model we extract the OOB predictions
  oob.predictions <- model.views %>%
    purrr::map(~ .x$predictions) %>%
    rlist::list.cbind() %>%
    tibble::as_tibble(.name_repair = make.names) %>%
    # And we also add again the target.vector that we want to predict
    dplyr::mutate(!!target := target.vector)

  # We use the OOB predictions to train the linear model
  combined.views <- stats::lm(
    # This is a neat trick to use all other columns except for target as predictors
    stats::as.formula(paste0(target, "~.")),
    oob.predictions
  )

  # To assess the linear model we will use k-fold cross validation
  # In the first step we create k folds
  test.folds <- withr::with_seed(
    seed,
    caret::createFolds(target.vector, k = cv.folds)
  )

  # And we also extract the intraview OOB predictions so we can compare
  # linear models trained with only these compared to models trained with
  # OOB predictions from all random forests (aka from all views)
  intra.view.only <-
    model.views[["intraview"]]$predictions %>%
    tibble::enframe(name = NULL) %>%
    dplyr::mutate(!!target := target.vector)

  # Now we loop (with map) over all folds
  performance.estimate <- test.folds %>% purrr::map_dfr(function(test.fold) {
    
    # Train linear model based on intraview RF OOB predictions
    meta.intra <- stats::lm(
      stats::as.formula(paste0(target, "~.")),
      intra.view.only %>% dplyr::slice(-test.fold)
    )
    
    # Train linear model based on all views RF OOB predictions
    meta.multi <- stats::lm(
      stats::as.formula(paste0(target, "~.")),
      oob.predictions %>% dplyr::slice(-test.fold)
    )

    # Get the predictions from the linear model to calculate stats 
    intra.prediction <- stats::predict(meta.intra, intra.view.only %>%
      dplyr::slice(test.fold))
    multi.view.prediction <- stats::predict(meta.multi, oob.predictions %>%
      dplyr::slice(test.fold))
    
    # Calcualte RMSE (root of mean squared error) and R2 (proportion of explained variance)
    # for intra and all_view (aka meta view)
    intra.RMSE <- caret::RMSE(intra.prediction, target.vector[test.fold])
    intra.R2 <- caret::R2(intra.prediction, target.vector[test.fold],
      formula = "traditional"
    )

    multi.RMSE <- caret::RMSE(multi.view.prediction, target.vector[test.fold])
    multi.R2 <- caret::R2(multi.view.prediction, target.vector[test.fold],
      formula = "traditional"
    )

    # Store everything in a tibble
    tibble::tibble(
      intra.RMSE = intra.RMSE, intra.R2 = intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = multi.R2
    )
  })

  final.model <- list(
    meta.model = combined.views,
    model.views = model.views,
    performance.estimate = performance.estimate
  )

  return(final.model)
}
