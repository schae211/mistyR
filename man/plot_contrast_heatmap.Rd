% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_contrast_heatmap}
\alias{plot_contrast_heatmap}
\title{Plot heatmap of local contrast between two views}
\usage{
plot_contrast_heatmap(misty.results, from.view, to.view, cutoff = 1)
}
\arguments{
\item{misty.results}{a results list generated by
\code{\link{collect_results}()}.}

\item{from.view, to.view}{abbreviated name of the view.}

\item{cutoff}{importance threshold. Importances below this value will
be colored white in the heatmap and considered as not relevant.}
}
\value{
The \code{misty.results} list (invisibly).
}
\description{
The heatmap shows the interactions that are present and have importance above
a \code{cutoff} value in the \code{to.view} but but not in the \code{from.view}.
}
\examples{
all.samples <- list.dirs("results", recursive = FALSE)

misty.results <- collect_results(all.samples)

misty.results \%>\%
  plot_contrast_heatmap("intra", "para.10")

misty.results \%>\%
  plot_contrast_heatmap("intra", "para.10", cutoff = 0.5)
}
\seealso{
\code{\link{collect_results}()} to generate a
    results list from raw results.

Other plotting functions: 
\code{\link{plot_contrast_results}()},
\code{\link{plot_improvement_stats}()},
\code{\link{plot_interaction_communities}()},
\code{\link{plot_interaction_heatmap}()},
\code{\link{plot_view_contributions}()}
}
\concept{plotting functions}
