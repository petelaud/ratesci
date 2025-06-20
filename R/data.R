#' Corticosteroids in acute traumatic brain injury: updated systematic review of
#' randomised controlled trials
#'
#' Data from systematic review of the effect on mortality of corticosteroids in
#' traumatic brain injury (reported with MRC CRASH trial results, Roberts et al.
#' 2001)
#'
#' @format A data frame with five variables:
#' \describe{
#'   \item{study}{Study author and year}
#'   \item{event.steroid}{Number of deaths in steroid-treated group}
#'   \item{n.steroid}{Number of patients in steroid-treated group}
#'   \item{event.control}{Number of deaths in control group}
#'   \item{n.control}{Number of patients in control group}
#'  }
#' @source <https://pubmed.ncbi.nlm.nih.gov/15474134>
"crash"

#' Meta-analysis of the effect of cisapride for treatment of non-ulcer dyspepsia
#'
#' Data from systematic review of the effect of cisapride for treatment of
#' non-ulcer dyspepsia (Hartung & Knapp 2001)
#'
#' @format A data frame with five variables:
#' \describe{
#'   \item{study}{Study author}
#'   \item{event.cisa}{Number of events (successes) in cisapride-treated group}
#'   \item{n.cisa}{Number of patients in cisapride-treated group}
#'   \item{event.plac}{Number of events (successes) in placebo group}
#'   \item{n.plac}{Number of patients in placebo group}
#'  }
#' @source \doi{10.1002/sim.1009}
"cisapride"

#' Systematic review of the effect of graduated compression stockings for
#' prevention of DVT
#'
#' Data from systematic review of the effect of graduated compression stockings
#' for prevention of DVT (Roderick et al. 2005)
#'
#' @format A data frame with five variables:
#' \describe{
#'   \item{study}{Study author}
#'   \item{event.gcs}{Number of events (DVTs) in GCS-treated group}
#'   \item{n.gcs}{Number of patients in GCS-treated group}
#'   \item{event.control}{Number of events (DVTs) in control group}
#'   \item{n.control}{Number of patients in control group}
#'  }
#' @source \doi{10.3310/hta9490}
"compress"
