################################# DEPENDENCIES #################################
packages_cran <-
  c("curl",
    "dplyr",
    "purrr",
    "readr",
    "rvest")
packages_bioconductor <- NULL
packages_github <- NULL

################################## VARIABLES ##################################

ZENODO_METADATA <- "https://doi.org/10.5281/zenodo.6378223"

PATTERN_METADATA_STRUCTURES <- "structure_metadata"

PATH_EXPORT_CFM <- "smiles4cfm.txt"

################################## FUNCTIONS ##################################

#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
check_and_load_packages <- function(cran = packages_cran,
                                    bioconductor = packages_bioconductor,
                                    github = packages_github) {
  installed_packages <- rownames(installed.packages())
  installed_packages_cran <- cran %in% installed_packages
  installed_packages_bioconductor <-
    bioconductor %in% installed_packages
  installed_packages_github <- github %in% installed_packages
  
  if (!is.null(bioconductor)) {
    cran <- cran |>
      append("BiocManager")
  }
  if (!is.null(github)) {
    cran <- cran |>
      append("remotes")
  }
  
  if (any(installed_packages_cran == FALSE)) {
    install.packages(cran[!installed_packages_cran])
  }
  if (any(installed_packages_bioconductor == FALSE)) {
    BiocManager::install(bioconductor[!installed_packages_bioconductor])
  }
  if (any(installed_packages_github == FALSE)) {
    lapply(X = github[!installed_packages_github], FUN = remotes::install_github)
  }
  
  return(lapply(c(
    cran,
    bioconductor,
    gsub(
      pattern = ".*/",
      replacement = "",
      x = github
    )
  ),
  require,
  character.only = TRUE) |>
    invisible())
}

#' Title
#'
#' @param session
#' @param text
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
follow_next <- function(session,
                        text = "Next",
                        ...) {
  link <- rvest::html_element(x = session,
                              xpath = sprintf("//*[text()[contains(.,'%s')]]", text))
  
  url <- rvest::html_attr(link, "href") |>
    trimws() |>
    gsub(pattern = "^\\.{1}/", replacement = "")
  
  message("Navigating to ", url)
  
  rvest::session_jump_to(session, url, ...)
}

#' Title
#'
#' @param url
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
get_last_file_from_zenodo <- function(url, pattern) {
  file <- rvest::session(url = url) |>
    follow_next(text = pattern) |>
    purrr::pluck("url") |>
    curl::curl_download(destfile = tempfile()) |>
    readr::read_delim()
  
  return(file)
}

################################### PROGRAM ###################################

start <- Sys.time()

#' (Down)load dependencies
check_and_load_packages()

#' Get last version from Zenodo
last_lotus_structures <-
  get_last_file_from_zenodo(url = ZENODO_METADATA,
                            pattern = PATTERN_METADATA_STRUCTURES)

#' Minimal cleaning and formatting
lotus_4cfm <- last_lotus_structures |>
  dplyr::filter(structureCleaned_exactMass > 50 &
                  structureCleaned_exactMass < 2000) |>
  dplyr::select(structure_inchikey_2D = structureCleaned_inchikey2D,
                structure_smiles_2D = structureCleaned_smiles2D) |>
  dplyr::filter(!grepl(pattern = "He", x = structure_smiles_2D)) |>
  dplyr::distinct()

#' Export
readr::write_delim(x = lotus_4cfm,
                   file = PATH_EXPORT_CFM,
                   col_names = FALSE)

end <- Sys.time()

message("Finished in ", format(end - start))
