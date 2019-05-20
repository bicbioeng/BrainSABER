#' ShinyBrainSABER
#' 
#' This function runs the Shiny app for the BrainSABER workflow.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   library(BrainSABER)
#'   runShinyBrainSABER()
#' }
runShinyBrainSABER <- function(){
  appDir <- system.file("shiny-examples", "shiny_brainsaber",
                        package = "BrainSABER")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `BrainSABER`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
