#' ShinyBrainSABER
#' 
#' This function runs the Shiny app for the BrainSABER workflow.
#' @return NULL (Invisibly)
#' @import shiny
#' @export
#'
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'     options(device.ask.default = FALSE)
#'     runShinyBrainSABER()
#' }
runShinyBrainSABER <- function(){
    appDir <- system.file("shiny-examples", "shiny_brainsaber",
                            package = "BrainSABER")
    if (appDir == "") {
        stop("Example directory not found. Try re-installing `BrainSABER`.",
            call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")
}
