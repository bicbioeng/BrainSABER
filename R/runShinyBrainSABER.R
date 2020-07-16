#' Run the shiny application for the BrainSABER package workflow.
#' 
#' @return NULL (Invisibly)
#' @export
#' 
#' @examples 
#' if (interactive()) {
#'     options(device.ask.default = FALSE)
#'     runShinyBrainSABER()
#' }

runShinyBrainSaber <- function(){
  # appDir <- system.file(file.path("inst","Shiny"), package = "BrainSABER")
  # if (appDir == "") {
  #   stop("Example directory not found. Try re-installing BrainSABER.",
  #        call. = FALSE)
  # }
  # shiny::runApp(appDir, display.mode = "normal")
  
  ui = shinyAppUI
  server= shinyAppServer
  options(shiny.maxRequestSize=120*1024^2)
  
  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}