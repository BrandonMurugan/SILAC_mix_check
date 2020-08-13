#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("united"),
                  
                  # Application title
                  
                  # titlePanel("SILAC Mixing Check"),
                  h1(id="big-heading", "SILAC/Dimethyl Mixing Check"),
                  tags$style(HTML("#big-heading{color: #E95420;}")),
                  tags$style(HTML("#sidebar-heading{color: #16bd5f;}")),
                  
                  # Sidebar with a slider input for number of bins 
                  sidebarLayout(
                      sidebarPanel(
                          h2(id="sidebar-heading", "Mixing Calculation"),
                          radioButtons("multiplicity",
                                       "SILAC Multiplicity:",
                                       choices = c("double", "triple"),
                                       selected = "double"),
                          
                          fileInput(inputId = "peptides_file",
                                    label = "Select peptides.txt file",
                                    multiple = F),
                          
                          fileInput(inputId = "mqpar_file",
                                    label = "Select mqpar.xml file",
                                    multiple = F),
                          
                          actionButton("calculate",
                                       "Calculate Mixing")
                          
                          # h2(id="sidebar-heading", "Experimental Layout"),
                          # 
                          # numericInput("total_protein",
                          #              "Input total protein (H+L)",
                          #              10,
                          #              1,
                          #              100),
                          # numericInput("target_protein",
                          #              "Input target protein amount (H+L)",
                          #              10,
                          #              1,
                          #              100),
                          # actionButton("experimentalDesign",
                          #              "Calculate Experimental Layout")
                          
                      ),
                      
                      # Show a plot of the generated distribution
                      mainPanel(
                          tabsetPanel(
                              # tabPanel("Peptides.txt", dataTableOutput("pepfile")),
                              tabPanel("Input",fluidRow(column(width = 12, HTML("<font size=4>Enter your volumes here (see help for more details)</font>")), 
                                                        hotable("hotable1")),
                                       fluidRow(column(width = 12,DT::dataTableOutput('tbl')))),
                              tabPanel("Mixing Calculation", DT::dataTableOutput("mixcalcTable"),
                                       downloadButton("download", "Download Mixing Calculations")),
                              tabPanel("Help", htmlOutput("helppage")),
                              tabPanel("About", htmlOutput("aboutpage"))
                          )
                          
                          # dataTableOutput("expDes")
                      )
                  )
)
)
