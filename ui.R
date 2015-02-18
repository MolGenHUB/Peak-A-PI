library(shiny)
shinyUI(fluidPage(
  titlePanel("Peak-A-PI"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        helpText(h4("Input:")),
#         textInput("wd", label = "Choose Directory"),
        uiOutput("uibam"),        
        uiOutput("uichr")
      ),
      wellPanel(
        helpText(h4("Options:")),   
        sliderInput(
          inputId = "bg", 
          label = "Minimal Ends", 
          min = 5, 
          max = 200, 
          value = 40, 
          step = 5),
        sliderInput(
          inputId = "se", 
          label = "Stacked End",
          min = 0.30, 
          max = 1.00, 
          value = 0.50,
          step = 0.05),
        sliderInput(
          inputId = "ft", 
          label = "Flat Top", 
          min = 0.50, 
          max = 1.00, 
          value = 0.80, 
          step = 0.05),
        selectInput(
          inputId = "end", 
          label = "select End", 
          choices = list("5' end" = 1, "3' end" = 2, "5' or 3' end" = 3, "5' and 3' ends" = 4), 
          selected = 3)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "About", 
          includeMarkdown("docs/README.md")
        ),
        
        tabPanel(
          title = "Result", 
          wellPanel(
            helpText("It may takes several minutes to get the coverage infomation."),
#             helpText("Currently the following parameters are used!"),
            downloadButton("dl", "Download")
          ),
          dataTableOutput("bed")
        ),
        
        tabPanel(
          title = "Plot", 
          wellPanel(
            helpText(strong("Note:")),
            helpText("It may takes several minutes to get the coverage infomation."),
            helpText("Be patient!")
          ),
          plotOutput("plot")
        )
        
#        tabPanel("nbed",tableOutput("table"))
      )
    )
  )
))
