
ui <- fluidPage(
  
  title = 'Gene Expression Profiler',
  
  
  h1('Expression'),
  
  fluidRow(
    column(4, plotOutput('x1', height = 300)),
    column(4, plotOutput('x2', height = 300)),
    column(4, plotOutput('x3', height = 300)),
  ),

  hr(),
 
  fluidRow(
    column(4, plotOutput('x4', height = 300)),
    column(4, plotOutput('x5', height = 300)),
    column(4, plotOutput('x6', height = 300)),
  ),
  
  hr(),
  
  fluidRow(
    column(12, DT::dataTableOutput('x7')),
  )
  
)