
# GNU Free Documentation License (GNU FDL)
# Copyright ©  2022 Alberto Joven Álvarez.
# Permission is granted to copy, distribute and/or modify this document under 
# the terms of the GNU Free Documentation License, Version 1.3 or any later
# version published by the Free Software Foundation; with no Invariant Sections, 
# no Front-Cover Texts, and no Back-Cover Texts. 
# A copy of the license is included in the section entitled "GNU Free Documentation License".


library(shiny)

shinyUI(fluidPage(

  tabsetPanel(
    tabPanel("TCGA diseasecodes classifier. Neural Network",
             
             sidebarLayout(
               sidebarPanel(
                 fileInput("betas", label="Text File with betas values", accept=c(".csv", ".txt")),
                 checkboxInput("divide", label="Your beta values should be divided by 1000?"),
                 checkboxInput("confussion", label="Ticks to get Confussion Matrix"),
                 fileInput("etiqueta", label="Labels of your samples"),
                 submitButton("Submit"),
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "author: Alberto Joven Álvarez" ),
                 tags$br(),
                 HTML(
                   " <a href='mailto:ajoven@uoc.edu'> mail: ajoven@uoc.edu</a> "
                 ),
                 tags$br(),
                 tags$br()
               ),
               
               mainPanel(
                 tags$html(
                   "This app uses a Neural Network to predict one of the 33 classes of tumors 
    codified by the TCGA project, using the methylation beta data obtained with the Illumina 
    450k array. Samples of morphologically normal adjacent tissue are grouped into an 
    additional class called Control."),
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "The classifier has been trained with 9267 samples obtained from the TCGA project and uses the 1000 
    probes with the highest standard deviation value as explanatory variables. " ),
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "The file with the betas values must have format .txt, The first column will contain the names 
    of the probes and will be named as cpg and the first row the sample identifiers. "),
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "If you include a file with the sample labels, this will also have .txt or .csv format 
    and the first column will contain the sample identifiers with a column name equal to id_clx." ),
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "The second column will contain the label: normal sample, or type of tumor, the name of this
    second column must be equal to type. ") ,
                 tags$br(),
                 tags$br(),
                 tags$html(  
                   "In the image below we see the formats that the files must have. The image on the left shows the beta 
    values of two samples, and the one on the right shows the values of the labels. " ) ,
                 tags$br(),
                 tags$br(),
                 tags$html(
                   "By default, the classifier displays the tumor class (TCGA category) and the probability of success.") ,
                 tags$br(),
                 tags$br(),
                 tags$html(  
                   "If the Ticks to get Confusion Matrix check is marked, the label file must be indicated 
    and the classifier will present the confusion matrix." ) ,
                 tags$br(),
                 tags$br(),
                 tags$img(src='betas.png'),
                 tags$br(),
                 tags$html(
                   "If the sample file contains only one sample, the sensitivity analysis is offered. A table is presented
               in which, for each type of probes, according to the classification proposed by Price et all 2013:
               HC high-density CpG Island, IC Intermediate density island, ICshore regions of intermediate-density 
               CpG island that border HC and LC probes low-density CpG regions. It is calculated what is the 
               probability of the initially chosen class if all the probes of that group had a value of zero or one? 
               and the newly selected class is displayed.Likewise, the analysis is carried out considering 
               the probes grouped by chromosome."
                 ),
                 tags$br(),
                 tableOutput("table1"),
                 tags$br(),
                 tableOutput("table2"),
                 tags$br(),
                 tableOutput("table3"),
                 tags$br(),
                 tableOutput("table4")
               )
             )
    ),
    tabPanel("Detail of the probes used in the classifier",
             mainPanel(
               tags$html("Table with probes by type of probe"),
               tags$br(),
               verbatimTextOutput("table5"),
               tags$br(),
               tags$html("Table with probes by cromosome"),
               tags$br(),
               verbatimTextOutput("table6"),
               tags$br(),
               tags$html("Detail of the probes used in the classifier"),
               dataTableOutput("table7")
             )
    ),
    tabPanel("Diseases Codes TCGA",
             mainPanel(
               tags$html("Table with diseases codes TCGA"),
               tags$br(),
               tableOutput("table8")

               )
             )
  )
)
)

