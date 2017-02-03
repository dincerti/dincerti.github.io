---
layout: apps
title: 
---
* TOC
{:toc}

### Overview
This page provides additional information about our R Shiny skeletal dysplasias web application. It is a supplement to the education exhibit "A Shiny New World: Creating Your Own Radiology Decision App Using R" presented at the Radiological Society of North America 2016 Annual Meeting in Chicago, Illinois. It provides additional details on how to use R to create interactive HTML tables. The hope is that this will be a useful example for those new to Shiny, although readers are definitely encouraged to read through Shiny's official [tutorials](http://shiny.rstudio.com/tutorial/). Our example describes three files, `global.R`, `server.R`, and `ui.R`, that are used to create a Shiny app. 

### Global Environment
The `global.R` script creates R ojects that are shared accross both `server.R` and `ui.R`. This script is not strictly necessary, but it has a number of advantages. First, it allows you to keep all code that creates global objects in one place. Second, in many cases, you will want to reference the same objects in both `server.R` and `ui.R`. For example, this might be the case if you want to use an object for both an R computation and to set input values in a Shiny widget. Third, objects created in `global.R` are not loaded each time a new user uses the app, which is especially important if you need to load a large dataset. 

Lets now take a look at `global.R` in our Skeletal Dysplasia app. 

{% highlight r %}
library("shiny")
library("data.table")
diseases <- fread("data/lachman_table.txt")
features <- sort(unique(diseases$feature))
{% endhighlight %}

The first two lines load the shiny and data.table packages respectively. The third line uses the data.table function `fread`--which is a much faster version of base R's `read.csv` and `read.txt`--to load in our dataset into a data.table object named diseases. Let's examine the first 6 rows of diseases. 

{% highlight text %}
                                      feature                               disease
1: Acro-osteolysis, Distal Phalangeal Erosion             Acrodystrophic neuropathy
2: Acro-osteolysis, Distal Phalangeal Erosion                             Acrogeria
3: Acro-osteolysis, Distal Phalangeal Erosion Acroosteolysis-proximal symphangilism
4: Acro-osteolysis, Distal Phalangeal Erosion           Arteriosclerosis obliterans
5: Acro-osteolysis, Distal Phalangeal Erosion              Bureau-Barriere syndrome
6: Acro-osteolysis, Distal Phalangeal Erosion                Carpal tunnel syndrome
{% endhighlight %}

This is the dataset that we will use in our app to match features to diseases. We also want a vector of the unique features in the dataset for our users to select from. We create this vector, in alphabetical order, with the fourth line: `features <- sort(unique(diseases$feature))`.

### User Interface
`ui.R` contains the code to create the user interface. The Shiny package provides a number of R functions that convert R code into html. To illustrate, consider our `ui.R` script below. 

{% highlight r %}
# title
title <- titlePanel("Skeletal Dysplasias")

# select inputs
select <- selectInput("features", "Features", 
                       choices = features,
                       multiple = TRUE,
                       selected = NULL)
sidebar <- sidebarPanel(select)

# main panel
main <- mainPanel(dataTableOutput("Table")) 

# UI
fluidPage(title, sidebar, main)
{% endhighlight %}

The `titlePanel` function creates a header for our webpage by transating our R code to the HTML syntax <h2>Skeletal Dysplasias</h2>. 

One major advantage of Shiny is that it comes with a number of predesigned widgets. In our app, we use the `selectInput` widget, which allows users to select inputs from a list of options. In the `selectInput` function above, use the object features that we created in `global.R` to populate the list of options (`choices = diseases`). Again, this is an advantage of using `global.R`. We also set `multiple = TRUE`, which means that users can select more than one feature at a time, and `selected = NULL`, because we don't want any features to be initially selected for a user. 

Once the inputs are selected by the user, we need to display some output. This is the purpose of the `dataTableOutput(("Table"))` function, which takes an HTML table created in `server.R` named "Table" and displays it to the user. 

Three other functions, `sidebarPanel`, `mainPanel`, and `fluidPage` simply tell our web browser how to display our selectInput widget and our HTML table on the page. `fluidPage` is especially useful because it automatically scales the components on the page to fit the browser width.

### Server
The final piece of the puzzle is `server.R`, which is where all the R computations lie. This is one of the main advantages of Shiny, because it makes it very easy to use R to make server side calculations on a web page. R has a vast array of statistical and compuational algorithms, so this integration prevents the web developer from having to write code for an algorithm in another language such as javascript. 

Even in our simple app, we are able to use R's data.table package to quickly manipulate a dataset based on a user's selections. Our aim was to help radiologists quickly determine whether a patient would be likely to have a particular given their clinical findings (which we refer to as features.) To do this, we used the `server.R` code below to count the number of (user chosen) features associated with each disease. 

{% highlight r %}
shinyServer(function(input, output) {
  
  output$Table <- renderDataTable({  
    x <- diseases[feature %in% input$features]
    x <- x[, .N, by = disease]
    x <- x[order(-N)]
    setnames(x, names(x), c("Disease", "# of Features"))
    return(x)
  }, options = list(pageLength = 10, searching = TRUE, searching = 1,
                    language = list(emptyTable = "No disease match found")))
})
{% endhighlight %}

The workhorse function in the code above is `renderDataTable`, which converts an R data.frame, matrix, or data.table into an HTML table from the JQuery Javascript plugin DataTable. We named this table "Table", which recall from above, we display in our `ui.R` using `dataTableOutput("Table")`. 

To create the "Table" object, we manipulate the diseases data.table according to the features selected by the user (`input$features`). In particular, the line, `x <- diseases[feature %in% input$features]` subsets the dataset so that it only includes the user chosen features. The second line `x <- x[, .N, by = disease]`, counts the number of features associated with each disease. This allows us to determine which diseases are most likely to be associated with a particular set of features. 

The remaining lines of code make the data.table more presentable. To be exact we sort the table so the diseases with the most matching features appear at the top of the table and give the columns reasonable names. 

It's also worth noting that since we are using the DataTable JQuery plugin, we can pass arguments from R to the DataTable. For example, in our app, we set the default text to display in an empty table as "No disease match found". A full list of options can be found by visting the DataTable [website](https://datatables.net/).

### Concluding Thoughts
Our hope is that this app will help radiologists efficiently diagnose patients given their clinical findings. Our aim is to help others create similar apps to help clinicians. Shiny is particularly advantageous because it allows developers to integrate R's computational power with the web. We used R to quickly manipulate a dataset given user chosen options, but R's vast array of statistical algorithms could also be used for more complicated tasks, such as predicting the probability having a disease given clinical findings. 

