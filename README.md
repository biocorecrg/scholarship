# Shiny applications

## Heatmaps

This tool allows you to generate and customize heatmaps.<br>
The main input is a text file of numbers (e.g. gene expression), with at least one column of annotation. The tool is able to detect the non-numerical columns and proposes them as annotation<br>
The data frame should not contain more than 10000 rows. If it does, you can select an annotation's column to make a selection of genes to display: you can upload a file with a selection of genes in order to filter the data. You can also directly paste the selection in the text area.<br>
The "display" button will be disabled if both the data frame and the gene selection have more than 2000 rows.
Further customization:
* Choose the colors: lower, middle and upper
* Scale the data by row or by column
* Hide row label
* Change the size of the labels.
Once you click on display, just wait (depending on the size of the input, clustering can take a bit of time). Clicking again will hide the download options.<br>

### Pending
* Saving files not working. Is it browser dependent?

## Venn diagrams

This tool allows you to generate a Venn diagram using up to 4 lists.<br>
You can either upload a file or paste your selection in the text area (one element per row).<br>
To add a new file/selection, click the + button.<br>
After you display the plot (one file/selection is required for the button to appear), the graph appears along with the options to download, to select section, and to show a list with data selected. <br>
You can change the name of the lists and the colors of the different areas of the plot.<br>

## Scatter plots

This tool allows you to generate a scatter plot based on 2 columns of an input file.<br>
As the uploaded file is read, you can select the columns you wish to represent on the x and y axis. Only the columns containing numbers are available for selection.<br>
You can set up the limits, change the name and size of axis labels or title, choose the colour / shape / size of the points, and finally group the expressed points by colour or shape.
The button "Display" will be enable to click unless the file contains more than 5000 rows.<br>
You will be able to show the correlation and p-value inside the plot. 

### Pending
* Fitted line
* More options ?

## For all tools
You can select the format of the file you want to save the plot in, or directly download png of the graph in the plot widget.<br>
You can **save the session** and restore it later on, using the provided link from "Save Session" button.

## Pending

* Integrate PCA app
* Split app in different [partials](https://github.com/jcheng5/shiny-partials) for a clearer script.


