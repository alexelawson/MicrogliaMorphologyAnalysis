Project and Dataset Description:

Data used in analysis is contained in the data.csv file and includes morphometric descriptions of individual microglia. 
    Rows: individual microglia 
    Columns: attributes (e.g. area, number of branches, circularity...)

R-Package used with data set and relevant functions can be found at the following link: https://ciernialab.github.io/MicrogliaMorphologyR/reference/index.html

The process for this project was adapted from a pipeline developed by Kim J et al (2024). 

The R-script (analysis.R) imports the data, and performs PCA analysis, clustering, and subsequent statistical analysis. In order to be run on local computer, directories must be altered and relevant packages must be installed. 
Reference to the pipeline by Kim J. et al will be helpful if attempting to reproduce. 
See Wiki Page methods for more detailed description of process if attempting to implement. 