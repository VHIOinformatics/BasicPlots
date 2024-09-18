# BasicPlots
Basic plots we use in our projects

## Installation

```r
library(devtools)
install_github("VHIOinformatics/BasicPlots")
```

## Package contents

* **makePCA.R**: Makes a Principal Component Analysis from given counts table and metadata targets file and plots it. Can also be plotted interactively which may be useful for exploration when there are a lot of samples. Main function is `makePCA`, but old function to perform a 3D PCA is also included as `makePCA3D`.

* **makeHeatmap.R**: Makes a heatmap using pheatmap function, for plotting differentially expressed genes and how they cluster with condition.

* **makeVolcano**: Creates a volcano plot to visualize differentially expressed genes according to their significance and FC. Can also be plotted interactively, otherwise shows the top genes.

* **makeVenn.R**: Creates a Venn diagram from a list of named lists of genes to compare. Allows from 2 to 4 way-comparison.

* **manyClusters.R**: Creates a pdf with different dendrograms combining possible distances and linkage methods, to evaluate the most suited combination.

* **oneCluster.R**: Creates a dendrogram (hierchical clustering) with distance and linkage method of choice.

* **colorCluster.R**: Function used in `oneCluster` and `manyClusters` for colouring the dendrogram according to a vector of conditions.

* **setParameters.R**: Function to define different parameters of graphics. Used in other functions.

For more information on the parameters and usage of the functions, please check the documentation of each of the functions in R.
