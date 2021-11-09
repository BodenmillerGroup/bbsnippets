## Heatmaply
Explains how to generate an interactive heatmap from a `SingleCellExperiment` object.
Includes a `summarizeHeatmap` function that calculates the median counts by channel and by cluster (or any other variable in `colData(sce)`).  
For more information about heatmaply, see the [heatmaply vignette](https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html)
and the original [heatmaply publication](https://doi.org/10.1093/bioinformatics/btx657).

Requirements
- Packages: `SingleCellExperiment`, `heatmaply`, `data.table`.  
- An SCE object with the channels as `rownames(sce)`.
