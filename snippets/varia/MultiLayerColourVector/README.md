### Generation of multi-layered colour vectors

This script provides a function which generates a hierarchical colour vector for classes, sub-classes and sub-su-blasses.
To load the function, please clone this repository and source the file in your R script:

```{r}
source(~/path/to/bbsnippets/snippets/varia/MultiLayerColourVector.R)
```

The function takes three arguments: `layer1`, `layer2` and `layer3`.
`layer1` is a vector containing the major class names (e.g. tumour/stroma or Bcells/Tcells/Macrophages).
`layer2` is a named list, where each entry corresponds to an entry in `layer1`.
Furthermore, each entry is a vector containing the sub-classes (e.g. Tcells = c("CD4_Tcells", "CD8_Tcells"))
`layer3` is a named list, where each entry corresponds to an entry in one of the vectors of layer2.
Similar to `layer2`, each entry is a vector containing the sub-sub-classes (e.g. CD4_Tcells = c("Th1_Tcells", "Th2_Tcells"))

