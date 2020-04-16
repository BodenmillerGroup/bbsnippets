### Proper scatterplots with millions of cells

This gives strategies to properly plot scatterplots (eg TSNE plots) using millions of cells.

For ggplot2 it recommends "geom_bin2d" (or "geom_hexbin") to bin the data into tiny bins and visualize the
log counts as density. It also shows how to still use points instead of bins/hexbins to plot.

To visualize markers it recommends to use "stat_summary_2d" with either "mean" or a function to
randomly sample point values in each bin.

This approach is found to be much faster than using the approach used in 'geom_pointdensity'.
It also can be easily implemented in Matplolib or Plotnine when using Python. For Python additionally
Datashader (https://datashader.org/) is recommended.

This is a mirror of: https://github.com/votti/vz-varia/tree/master/heatscatter_example

Try this code also in Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/votti/vz-varia/master?urlpath=rstudio)

[filename](heatscatter_example.html ':include :type=iframe width=100% height=800px')

