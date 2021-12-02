## Antibody titration

### Goal
__The goal of this script is to facilitate antibody titration in IMC__

This script helps identifying, for each tested marker, the positive and negative cells (based on clustering or on gating). The signal and noise levels are then calculated to help choosing the ideal concentration for each antibody.

### Antibody titration
For titration, an antibody mix is generated and a serial dilution is performed. Different sections (or different areas of the same section) are then stained with the different dilutions (relative concentrations). Finally, these sections are imaged by IMC.

The present script requires that the mcd files and ROIs have been named in a consistent manner, in particular ___the mcd file names OR the ROI description must contain a number corresponding to the dilution factor / relative concentration of the antibody mix that was used to stain the section or section area___.  

### Image preprocessing
The `.mcd` files generated in IMC can be preprocessed using [steinbock](https://github.com/BodenmillerGroup/steinbock).  
Steinbock preprocessing allows to segment the images, measure single-cell marker inteisites and export the data to a readable format.  
This script assumes that IMC data preprocessing was performed with `steinbock` using default parameters. If other settings were used, part of this script may have to be adapted.
At the end of `steinbock` preprocessing, the data have to be exported in the `AnnData` format.

For reference, the following `steinbock` command can be used (adjust the path and steinbock version):
```
$ alias steinbock="docker run -v /path/to/your/data:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/home/steinbock/.Xauthority:ro -u $(id -u):$(id -g) -e DISPLAY ghcr.io/bodenmillergroup/steinbock:latest"
$ steinbock preprocess imc images --hpf 50
$ steinbock segment deepcell --minmax
$ steinbock measure intensities
$ steinbock measure regionprops
$ steinbock export anndata --x intensities --obs regionprops
```

Full `steinbock` documentation can be found here: https://bodenmillergroup.github.io/steinbock/latest.  

### Script usage
This should be run downstream of `steinbock` script performs the following tasks:

1. Import data (exported by `steinbock` as `AnnData` object(s)).
2. Calculate UMAP and generate some quality control plots.
3. Cluster the cells.
4. Facilitate the selection of clusters corresponding to positive and negative cells using either clustering (__4.1__) or gating (__4.2__).
5. Display images to make sure that the selection of positive and negative cells is correct.
6. Calculate the signal, noise, and signal-to-noise ratio.
7. Export the titration results.

__Paragraphs 4 to 6 have to be repeated for each marker__
