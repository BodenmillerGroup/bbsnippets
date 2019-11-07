## Cell shape description

Uses an autoencoder neural network to generate cell shape descriptors from cell/structure masks.

This is based on the following paper: https://www.biorxiv.org/content/10.1101/790162v1
Original code by Eliot McKinley (eliot.mckinley@vumc.org): https://github.com/Coffey-Lab/CellSegmentation

All code related to cell segmentation was removed to keep only the cell shape pre-processing and autoencoder sections.
The code was adapted to run directly from cell masks without the preliminary steps.

# How it works
- Put all your masks in a folder
- Open the RunAutoencoder.m file in MatLab (tested on MatLab R2018b, the Deep Learning Toolbox must be installed).
- Modify the input and output paths as indicated
- Run

# Output
A .csv file with one row per mask object and cell shape descriptors as columns.
This can be easily imported in R or Python for analysis (e.g. dimensionality reduction based on cell shapes)
