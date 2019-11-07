warning off
clear all
close all

% Directory containing the masks
dirmask='data\masks\';

% Output directory
outdir='data\shapes\';

% Common character string (used to identify mask files in dirmask)
maskstring='*isletmask.tiff';

%% Process cell shapes and run the autoencoder
% CellSegQuant2(mask directory, output directory, mask string, fraction of
% cells to train on)
CellSegQuant2(dirmask,outdir,maskstring,0.2)