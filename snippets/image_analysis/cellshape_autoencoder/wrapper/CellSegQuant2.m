function CellSegQuant2(dirmask, outdir, maskstring, fract)
%Eliot McKinley - 2018-11-26
%eliot.mckinley@vumc.org
%
% No returned output - function saves quantification of cell shapes
%
% dirmask: directory containing cell masks
% outdir: output directory
% maskstring: Common character string used to identify mask files
% fract: fraction of cells to train on

%% Parse Directory supplied for cell segmentation
if ~exist(outdir, 'dir')
    fprintf('No');
    mkdir(outdir);
end

MaskList=dir([dirmask maskstring]);
MaskList={MaskList.name};
% Masks=regexp(MaskList, '_', 'split'); %split characters before and after delimiter
% Masks=cat(3,Masks{:}); %reorganize cell
% Masks=unique(squeeze(Masks(1,1,:))); %get image names


%% Segmentation and Quantification for each position
for i=1:length(MaskList)
    %% Shape Pre-processing
    %file names
    mask=[dirmask char(MaskList(i))];
    fn=extractBefore(char(MaskList(i)),".");
    
    %load the final cell segmentation image, extract the cells, and save to .mat file
    if exist(mask, 'file') && ~exist([outdir fn '.mat' ], 'file')
        fprintf('\nCell Shape Pre-Processing; ');
        CellImages=imread(mask);
        CellImages=cell_shape_images(CellImages);
        save([outdir fn '.mat' ], 'CellImages')
    end
    
    %concatenate all the position cell shapes and
    if i==length(MaskList)
        if ~exist([outdir 'autoencoder.mat'], 'file') && ~exist([outdir 'encoded_cells.csv'], 'file')
            fprintf('\nTraining Autoencoder; ')
            FileList=dir([outdir '*.mat']);
            FileList=fullfile({FileList.folder}.', {FileList.name}.');
            %run autoencoder with specified percent of training data
            [autoencoder, trainList]=CellShapeAutoencoder(FileList, fract);
            
            save([outdir 'autoencoder.mat'], 'autoencoder')
            writetable(trainList, [outdir 'encoded_cells.csv']);
        end
    end
    fprintf('\n')
end

