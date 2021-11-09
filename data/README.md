## Example data

Folder to collect small example data sets used in the snippets.

### pancreas_example_sce.rds

`SingleCellExperiment` object containing a subset of the [pancreas dataset](https://doi.org/10.1016/j.cmet.2018.11.014).  
Contains single cell data and metadata for three images from three different donors.
- The `counts` **assay** represents raw IMC counts (spillover-compensated).
- **Row names** represent the target proteins and **column names** represent unique cell ids.
- **colData** contain cell-specific metadata:
  - `ImageName`, `ImageNumber`: image name and image number.
  - `CellNumber`, `id`: cell number and cell id.
  - `Pos_X`, `Pos_Y`: cell position on the image.
  - `ParentIslet`: id of the islet the cell belongs to (0 if the cell is not located inside an islet).
  - `ClosestIset`: id of the islet located closest to the cell.
  - `Area`: cell area.
  - `NbNeighbours`: number of neighbours (with a pixel expansion of 4).
  - `width`, `height`: image size.
  - `CellCat`, `CellType`: Cell category (e.g. islet, exocrine, immune) and cell type.
  - `case`, `Age`, `Gender`: donor id, donor age and donor gender.
  - `stage`: donor T1D stage (non-diabetic, recent-onset or long-duration).
  - `part`, `slide`: pancreas region (tail, body or head) and slide name.
  - `group`: matched donor group (donors stratified in three groups matched by BMI, gender age and ethnicity)
- **rowData** contain target-specific metadata:
  - `channel`: channel number (from 1 to 38).
  - `metal`: metal tag (format: "In113").
  - `target`: short name of the target protein.

### survival_example_data.csv

The survival example data set contains overall and disease free 
survival information and the pathology grading of 285 breast cancer 
patients.The pathology grades split breast cancer patients into 3 
groups associated with distinct outcomes.

Column description:  
PID: patient identifier  
grade: tumor pathology grade  
OSmonth: overall suvival from time of diagnosis [months]  
DFSmonth: disease free survival from time of diagnosis [months]  
Patientstatus: life/death status   
"death by primary disease" = death due to disease  
"death" = disease unrelated death  
"alive" = alive at given survival time  
"alive w metastases" = alive with metasatses at given survival time  

