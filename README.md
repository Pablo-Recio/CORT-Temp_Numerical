# CORT-Temp_Numerical

#### Research paper

This repository contains the final code, data, and figures used in the following manuscript:

Quantity discrimination, decision-making, and the role of early-life conditions in a lizard  

Supplementary Materials: The supplementary materials associated with this paper is integrated within the ms.docx or ms.qmd files. 

#### How to use this repository?

Users can download a zip file of the entire repository by clicking on the green code tab at the top of the page and then clicking Download ZIP. Users who already have a GitHub account can fork the repository.

The key file in this repository is the 📄 ms.qmd. This file can be rendered in R with Quarto to reproduce the entire paper. Code chunks within the file provide the code used to reproduce figures and analyses along with supporting statements within the text. Note that inline code chunks use specific objects which are then rendered.

The 📄 ms.qmd file makes use of files within a number of folders that are identified in the code chunks. There are a number of important folders in the repository.

  📂 data folder contains all the raw data used in files. Note that there are different files, but the main one employed in our analyses is lizard_database_number.csv. The file orientation_test.csv contains the results of the orientation test (see Suppl. Mat.). 
  📂 output/figs/ Folder contains all the figures for the paper that are read and included in the paper. 
  📂 R The R folder contains three files, two of them used to clean and process data to prepare it for use in the ms.qmd file. Note that readers do not need to open and run these files, but they are simply here to document the workflow and code used to clean up data to be used. These include:
        📄 data_process.R, which is used to clean the data;
        📄 func.R, which contains all the functions called later in the ms.qmd file. 
  📂 bib The bib folder contains:
        📄 refs.bib the references;
        📄 animal-cognition.csl the journal formatting style file;
        📄 templ.docx a template docx file to format the resulting rendered files.
  📂 ECOLEVOL This folder contains all the material submitted to the journal and all the information from the reviews.

#### Reporting Issues or Asking Questions
If anything is unclear or you require further detail please don't hesitate to lodge an issue.
