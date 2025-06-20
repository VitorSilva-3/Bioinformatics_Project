# Bioinformatics Project

Repository created as part of a Bioinformatics Project.

Created by [Vítor Silva](https://github.com/VitorSilva-3) PG55538

---

## Enzymes in Algae Platform

### Scripts Overview

- **app.py**  
  Main Streamlit interface with 5 tabs for data visualization, article search, taxonomy info, and enzyme databases (KEGG/UniProt)

- **data.py**  
  Contains enzyme configurations: EC numbers, target sugars, search queries, and algae taxonomic groups

- **data_processing.py**  
  Fetches protein data from NCBI, processes GenBank records, and builds the main dataset with enzyme classification

- **search.py**  
  Integrates external APIs for PubMed articles, KEGG enzyme pathways, and UniProt biochemical data

- **utils.py**  
  Helper functions for text normalization and identifying hypothetical proteins

- **requirements.txt**  
  Python packages that were used in this project
