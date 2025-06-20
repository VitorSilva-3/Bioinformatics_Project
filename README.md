# Bioinformatics Project

Repository created as part of a Bioinformatics Project.

Created by [Vítor Silva](https://github.com/VitorSilva-3) PG55538  

Access the live platform at: [https://enzymesalgae.streamlit.app/](https://enzymesalgae.streamlit.app/)

---

## Enzymes in Algae Platform

### Files overview

- **Main results/**  
  A folder with images of the platform, containing the main results.
  
- **Presentation.pdf**  
  Used to show the main results of the project.
  
- **Vítor_Silva_PG55538.pdf**  
  Article of the project

- **app.py**  
  Main Streamlit interface with 5 tabs for data visualisation, article search, taxonomy info, and enzyme databases (KEGG/UniProt).

- **data.py**  
  Contains enzyme configurations: EC numbers, target sugars, search queries, and algae taxonomic groups.

- **data_processing.py**  
  Fetches protein data from NCBI, processes GenBank records, and builds the main dataset with enzyme classification.

- **search.py**  
  Integrates external APIs for PubMed articles, KEGG enzyme pathways, and UniProt biochemical data.

- **utils.py**  
  Helper functions for text normalization and identifying hypothetical proteins.

- **requirements.txt**  
  Python packages that were used in this project.
