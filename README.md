# Bioinformatics Project: Enzymes in Algae

This repository contains the code and resources for a Bioinformatics Project.

---

## Project Overview

* **Developer:** [Vítor Silva](https://github.com/VitorSilva-3) (PG55538)  
* **Live Application:** [https://enzymesalgae.streamlit.app/](https://enzymesalgae.streamlit.app/)  

---

## Objectives

- **Link waste derived sugars to microalgal enzymes**  
  Develop a bioinformatics platform that maps specific sugars from agro‑industrial waste to microalgal strains encoding the enzymes necessary to hydrolyse them into simple sugars (e.g. glucose).

- **Integrate enzyme profiles**  
  Consolidate enzyme data from diverse microalgal species into a single, user‑friendly interface for rapid querying and comparison.

- **Match byproducts to algal strains**  
  Use the platform to identify which agro‑industrial by‑products can serve as substrates and pinpoint the microalgae species with the greatest cultivation potential for each.  

---

## The Enzymes in Algae Platform  

This application offers a 5 tab interface for:  
* Data visualisation  
* Article search  
* Taxonomy information  
* Database information (KEGG/UniProt)   

### Repository Files  

* **`Main results/`**: Visual outputs and key results from the platform.
   
* **`Presentation.pdf`**: Project presentation.
  
* **`Vítor_Silva_PG55538.pdf`**: Full project article.
    
* **`app.py`**: Main Streamlit interface.

* **`data.py`**: Enzyme configurations (EC numbers, sugars, queries, taxonomy).
    
* **`data_processing.py`**: Processes NCBI protein data and builds the enzyme dataset.
    
* **`search.py`**: Integrates PubMed, KEGG, and UniProt APIs.
  
* **`utils.py`**: Helper functions for text normalisation and protein identification.
    
* **`requirements.txt`**: Lists all project dependencies.    
