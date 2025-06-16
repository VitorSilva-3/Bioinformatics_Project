
import streamlit as st
from Bio import Entrez
import pandas as pd
from data_processing import build_df
from search import (fetch_kegg_enzyme_info, fetch_uniprot_enzyme_info, search_articles_by_organism_with_enzymes)
import time
from data import enzymes

st.set_page_config(
    page_title="Enzymes in Algae",
    layout="wide",
    initial_sidebar_state="expanded"
)
st.title("Enzymes in Algae")

st.sidebar.header("Information & Filters")
with st.sidebar.expander("About this platform", expanded=True):
    st.markdown(
        "This platform integrates data on agro-industrial waste and the enzymes produced by various algae, "
        "organised by species. It identifies which waste materials can serve as substrates for each alga, "
        "indicating those with the greatest growth potential."
    )

Entrez.email = "vtsilva3@gmail.com"

def init_session_state():
    st.session_state.setdefault("data_loaded", False)
    
    st.session_state.setdefault("applied_alga", "All")
    st.session_state.setdefault("applied_sugar", "All")
    st.session_state.setdefault("filters_applied", False)
    
    st.session_state.setdefault("tab2_search_initiated", False)
    st.session_state.setdefault("tab2_combined_articles", {})
    
    st.session_state.setdefault("tab3_taxonomy_loaded", False)
    st.session_state.setdefault("tab3_taxonomy_data", {})
    
    st.session_state.setdefault("tab4_kegg_loaded", False)
    st.session_state.setdefault("tab4_kegg_data", {})
    
    st.session_state.setdefault("tab5_uniprot_loaded", False)
    st.session_state.setdefault("tab5_uniprot_data", {})

init_session_state()

@st.cache_data(ttl=3600)
def load_data() -> pd.DataFrame:
    with st.spinner("Loading data..."):
        time.sleep(0.5)  
        return build_df()

@st.cache_data(ttl=7200) 
def get_all_enzyme_info_cached() -> dict:
    enzyme_info = {}
    for enzyme_name in enzymes.keys():
        info = fetch_kegg_enzyme_info(enzyme_name)
        if info:  
            enzyme_info[enzyme_name] = info
        time.sleep(0.1)
    return enzyme_info

@st.cache_data(ttl=7200)
def get_all_uniprot_enzyme_info_cached() -> dict:
    uniprot_info = {}
    for enzyme_name in enzymes.keys():
        info = fetch_uniprot_enzyme_info(enzyme_name, max_entries=30)
        if info:
            uniprot_info[enzyme_name] = info
        time.sleep(0.5) 
    return uniprot_info

if not st.session_state.data_loaded:
    with st.status("Initializing application...", expanded=True) as status:
        st.write("Loading enzyme data...")
        df_full = load_data()
        st.write("Preparing filter options...")
        time.sleep(0.3) 
        st.write("Application ready!")
        status.update(label="Application initialized successfully!", state="complete", expanded=False)
        st.session_state.data_loaded = True
        st.session_state.df_full = df_full
else:
    df_full = st.session_state.df_full

if df_full.empty:
    st.error("No records available for filtering.")
    st.stop()

col_micro = "Algae"
alga_options = ["All"] + sorted(df_full[col_micro].unique())
sugar_options = ["All"] + sorted(df_full["Target sugar"].unique())

selected_alga = st.sidebar.selectbox("Algae", alga_options, 
                                   index=alga_options.index(st.session_state.get("applied_alga", "All")))
selected_sugar = st.sidebar.selectbox("Target sugar", sugar_options, 
                                    index=sugar_options.index(st.session_state.get("applied_sugar", "All")))

def reset_filters():
    st.session_state.applied_alga = "All"
    st.session_state.applied_sugar = "All"
    st.session_state.filters_applied = False

col_btn1, col_btn2 = st.sidebar.columns(2)
with col_btn1:
    if st.button("Apply Filters", type="primary"):
        st.session_state.applied_alga = selected_alga
        st.session_state.applied_sugar = selected_sugar
        st.session_state.filters_applied = True
        st.rerun()
with col_btn2:
    st.button("Clear Filters", on_click=reset_filters)

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Data",
    "Articles by Alga and Enzyme (NCBI)",
    "Taxonomy",
    "KEGG Enzyme Database",
    "UniProt Database"
])

with tab1:
    st.subheader("Filtered Data Results")
    
    filtered = df_full.copy()
    if st.session_state.get("applied_alga", "All") != "All":
        filtered = filtered[filtered[col_micro] == st.session_state.applied_alga]
    if st.session_state.get("applied_sugar", "All") != "All":
        filtered = filtered[filtered["Target sugar"] == st.session_state.applied_sugar]
    
    filtered = filtered.sort_values(by=col_micro).reset_index(drop=True)
    filtered[col_micro] = filtered[col_micro].astype(str)
    
    count = len(filtered)
    st.write(f"Results ({count} {'record' if count == 1 else 'records'})")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Unique algae", filtered[col_micro].nunique(),
                 help="Number of different algae species in the results")
    with col2:
        st.metric("Displayed records", count,
                 help="Total number of records matching your filters")
    with col3:
        if count > 0:
            completion_rate = (filtered["Status"].str.contains("🟢").sum() / count * 100)
            st.metric("Confirmed enzymes", f"{completion_rate:.1f}%",
                     help="Percentage of records with confirmed enzyme activity")
    
    if count > 100:
        with st.spinner("Rendering large dataset..."):
            st.dataframe(filtered, use_container_width=True)
    else:
        st.dataframe(filtered, use_container_width=True)
    
    st.markdown(
        """
        Status legend:
        - 🟢 Enzyme confirmed
        - 🟡 Probable enzyme
        - 🔴 Enzyme not confirmed
        """
    )

with tab2:
    st.subheader("Scientific Articles by Algae and Enzyme Combination")
    
    organisms = sorted(df_full["Algae"].str.replace(r"\*", "", regex=True).unique())
    
    col_search, col_clear = st.columns([3, 1])
    with col_search:
        if st.button("Search Combined Articles", key="tab2_search_combined"):
            if not organisms:
                st.warning("No algae species found in the dataset.")
            else:
                st.session_state.tab2_search_initiated = True
                st.info(f"Searching for articles combining {len(organisms)} algae species with all enzyme types. This may take several minutes...")
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                all_combination_results = {}
                
                for i, organism in enumerate(organisms):
                    status_text.text(f"Searching combined articles for {organism}... ({i+1}/{len(organisms)})")
                    progress_bar.progress((i + 1) / len(organisms))
                
                    try:
                        organism_enzyme_results = search_articles_by_organism_with_enzymes(organism, max_articles_per_enzyme=20)
                        all_combination_results[organism] = organism_enzyme_results
                        time.sleep(1)  
                    except Exception as e:
                        st.error(f"Error searching combinations for {organism}: {str(e)}")
                        all_combination_results[organism] = {}
                
                progress_bar.empty()
                status_text.empty()
                
                st.session_state.tab2_combined_articles = all_combination_results
    
    with col_clear:
        if st.button("Clear Results", key="tab2_clear"):
            st.session_state.tab2_search_initiated = False
            st.session_state.tab2_combined_articles = {}
            st.rerun()
    
    if st.session_state.tab2_search_initiated and st.session_state.tab2_combined_articles:
        all_combination_results = st.session_state.tab2_combined_articles
        total_articles = 0
        combinations_with_articles = 0
        for organism_results in all_combination_results.values():
            for enzyme_articles in organism_results.values():
                total_articles += len(enzyme_articles)
                if enzyme_articles:
                    combinations_with_articles += 1
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Algae Species", len(organisms))
        with col2:
            st.metric("Total Articles Found", total_articles)
        with col3:
            st.metric("Successful Combinations", combinations_with_articles)
        
        for organism, organism_results in all_combination_results.items():
            organism_article_count = sum(len(articles) for articles in organism_results.values())
            
            with st.expander(f"{organism} - {organism_article_count} combined articles found", 
                            expanded=len(organisms) == 1):
                if organism_article_count == 0:
                    st.warning("No articles found combining this algae species with any enzymes.")
                    continue
                
                for enzyme_name, articles in organism_results.items():
                    if not articles:
                        continue
                    
                    ec_number = enzymes.get(enzyme_name, "Unknown EC")
                    st.markdown(f"### {enzyme_name.title()} (EC {ec_number}) - {len(articles)} articles")
                    
                    for idx, article in enumerate(articles, 1):
                        with st.container():
                            st.markdown(f"{idx}. {article['title']}")
                            
                            authors_str = ", ".join(article['authors'][:3])
                            if len(article['authors']) > 3:
                                authors_str += f" and {len(article['authors']) - 3} others"
                            st.markdown(f"Authors: {authors_str}")
                            
                            st.markdown(f"Journal: {article['journal']} ({article['year']})")
                            st.markdown(f"Search focus: {organism} + {enzyme_name} (EC {ec_number})")
                            
                            col_a, col_b = st.columns(2)
                            with col_a:
                                if article['pmid'] != 'N/A':
                                    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{article['pmid']}"
                                    st.markdown(f"[View on PubMed (ID: {article['pmid']})]({pubmed_url})")
                            with col_b:
                                if article['doi'] != 'No DOI' and 'doi' in article['doi'].lower():
                                    st.markdown(f"DOI: {article['doi']}")
                            
                            if article['abstract']:
                                st.markdown(f"*Abstract:* {article['abstract']}")
                            
                            st.divider()
                    
                    st.markdown("---")

with tab3:
    st.subheader("Taxonomy Details")
    
    organisms = sorted(df_full[col_micro].str.replace(r"\*", "", regex=True).unique())
    
    col_load, col_clear = st.columns([3, 1])
    with col_load:
        if st.button("Load Taxonomy Data", key="tab3_load_taxonomy"):
            st.session_state.tab3_taxonomy_loaded = True
            taxonomy_data = {}
        
            if len(organisms) > 1:
                taxonomy_progress = st.progress(0)
                taxonomy_status = st.empty()
                taxonomy_status.text(f"Loading taxonomy data for {len(organisms)} organisms...")
        
            for i, organism in enumerate(organisms):
                if len(organisms) > 1:
                    taxonomy_status.text(f"Loading taxonomy for {organism}... ({i+1}/{len(organisms)})")
                    taxonomy_progress.progress((i + 1) / len(organisms))
            
                try:
                    with st.spinner(f"Fetching taxonomy data...") if len(organisms) <= 5 else st.empty():
                        handle = Entrez.esearch(db="taxonomy", term=f"{organism}[Scientific Name]")
                        ids = Entrez.read(handle)["IdList"]
                        if ids:
                            tax_handle = Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml")
                            record = Entrez.read(tax_handle)[0]
                            lineage = record.get("Lineage", "")
                        else:
                            lineage = ""
                        taxonomy_data[organism] = lineage
                            
                except Exception as e:
                    taxonomy_data[organism] = f"Error: {str(e)}"
            
            if len(organisms) > 1:
                taxonomy_progress.empty()
                taxonomy_status.empty()

            st.session_state.tab3_taxonomy_data = taxonomy_data

            successful_lookups = sum(1 for lineage in taxonomy_data.values() 
                            if lineage and not lineage.startswith("Error:"))
            st.success(f"Taxonomy lookup completed successfully for {successful_lookups} out of {len(organisms)} organisms")
    
    with col_clear:
        if st.button("Clear Taxonomy", key="tab3_clear"):
            st.session_state.tab3_taxonomy_loaded = False
            st.session_state.tab3_taxonomy_data = {}
            st.rerun()
    
    if st.session_state.tab3_taxonomy_loaded and st.session_state.tab3_taxonomy_data:
        for organism, lineage in st.session_state.tab3_taxonomy_data.items():
            with st.expander(f"{organism}"):
                if lineage and not lineage.startswith("Error:"):
                    st.success("Taxonomy data retrieved successfully")
                    for taxon in lineage.split("; "):
                        st.markdown(f"- {taxon}")
                elif lineage.startswith("Error:"):
                    st.error("Failed to retrieve taxonomy data")
                    st.caption(lineage)
                else:
                    st.warning("No taxonomy data available for this organism")

with tab4:
    st.subheader("KEGG Enzyme Database Information")
    
    col_load, col_clear = st.columns([3, 1])
    with col_load:
        if st.button("Load KEGG Data", key="tab4_load_kegg"):
            st.session_state.tab4_kegg_loaded = True
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            enzyme_info = {}
            total_enzymes = len(enzymes)
            
            for i, enzyme_name in enumerate(enzymes.keys()):
                status_text.text(f"Fetching KEGG data for {enzyme_name}...")
                progress_bar.progress((i + 1) / total_enzymes)
                
                info = fetch_kegg_enzyme_info(enzyme_name)
                if info:  
                    enzyme_info[enzyme_name] = info
                
                time.sleep(0.1)
            
            progress_bar.empty()
            status_text.empty()
            
            st.session_state.tab4_kegg_data = enzyme_info
    
    with col_clear:
        if st.button("Clear KEGG Data", key="tab4_clear"):
            st.session_state.tab4_kegg_loaded = False
            st.session_state.tab4_kegg_data = {}
            st.rerun()
    
    if st.session_state.tab4_kegg_loaded:
        all_enzyme_info = st.session_state.tab4_kegg_data
        
        if not all_enzyme_info:
            st.warning("No enzyme information could be retrieved from KEGG database.")
        else:
            st.success(f"Successfully retrieved information for {len(all_enzyme_info)} enzymes from KEGG")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Enzymes", len(enzymes))
            with col2:
                st.metric("Successfully Retrieved", len(all_enzyme_info))
            with col3:
                unique_ec_numbers = len(set(enzymes.values()))
                st.metric("Unique EC Numbers", unique_ec_numbers)
            
            for enzyme_name, info in all_enzyme_info.items():
                with st.expander(f"{enzyme_name.title()} (EC {info.get('ec_number', 'N/A')})"):
                    info_col1, info_col2 = st.columns([2, 1])
                    
                    with info_col1:
                        if info.get('name'):
                            st.markdown(f"Official Name: {info['name']}")
                        
                        if info.get('definition'):
                            st.markdown(f"Definition: {info['definition']}")
                        
                        if info.get('reaction'):
                            st.markdown(f"Reaction: {info['reaction']}")
                    
                    with info_col2:
                        ec_number = info.get('ec_number', '')
                        if ec_number:
                            st.markdown(f"EC Number: {ec_number}")
                            ec_parts = ec_number.split('.')
                            if len(ec_parts) >= 4:
                                st.markdown(f"Class: {ec_parts[0]} (Main enzyme class)")
                                st.markdown(f"Subclass: {ec_parts[1]}")
                                st.markdown(f"Sub-subclass: {ec_parts[2]}")
                                st.markdown(f"Serial number: {ec_parts[3]}")
                    
                    if info.get('pathways'):
                        st.markdown("Associated Metabolic Pathways:")
                        for pathway in info['pathways']:
                            st.markdown(f"- {pathway['code']}: {pathway['description']}")
                    else:
                        st.markdown("No specific pathways found in KEGG database")
                    
                    if info.get('ec_number'):
                        kegg_url = f"https://www.genome.jp/entry/ec:{info['ec_number']}"
                        st.markdown(f"[View detailed information on KEGG database]({kegg_url})")

with tab5:
    st.subheader("UniProt Database Information")
    
    col_load, col_clear = st.columns([3, 1])
    with col_load:
        if st.button("Load UniProt Data", key="tab5_load_uniprot"):
            st.session_state.tab5_uniprot_loaded = True
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            uniprot_info = {}
            total_enzymes = len(enzymes)
            
            for i, enzyme_name in enumerate(enzymes.keys()):
                status_text.text(f"Fetching UniProt data for {enzyme_name}...")
                progress_bar.progress((i + 1) / total_enzymes)
                
                info = fetch_uniprot_enzyme_info(enzyme_name, max_entries=30)
                if info:
                    uniprot_info[enzyme_name] = info
                
                time.sleep(0.5) 
            
            progress_bar.empty()
            status_text.empty()
            
            st.session_state.tab5_uniprot_data = uniprot_info
    
    with col_clear:
        if st.button("Clear UniProt Data", key="tab5_clear"):
            st.session_state.tab5_uniprot_loaded = False
            st.session_state.tab5_uniprot_data = {}
            st.rerun()
    
    if st.session_state.tab5_uniprot_loaded:
        all_uniprot_info = st.session_state.tab5_uniprot_data
        
        if not all_uniprot_info:
            st.warning("No enzyme information could be retrieved from UniProt database.")
        else:
            enzymes_with_data = [name for name, info in all_uniprot_info.items() if info.get('general_info')]
            st.success(f"Successfully retrieved information for {len(enzymes_with_data)} enzymes from UniProt")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Enzymes Queried", len(enzymes))
            with col2:
                st.metric("Enzymes with UniProt Data", len(enzymes_with_data))
            with col3:
                total_entries_analyzed = sum(info.get('total_entries_analyzed', 0) for info in all_uniprot_info.values())
                st.metric("Total UniProt Entries Analyzed", total_entries_analyzed)
            
            for enzyme_name, info in all_uniprot_info.items():
                general_info = info.get('general_info', {})
                if not general_info:
                    continue
                    
                entries_analyzed = info.get('total_entries_analyzed', 0)
                with st.expander(f"{enzyme_name.title()} (EC {info.get('ec_number', 'N/A')}) - Analyzed {entries_analyzed} entries"):
                    
                    enzyme_col1, enzyme_col2 = st.columns([3, 2])
                    
                    with enzyme_col1:
                        if general_info.get('protein_name') and general_info['protein_name'] != 'Unknown':
                            st.markdown(f"Protein Name: {general_info['protein_name']}")
                        
                        st.markdown("Biochemical Functions:")
                        functions = general_info.get('functions', [])
                        if functions and functions != ['No function information available']:
                            for func in functions[:3]:
                                st.markdown(f"• {func}")
                        else:
                            st.markdown("• No specific function information available")
                        
                        if general_info.get('catalytic_activities'):
                            st.markdown("Catalytic Activities:")
                            for activity in general_info['catalytic_activities']:
                                reaction = activity.get('reaction', 'Unknown reaction')
                                st.markdown(f"• {reaction}")
                    
                    with enzyme_col2:
                        ec_number = info.get('ec_number', '')
                        if ec_number:
                            st.markdown(f"EC Number: {ec_number}")
                            ec_parts = ec_number.split('.')
                            if len(ec_parts) >= 4:
                                st.markdown(f"Enzyme Class: {ec_parts[0]}")
                        
                        if general_info.get('cofactors'):
                            st.markdown("Required Cofactors:")
                            for cofactor in general_info['cofactors']:
                                st.markdown(f"• {cofactor}")
                        
                        if general_info.get('subunit_structure'):
                            st.markdown("Subunit Structure:")
                            for subunit in general_info['subunit_structure']:
                                st.markdown(f"• {subunit}")
                    
                    if general_info.get('pathways'):
                        st.markdown("Associated Metabolic Pathways:")
                        for pathway in general_info['pathways']:
                            st.markdown(f"• {pathway}")
                    
                    if ec_number:
                        uniprot_search_url = f"https://www.uniprot.org/uniprotkb?query=ec%3A{ec_number}"
                        st.markdown(f"[Search this enzyme on UniProt database]({uniprot_search_url})")
            
            no_data_enzymes = [name for name, info in all_uniprot_info.items() if not info.get('general_info')]
            if no_data_enzymes:
                with st.expander(f"Enzymes with no UniProt biochemical data ({len(no_data_enzymes)})"):
                    st.markdown("The following enzymes did not return sufficient biochemical information from UniProt database:")
                    for enzyme_name in no_data_enzymes:
                        ec_number = enzymes.get(enzyme_name, 'Unknown')
                        st.markdown(f"• {enzyme_name.title()} (EC {ec_number})")