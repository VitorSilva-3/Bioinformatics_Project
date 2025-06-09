
import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import Entrez, SeqIO, Medline
import re
import time
import urllib.error
from data import ec_manual, sugar_manual, taxa, queries
from collections import defaultdict

def normalize_text(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())

normalised_keys = []
for full_key, ec_number in ec_manual.items():
    subkeys = [sub.strip() for sub in full_key.split("/")]
    normalised_subkeys = [normalize_text(sub) for sub in subkeys]
    normalised_keys.append((full_key, normalised_subkeys, ec_number))

future_terms = ["hypothetical", "similar", "putative", "uncharacterized", "probable"]

st.set_page_config(
    page_title="Enzymes in Algae",
    layout="wide",
    initial_sidebar_state="expanded"
)
st.title("Enzymes in Algae")

st.sidebar.header("Information & Filters")
with st.sidebar.expander("About this platform", expanded=True):
    st.markdown(
        "This platform integrates data on agro-industrial waste and the enzymes produced by various algae, organised by species. "
        "It identifies which waste materials can serve as substrates for each alga, indicating those with the greatest growth potential."
    )

Entrez.email = "vtsilva3@gmail.com"

def fetch_ids():
    taxa_query = " OR ".join(f'"{t}"[Organism]' for t in taxa)
    all_ids = set()
    for q in queries:
        full_query = f"{q} AND ({taxa_query})"
        handle = Entrez.esearch(db="protein", term=full_query, retmax=10000)
        all_ids.update(Entrez.read(handle)["IdList"])
    return list(all_ids)

def search_articles_by_organism(organism, max_articles=8):
    articles = []
    clean_organism = re.sub(r'[^\w\s]', '', organism).strip()
    search_strategies = [
        f'("{clean_organism}"[Organism] OR "{clean_organism}"[Title/Abstract]) AND (metabolism[MeSH Terms] OR enzyme*[Title/Abstract] OR biochemical[Title/Abstract])',
        f'"{clean_organism}"[Title/Abstract] AND (biotechnology[Title/Abstract] OR industrial[Title/Abstract] OR biofuel*[Title/Abstract])',
        f'"{clean_organism}"[Organism] AND (cultivation[Title/Abstract] OR growth[Title/Abstract] OR production[Title/Abstract])',
        f'"{clean_organism}"[Title/Abstract]'
    ]
    
    for strategy_idx, query in enumerate(search_strategies):
        try:
            handle = Entrez.esearch(
                db="pubmed", 
                term=query, 
                retmax=max_articles//len(search_strategies) + 2,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            if search_results["IdList"]:
                pmids = search_results["IdList"]
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(pmids),
                    rettype="medline",
                    retmode="text"
                )
                
                medline_records = Medline.parse(handle)
                
                for record in medline_records:
                    pmid = record.get('PMID', 'N/A')
                    if not any(art['pmid'] == pmid for art in articles):
                        article_info = {
                            'pmid': pmid,
                            'title': record.get('TI', 'No title available'),
                            'authors': record.get('AU', ['No authors listed']),
                            'journal': record.get('JT', 'Unknown journal'),
                            'year': record.get('DP', 'Unknown year'),
                            'abstract': record.get('AB', 'No abstract available'),
                            'doi': record.get('LID', ['No DOI'])[0] if record.get('LID') else 'No DOI',
                            'search_focus': 'organism',
                            'search_strategy': strategy_idx + 1
                        }
                        articles.append(article_info)
                
                handle.close()
                
                if len(articles) >= max_articles:
                    break
                    
        except Exception as e:
            print(f"Error searching for alga (strategy {strategy_idx + 1}): {e}")
            continue
    
    return articles[:max_articles]

def search_articles_by_enzyme_type(enzyme_name, ec_number=None, max_articles=8):
    articles = []
    clean_enzyme = re.sub(r'[^\w\s]', '', enzyme_name).strip()
    search_strategies = []
    
    if ec_number and ec_number != "unknown":
        search_strategies.extend([
            f'"{ec_number}"[Title/Abstract]',
            f'("{ec_number}"[Title/Abstract] OR "{clean_enzyme}"[Title/Abstract]) AND enzyme*[Title/Abstract]'
        ])
    
    search_strategies.extend([
        f'"{clean_enzyme}"[Title/Abstract] AND (industrial[Title/Abstract] OR biotechnology[Title/Abstract] OR application*[Title/Abstract])',
        f'"{clean_enzyme}"[Title/Abstract] AND (characterization[Title/Abstract] OR purification[Title/Abstract] OR activity[Title/Abstract])',
        f'"{clean_enzyme}"[Title/Abstract] AND enzyme*[Title/Abstract]',
        f'"{clean_enzyme}"[Title/Abstract]'
    ])
    
    for strategy_idx, query in enumerate(search_strategies):
        try:
            handle = Entrez.esearch(
                db="pubmed", 
                term=query, 
                retmax=max_articles//len(search_strategies) + 2,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            if search_results["IdList"]:
                pmids = search_results["IdList"]
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(pmids),
                    rettype="medline",
                    retmode="text"
                )
                
                medline_records = Medline.parse(handle)
                
                for record in medline_records:
                    pmid = record.get('PMID', 'N/A')
                    if not any(art['pmid'] == pmid for art in articles):
                        article_info = {
                            'pmid': pmid,
                            'title': record.get('TI', 'No title available'),
                            'authors': record.get('AU', ['No authors listed']),
                            'journal': record.get('JT', 'Unknown journal'),
                            'year': record.get('DP', 'Unknown year'),
                            'abstract': record.get('AB', 'No abstract available'),
                            'doi': record.get('LID', ['No DOI'])[0] if record.get('LID') else 'No DOI',
                            'search_focus': 'enzyme',
                            'search_strategy': strategy_idx + 1,
                            'ec_number': ec_number
                        }
                        articles.append(article_info)
                
                handle.close()
                
                if len(articles) >= max_articles:
                    break
                    
        except Exception as e:
            print(f"Error searching for enzyme (strategy {strategy_idx + 1}): {e}")
            continue
    
    return articles[:max_articles]

def display_articles(articles, focus_type):
    if not articles:
        st.info(f"No relevant articles found for this {focus_type}.")
        return
    
    st.markdown(f"**Found {len(articles)} relevant articles:**")
    
    for article in articles:
        with st.container():
            st.markdown(f"### {article['title']}")
            
            info_parts = []
            
            if article['authors'] and article['authors'] != ['No authors listed']:
                authors_str = ", ".join(article['authors'][:3])
                if len(article['authors']) > 3:
                    authors_str += " et al."
                info_parts.append(f"**Authors:** {authors_str}")
            
            journal_year = f"**Journal:** {article['journal']}"
            if article['year'] != 'Unknown year':
                journal_year += f" ({article['year']})"
            info_parts.append(journal_year)
            
            if article['pmid'] != 'N/A':
                info_parts.append(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{article['pmid']}/)")
            
            st.markdown(" • ".join(info_parts))
            
            if article['abstract'] != 'No abstract available':
                st.markdown("**Abstract:**")
                st.markdown(f"> {article['abstract']}")
            
            st.markdown("---")

def display_organism_articles_tab(filtered_df, show_articles_flag):
    if show_articles_flag:
        st.subheader("Scientific articles by algae")
        st.markdown("*Articles focused on algae.*")
        
        unique_organisms = sorted(filtered_df['Algae'].unique())
        
        for organism in unique_organisms:
            with st.expander(f"**{organism}**", expanded=False):
                with st.spinner(f"Searching articles for {organism}..."):
                    articles = load_organism_articles(organism)
                
                display_articles(articles, "organism")
    else:
        st.info("Enable *Show scientific articles* in the sidebar to view relevant articles.")

def display_enzyme_articles_tab(filtered_df, show_articles_flag):
    if show_articles_flag:
        st.subheader("Scientific articles by enzymes")
        st.markdown("*Articles focused on enzymes.*")
        
        confirmed_enzymes_df = filtered_df[filtered_df['Status'] == '🟢']
        
        if confirmed_enzymes_df.empty:
            st.info("No confirmed enzymes (🟢 status) found with the current filters.")
            return
        
        enzyme_groups = confirmed_enzymes_df.groupby(['Enzyme', 'EC number']).size().reset_index(name='count')
        enzyme_groups = enzyme_groups.sort_values('count', ascending=False)
        
        for _, enzyme_row in enzyme_groups.iterrows():
            enzyme_name = enzyme_row['Enzyme']
            ec_number = enzyme_row['EC number']
            
            expander_title = f"**{enzyme_name}**"
            if ec_number and ec_number != "unknown":
                expander_title += f" (EC: {ec_number})"
            
            with st.expander(expander_title, expanded=False):
                with st.spinner(f"Searching articles for {enzyme_name}..."):
                    articles = load_enzyme_articles(enzyme_name, ec_number)
                
                display_articles(articles, "enzyme")
    else:
        st.info("Enable *Show scientific articles* in the sidebar to view relevant articles.")

def fetch_records(id_list, batch_size=500, max_retries=3):
    records = []
    if not id_list:
        return records

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    total_batches = (len(id_list) + batch_size - 1) // batch_size

    for batch_idx, batch in enumerate(chunks(id_list, batch_size), start=1):
        attempt = 0
        while attempt < max_retries:
            try:
                ids_str = ",".join(batch)
                handle = Entrez.efetch(
                    db="protein",
                    id=ids_str,
                    rettype="gb",
                    retmode="text"
                )
                batch_records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                records.extend(batch_records)
                break

            except urllib.error.HTTPError as e:
                attempt += 1
                if e.code == 500 and attempt < max_retries:
                    wait_time = 2 ** attempt
                    print(f"[WARNING] Batch {batch_idx}/{total_batches} returned HTTP 500. Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"[ERROR] Batch {batch_idx} failed with HTTP {e.code}: {e}")
                    raise

            except Exception as e:
                attempt += 1
                if attempt < max_retries:
                    wait_time = 2 ** attempt
                    print(f"[WARNING] Unexpected error in batch {batch_idx}. Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"[ERROR] Critical failure in batch {batch_idx}: {e}")
                    raise

        time.sleep(0.5)

    return records

def build_df():
    rows = []
    ids = fetch_ids()
    records = fetch_records(ids)

    for rec in records:
        organism = rec.annotations.get("organism", "unknown")
        cds = next((f for f in rec.features if f.type == "CDS"), None)
        product_raw = cds.qualifiers.get("product", [rec.description])[0]
        product = re.split(r"\s*\[", product_raw)[0].strip()
        lower_prod = product.lower()

        has_term = any(t in lower_prod for t in future_terms)

        matched_ec = None
        matched_sugar = None

        if not has_term:
            normalised_product = normalize_text(product)
            for full_key, subkeys_norm, ec_number in normalised_keys:
                if all(sub in normalised_product for sub in subkeys_norm):
                    matched_ec = ec_number
                    matched_sugar = sugar_manual.get(full_key, "unknown")
                    break

        if has_term:
            if "uncharacterized" in lower_prod:
                status = "🔴"
            else:
                status = "🟡"
        elif matched_ec is not None:
            status = "🟢"
        else:
            continue

        if has_term and matched_ec is None:
            matched_ec = "unknown"
            matched_sugar = "unknown"

        rows.append({
            "Algae":         organism,
            "Enzyme":        product,
            "EC number":     matched_ec,
            "Target sugar":  matched_sugar,
            "Description":   product_raw,
            "Status":        status
        })

    df = pd.DataFrame(rows)
    return df

@st.cache_data(ttl=3600)
def load_data():
    return build_df()

@st.cache_data(ttl=7200)
def load_organism_articles(organism):
    return search_articles_by_organism(organism, max_articles=8)

@st.cache_data(ttl=7200)
def load_enzyme_articles(enzyme_name, ec_number):
    return search_articles_by_enzyme_type(enzyme_name, ec_number, max_articles=8)

st.session_state.setdefault("show_results", False)
st.session_state.setdefault("sel_alga", "All")
st.session_state.setdefault("sel_sugar", "All")
st.session_state.setdefault("show_articles", False)

df_temp = load_data()
if df_temp.empty:
    st.warning("No records available for filtering.")
    st.stop()

col_micro = "Algae"
alga_options  = ["All"] + sorted(df_temp[col_micro].unique())
sugar_options = ["All"] + sorted(df_temp["Target sugar"].unique())
st.sidebar.selectbox("Algae", alga_options, key="sel_alga")
st.sidebar.selectbox("Target sugar", sugar_options, key="sel_sugar")

st.sidebar.checkbox("Show scientific articles", key="show_articles", 
                   help="Enable this to search for relevant scientific articles organized by algae and enzymes.")

if st.sidebar.button("Show results"):
    st.session_state.show_results = True

def reset_filters():
    st.session_state.sel_alga     = "All"
    st.session_state.sel_sugar    = "All"
    st.session_state.show_results = False
    st.session_state.show_articles = False

st.sidebar.button("Clear filters", on_click=reset_filters)

if st.session_state.show_results:
    df = load_data()
    if df.empty:
        st.warning("No records found for the specified query.")
        st.stop()

    filtered = df.copy()
    if st.session_state.sel_alga != "All":
        filtered = filtered[filtered[col_micro] == st.session_state.sel_alga]
    if st.session_state.sel_sugar != "All":
        filtered = filtered[filtered["Target sugar"] == st.session_state.sel_sugar]

    filtered = filtered.sort_values(by=col_micro, ascending=True).reset_index(drop=True)
    filtered[col_micro] = filtered[col_micro].apply(lambda x: f"{x}")

    tab1, tab2, tab3, tab4 = st.tabs(["Data", "Articles by Organism", "Articles by Enzyme", "Taxonomy"])

    with tab1:
        count = len(filtered)
        rec_str = "record" if count == 1 else "records"
        c1, c2 = st.columns(2)
        c1.metric("Unique algae", filtered[col_micro].nunique())
        c2.metric("Displayed records", count)
        st.subheader(f"Results ({count} {rec_str})")
        st.dataframe(filtered, use_container_width=True)
        st.markdown(
            """
            Status:

            - 🟢 Enzyme confirmed
            - 🟡 Probable enzyme
            - 🔴 Enzyme not confirmed
            """
        )

    with tab2:
        display_organism_articles_tab(filtered, st.session_state.show_articles)

    with tab3:
        display_enzyme_articles_tab(filtered, st.session_state.show_articles)

    with tab4:
        st.subheader("Taxonomy details")
        for org in sorted(filtered[col_micro].str.replace(r"\*", "", regex=True).unique()):
            try:
                handle = Entrez.esearch(db="taxonomy", term=f"{org}[Scientific Name]")
                ids = Entrez.read(handle)["IdList"]
                lineage = ""
                if ids:
                    tax_handle = Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml")
                    tax_record = Entrez.read(tax_handle)[0]
                    lineage = tax_record.get("Lineage", "")
                with st.expander(org, expanded=False):
                    if lineage:
                        for taxon in lineage.split("; "):
                            st.markdown(f"- {taxon}")
                    else:
                        st.write("No taxonomy data available.")
            except Exception:
                with st.expander(org, expanded=False):
                    st.write("No taxonomy data available.")
else:
    st.info("Select your filters and click *Show results* to view the data.")