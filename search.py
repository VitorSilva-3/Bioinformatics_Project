
import re
import requests
from Bio import Entrez, Medline
from Bio.KEGG import REST
from data import enzymes

Entrez.email = "vtsilva3@gmail.com"

def search_articles_by_organism_with_enzymes(organism: str, max_articles_per_enzyme: int = 20) -> dict:
    results = {}
    clean_org = re.sub(r"[^\w\s]", "", organism).strip()
    
    for enzyme_name, ec_number in enzymes.items():
        articles = []
        clean_enzyme = re.sub(r"[^\w\s]", "", enzyme_name).strip()
        
        query_parts = [
            f'"{clean_org}"[Organism]',  
            f'("{clean_enzyme}"[Title/Abstract] OR "{ec_number}"[Title/Abstract])'
        ]
        query = f'({query_parts[0]}) AND ({query_parts[1]})'
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_articles_per_enzyme, sort="relevance")
            search_results = Entrez.read(handle)
            pmids = search_results.get("IdList", [])
            handle.close()
            
            if pmids:
                handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="medline", retmode="text")
                
                for rec in Medline.parse(handle):
                    pmid = rec.get('PMID', 'N/A')
                    
                    if any(a['pmid'] == pmid for a in articles):
                        continue
                    
                    articles.append({
                        'pmid': pmid,
                        'title': rec.get('TI', 'No title available'),
                        'authors': rec.get('AU', []),
                        'journal': rec.get('JT', 'Unknown journal'),
                        'year': rec.get('DP', 'Unknown year'),
                        'abstract': rec.get('AB', ''),
                        'doi': rec.get('LID', ['No DOI'])[0] if rec.get('LID') else 'No DOI',
                        'search_focus': 'organism_enzyme',
                        'organism': organism,
                        'enzyme_name': enzyme_name,
                        'ec_number': ec_number
                    })
                
                handle.close()
                print(f"[search_articles_by_organism_with_enzymes] Found {len(articles)} articles for {organism} + {enzyme_name}")
            else:
                print(f"[search_articles_by_organism_with_enzymes] No articles found for {organism} + {enzyme_name}")
                
        except Exception as e:
            print(f"[search_articles_by_organism_with_enzymes] Error searching {organism} + {enzyme_name}: {e}")
        
        results[enzyme_name] = articles[:max_articles_per_enzyme]
    
    return results


def fetch_kegg_enzyme_info(enzyme_name: str) -> dict:
    if enzyme_name not in enzymes:
        print(f"[fetch_kegg_enzyme_info] Enzyme '{enzyme_name}' not in allowed list.")
        return {}

    ec_number = enzymes[enzyme_name]
    entry_id = f"ec:{ec_number}"
    try:
        raw = REST.kegg_get(entry_id).read()
    except Exception as e:
        print(f"[fetch_kegg_enzyme_info] Error fetching KEGG entry {entry_id}: {e}")
        return {}

    info = {
        'ec_number': ec_number,
        'name': None,
        'definition': None,
        'reaction': None,
        'pathways': []
    }

    for line in raw.splitlines():
        if line.startswith("NAME"):
            info['name'] = line.split("NAME")[1].strip().rstrip(";")
        elif line.startswith("DEFINITION"):
            info['definition'] = line.split("DEFINITION")[1].strip()
        elif line.startswith("REACTION"):
            info['reaction'] = line.split("REACTION")[1].strip()
        elif line.startswith("PATHWAY"):
            parts = line.split()
            code = parts[1]
            desc = " ".join(parts[2:])
            info['pathways'].append({'code': code, 'description': desc})

    return info


def fetch_uniprot_enzyme_info(enzyme_name: str, max_entries: int = 30) -> dict:
    if enzyme_name not in enzymes:
        print(f"[fetch_uniprot_enzyme_info] Enzyme '{enzyme_name}' not in allowed list.")
        return {}

    ec_number = enzymes[enzyme_name]
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    params = {
        'query': f'ec:{ec_number}',
        'format': 'json',
        'size': max_entries,
        'fields': 'accession,protein_name,cc_function,cc_catalytic_activity,cc_pathway,cc_cofactor,cc_subunit'
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        if not data.get('results'):
            print(f"[fetch_uniprot_enzyme_info] No UniProt entries found for {enzyme_name} (EC {ec_number})")
            return {'enzyme_name': enzyme_name, 'ec_number': ec_number, 'general_info': {}}
        
        function_descriptions = []
        catalytic_activities = []
        pathways = []
        cofactors = []
        subunit_info = []
        
        for entry in data['results']:
            for comment in entry.get('comments', []):
                comment_type = comment.get('commentType')
                
                if comment_type == 'FUNCTION':
                    for text in comment.get('texts', []):
                        if text.get('value') and text['value'] not in function_descriptions:
                            function_descriptions.append(text['value'])
                
                elif comment_type == 'CATALYTIC ACTIVITY':
                    reaction_info = comment.get('reaction', {})
                    if reaction_info:
                        catalytic_activities.append({
                            'reaction': reaction_info.get('name', 'Unknown reaction'),
                            'ec_number': reaction_info.get('ecNumber', ec_number)
                        })
                
                elif comment_type == 'PATHWAY':
                    for text in comment.get('texts', []):
                        if text.get('value') and text['value'] not in pathways:
                            pathways.append(text['value'])
                
                elif comment_type == 'COFACTOR':
                    for text in comment.get('texts', []):
                        if text.get('value') and text['value'] not in cofactors:
                            cofactors.append(text['value'])
                
                elif comment_type == 'SUBUNIT':
                    for text in comment.get('texts', []):
                        if text.get('value') and text['value'] not in subunit_info:
                            subunit_info.append(text['value'])
        
        protein_name = 'Unknown'
        if data['results']:
            first_entry = data['results'][0]
            protein_desc = first_entry.get('proteinDescription', {})
            recommended_name = protein_desc.get('recommendedName', {})
            if recommended_name:
                protein_name = recommended_name.get('fullName', {}).get('value', 'Unknown')
        
        result = {
            'enzyme_name': enzyme_name,
            'ec_number': ec_number,
            'general_info': {
                'protein_name': protein_name,
                'functions': function_descriptions[:3] if function_descriptions else ['No function information available'],
                'catalytic_activities': catalytic_activities[:3] if catalytic_activities else [],
                'pathways': pathways[:3] if pathways else [],
                'cofactors': cofactors[:3] if cofactors else [],
                'subunit_structure': subunit_info[:2] if subunit_info else []
            },
            'total_entries_analyzed': len(data['results'])
        }
        
        print(f"[fetch_uniprot_enzyme_info] Analyzed {len(data['results'])} entries for {enzyme_name}")
        return result
        
    except Exception as e:
        print(f"[fetch_uniprot_enzyme_info] Error fetching UniProt data for {enzyme_name}: {e}")
        return {}