
import time
import urllib.error
import re
import pandas as pd
from Bio import Entrez, SeqIO
from data import ec_manual, sugar_manual, taxa, queries
from utils import normalize_text, future_terms

normalised_keys = []
for full_key, ec_number in ec_manual.items():
    subkeys = [sub.strip() for sub in full_key.split("/")]
    normalised_subkeys = [normalize_text(sub) for sub in subkeys]
    normalised_keys.append((full_key, normalised_subkeys, ec_number))


def fetch_ids() -> list[str]:
    taxa_query = " OR ".join(f'"{t}"[Organism]' for t in taxa)
    all_ids = set()
    for q in queries:
        full_query = f"{q} AND ({taxa_query})"
        handle = Entrez.esearch(db="protein", term=full_query, retmax=10000)
        all_ids.update(Entrez.read(handle)["IdList"])
    return list(all_ids)


def fetch_records(id_list: list[str], batch_size=500, max_retries=3) -> list:
    records = []
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i+n]

    total = len(id_list)
    for batch in chunks(id_list, batch_size):
        attempt = 0
        while attempt < max_retries:
            try:
                handle = Entrez.efetch(db="protein", id=','.join(batch), rettype="gb", retmode="text")
                records.extend(list(SeqIO.parse(handle, "genbank")))
                handle.close()
                break
            except urllib.error.HTTPError as e:
                attempt += 1
                if e.code == 500 and attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                else:
                    raise
            except Exception:
                attempt += 1
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                else:
                    raise
        time.sleep(0.5)
    return records


def build_df() -> pd.DataFrame:
    rows = []
    ids = fetch_ids()
    records = fetch_records(ids)

    for rec in records:
        organism = rec.annotations.get("organism", "unknown")
        cds = next((f for f in rec.features if f.type == "CDS"), None)
        try:
            product_raw = cds.qualifiers.get("product", [rec.description])[0]
        except Exception:
            product_raw = rec.description
        product = re.split(r"\s*\[", product_raw)[0].strip()
        lower = product.lower()
        has_term = any(t in lower for t in future_terms)

        matched_ec = None
        matched_sugar = None
        if not has_term:
            norm = normalize_text(product)
            for full_key, subkeys, ec in normalised_keys:
                if all(sub in norm for sub in subkeys):
                    matched_ec = ec
                    matched_sugar = sugar_manual.get(full_key, "unknown")
                    break

        if has_term:
            status = "🔴" if "uncharacterized" in lower else "🟡"
        elif matched_ec:
            status = "🟢"
        else:
            continue

        if has_term and not matched_ec:
            matched_ec, matched_sugar = "unknown", "unknown"

        rows.append({
            "Algae": organism,
            "Enzyme": product,
            "EC number": matched_ec,
            "Target sugar": matched_sugar,
            "Description": product_raw,
            "Status": status
        })
    return pd.DataFrame(rows)