
import re

future_terms = ["hypothetical", "similar", "putative", "uncharacterized", "probable", "partial"]

def normalize_text(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())