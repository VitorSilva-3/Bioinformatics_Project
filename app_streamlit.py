
import os
import streamlit as st
import pandas as pd

# --- Configuração da página ---
st.set_page_config(
    page_title="Enzimas em microalgas",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Helper para caminhos ---
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))
BASE_DIR = get_base_dir()

# --- Carregar dados com cache ---
@st.cache_data(show_spinner=False)
def load_data():
    return pd.read_csv(os.path.join(BASE_DIR, "results.csv"))

df = load_data()

# --- Usar flag para resetar filtros ---
if "reset_filters" not in st.session_state:
    st.session_state.reset_filters = False

# --- Inicializar filtros apenas no primeiro carregamento ou após reset ---
if "selected_alga" not in st.session_state or st.session_state.reset_filters:
    st.session_state.selected_alga = "Todas"
if "selected_sugar" not in st.session_state or st.session_state.reset_filters:
    st.session_state.selected_sugar = "Todos"
    st.session_state.reset_filters = False  # limpar a flag após uso

# --- Barra lateral: informações e filtros ---
st.sidebar.title("Informações/Filtros")

with st.sidebar.expander("Plataforma", expanded=True):
    st.markdown("")

# Filtros
alga_options = ["Todas"] + sorted(df["Microalga"].unique())
sugar_options = ["Todos"] + sorted(df["Açúcar alvo"].astype(str).unique())

selected_alga = st.sidebar.selectbox("Microalga", alga_options, key="selected_alga")
selected_sugar = st.sidebar.selectbox("Açúcar alvo", sugar_options, key="selected_sugar")

# Botão para limpar filtros
if st.sidebar.button("Limpar filtros"):
    st.session_state.reset_filters = True
    st.rerun()

# --- Aplicar filtros ---
filtered = df.copy()
if selected_alga != "Todas":
    filtered = filtered[filtered["Microalga"] == selected_alga]
if selected_sugar != "Todos":
    filtered = filtered[filtered["Açúcar alvo"] == selected_sugar]

# --- Cabeçalho principal ---
st.title("Enzimas em microalgas")

# Métricas resumidas
c1, c2 = st.columns(2)
c1.metric("Microalgas únicas", filtered["Microalga"].nunique())
c2.metric("Registos mostrados", len(filtered))

# --- Adicionar coluna de status ---
filtered["Estado"] = "🟢"

# --- Formatação de nomes científicos ---
filtered["Microalga"] = filtered["Microalga"].apply(lambda x: f"*{x}*")

# --- Exibir resultados com título correto ---
registo_str = "registo" if len(filtered) == 1 else "registos"
st.subheader(f"Resultados ({len(filtered)} {registo_str})")

# --- Tabela em Markdown para manter itálico e emojis ---
table_md = filtered.to_markdown(index=False)
st.markdown(table_md, unsafe_allow_html=True)