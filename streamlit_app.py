
# streamlit_app.py
import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import seq3, molecular_weight
import matplotlib.pyplot as plt

# Configuración inicial
st.title("Dashboard Bioinformático Interactivo")
st.sidebar.header("Configuración")
analysis_type = st.sidebar.selectbox("Selecciona el tipo de análisis:", ["Análisis de ADN/ARN", "Análisis de Proteínas"])
sequence_input = st.sidebar.text_area("Introduce la secuencia:", placeholder="Ejemplo: ATGCGA para ADN o MGLSD para proteína")

# Función para análisis de ADN/ARN
def analyze_dna(sequence):
    sequence = Seq(sequence.strip().upper())
    if not all(base in "ACGTU" for base in sequence):
        st.error("La secuencia contiene caracteres no válidos. Solo se permiten A, C, G, T, U.")
        return
    
    length = len(sequence)
    gc_content = 100 * (sequence.count("G") + sequence.count("C")) / length
    complement = sequence.complement()
    transcribed = sequence.transcribe() if "T" in sequence else "No aplica (ARN)"
    
    st.subheader("Resultados del Análisis de ADN/ARN")
    st.markdown(f"""
    - **Longitud de la secuencia:** {length} bases  
    - **Contenido GC:** {gc_content:.2f}%  
    - **Complemento:** {complement}  
    - **Transcripción:** {transcribed}  
    """)
    
    # Gráfico del contenido GC
    fig, ax = plt.subplots()
    ax.bar(["G", "C", "Otros"], [sequence.count("G"), sequence.count("C"), length - sequence.count("G") - sequence.count("C")])
    ax.set_title("Distribución de bases")
    ax.set_ylabel("Frecuencia")
    st.pyplot(fig)

# Función para análisis de proteínas
def analyze_protein(sequence):
    sequence = Seq(sequence.strip().upper())
    if not all(residue in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
        st.error("La secuencia contiene caracteres no válidos. Usa el formato de una letra para aminoácidos.")
        return

    length = len(sequence)
    hydrophobic = sum(sequence.count(res) for res in "AILMFWV")
    hydrophilic = sum(sequence.count(res) for res in "RNDQEGKH")
    seq_three_letter = seq3(str(sequence))
    mw = molecular_weight(sequence)

    st.subheader("Resultados del Análisis de Proteínas")
    st.markdown(f"""
    - **Longitud de la secuencia:** {length} aminoácidos  
    - **Residuos hidrofóbicos:** {hydrophobic} ({100 * hydrophobic / length:.2f}%)  
    - **Residuos hidrofílicos:** {hydrophilic} ({100 * hydrophilic / length:.2f}%)  
    - **Peso molecular:** {mw:.2f} Da  
    - **Secuencia en formato de tres letras:** {seq_three_letter}  
    """)

    # Gráfico de composición de residuos
    fig, ax = plt.subplots()
    ax.pie([hydrophobic, hydrophilic, length - hydrophobic - hydrophilic], labels=["Hidrofóbicos", "Hidrofílicos", "Otros"], autopct='%1.1f%%')
    ax.set_title("Composición de residuos de la proteína")
    st.pyplot(fig)

# Ejecución del análisis
if sequence_input:
    if analysis_type == "Análisis de ADN/ARN":
        analyze_dna(sequence_input)
    else:
        analyze_protein(sequence_input)
