# logic_engine.py
# VERSIÓN CON 'combinar_gen' SIN ORDENACIÓN FINAL 'sorted()'

import pandas as pd
import json
import os

# === 1. Funciones del motor de análisis (de tu script) ===

def mapear_nts_a_alelos(serie_snp: pd.Series, reglas_map: dict):
    """
    Toma una columna (Serie) de genotipos (ej: 'C/G') y la convierte en 
    una columna de alelos mapeados (ej: ['*1', '*2A']) usando un 
    diccionario de reglas (cargado desde el JSON).
    """
    nombre_columna_snp = serie_snp.name
    
    def convertir_celda_a_alelos(celda_genotipo: str, nombre_col: str):
        """
        Función interna. Convierte una sola celda (ej: 'T/T') 
        a una lista de alelos (ej: ['*3', '*3']).
        """
        # Comprobación de seguridad por si llegan datos no-string
        if not isinstance(celda_genotipo, str):
            celda_genotipo = str(celda_genotipo)
            
        nucleotidos = celda_genotipo.split('/')
        alelos_resultantes = []
        
        for nt in nucleotidos:
            # LÓGICA ORIGINAL: Si no se encuentra el nucleótido, se asume '*1'
            alelo_mapeado = reglas_map.get(nombre_col, dict()).get(nt, '*1')
            alelos_resultantes.append(alelo_mapeado)
            
        return alelos_resultantes
    
    # Aplica la conversión celda por celda a toda la columna
    serie_mapeada = serie_snp.apply(lambda celda: convertir_celda_a_alelos(celda, nombre_columna_snp))
    return serie_mapeada


def combinar_gen(lista_de_alelos_por_snp: list, gen: str):
    """
    Combina los alelos de múltiples SNPs de un gen en un diplotipo final.
    (Versión actualizada con lógica de prioridad CYP2D6 detallada).
    """
    # Tu nueva lógica de separación de alelos
    alelos_unicos = list(
        a for sublist in lista_de_alelos_por_snp 
        for a in sublist 
        if not a in ['*1','*4','*10','*10*4']
    )

    alelos_especiales = list(
        a for sublist in lista_de_alelos_por_snp 
        for a in sublist 
        if  a in ['*4','*10','*10*4']
    )

    finales = []
    
    if gen == 'CYP2D6':
        # Tu nueva lógica de bucles 'while'
        tiene_10_4 = '*10*4' in alelos_especiales
        tiene_10 = '*10' in alelos_especiales
        tiene_4 = '*4' in alelos_especiales
        
        flag1=True
        while flag1:
            if tiene_10_4 and tiene_10 and tiene_4:
                alelos_especiales.remove('*10')
                alelos_especiales.remove('*10*4')
                alelos_especiales.remove('*4')
                finales.append('*4')
                tiene_10_4 = '*10*4' in alelos_especiales
                tiene_10 = '*10' in alelos_especiales
                tiene_4 = '*4' in alelos_especiales
            else: flag1 = False
        
        flag2 = True
        while flag2:
            if tiene_10_4 and tiene_10:
                alelos_especiales.remove('*10')
                alelos_especiales.remove('*10*4')
                tiene_10_4 = '*10*4' in alelos_especiales
                tiene_10 = '*10' in alelos_especiales      
                finales.append('*10')
            else: flag2 = False
        
        for i in range(alelos_especiales.count('*10')):
            alelos_especiales.remove('*10')
            
    restantes = alelos_unicos + alelos_especiales
    
    # --- LÓGICA DE RETORNO ORIGINAL (SIN 'sorted()') ---
    
    if len(finales) == 0:
        if len(restantes)==0:
            return '*1/*1'
        elif len(restantes) == 1:
            return f'*1/{restantes[0]}'
        else: 
            return f'{restantes[0]}/{restantes[1]}'
            
    elif len(finales) == 1:
        if len(restantes) == 1:
            return f'{finales[0]}/{restantes[0]}'
        elif len(restantes) == 0:
            return f'*1/{finales[0]}'
    
    elif len(finales)==2:
         return f'{finales[0]}/{finales[1]}'

    # Fallback (si la lógica no cubre un caso, ej: finales=1, restantes=2)
    # Devolvemos un valor que probablemente falle el mapeo
    return 'Error/Logica'


def fenotipo_dpyd(geno: str) -> str:
    """Asigna fenotipo de DPYD basado en la 'dosis' de alelos *1."""
    if geno == '*1/*1':
        return 'Metabolizador normal'
    elif geno.count('*1') == 1:
        return 'Metabolizador intermedio'
    else:
        return 'Metabolizador lento'

def fenotipo_ugt1a1(geno: str) -> str:
    """Asigna fenotipo de UGT1A1 basado en la 'dosis' de alelos *1."""
    if geno == '*1/*1':
        return 'Metabolizador normal'
    elif geno.count('*1') == 1:
        return 'Metabolizador intermedio'
    else:
        return 'Metabolizador lento'


# === 2. Función principal (wrapper) que la GUI llamará ===

def run_full_analysis(df_genotipos_raw, cyp2d6_phenotype_map):
    """
    Función principal que ejecuta todo el pipeline de análisis de pandas.
    Toma el DataFrame crudo y el mapa de fenotipos de CYP2D6.
    Devuelve un DataFrame final con todos los resultados.
    """
    
    # === 1. CARGAR REGLAS ===
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, 'reglas_alelos.json')
    try:
        with open(json_path, 'r') as f:
            mapa_reglas_alelos = json.load(f)
    except FileNotFoundError:
        return None, "Error: No se encontró 'reglas_alelos.json'. Asegúrate de que está en la misma carpeta."
    except Exception as e:
        return None, f"Error al leer 'reglas_alelos.json': {e}"

    # === 2. PROCESAR GENOTIPOS RAW ===
    df_genotipos_para_procesar = df_genotipos_raw.copy()
    
    nombres_columnas_limpios = list(map(lambda x: x.replace('*','_'), df_genotipos_para_procesar.columns))
    df_genotipos_para_procesar.columns = nombres_columnas_limpios
    
    df_alelos_mapeados = df_genotipos_para_procesar.apply(
        lambda columna: mapear_nts_a_alelos(columna, mapa_reglas_alelos)
    )
    
    # === 3. COMBINAR ALELOS POR GEN ===
    df_resultados_finales = pd.DataFrame(index=df_genotipos_para_procesar.index)
    
    df_resultados_finales['DPYD'] = df_alelos_mapeados[
        [c for c in df_alelos_mapeados.columns if 'DPYD' in c]
    ].apply(lambda fila: combinar_gen(fila.tolist(), 'DPYD'), axis=1) # .tolist() es importante

    df_resultados_finales['UGT1A1'] = df_alelos_mapeados[
        [c for c in df_alelos_mapeados.columns if 'UGT1A1' in c]
    ].apply(lambda fila: combinar_gen(fila.tolist(), 'UGT1A1'), axis=1)

    df_resultados_finales['CYP2D6'] = df_alelos_mapeados[
        [c for c in df_alelos_mapeados.columns if 'CYP2D6' in c]
    ].apply(lambda fila: combinar_gen(fila.tolist(), 'CYP2D6'), axis=1)
    
    # === 4. ASIGNAR FENOTIPOS ===
    df_resultados_finales['Fenotipo_DPYD'] = df_resultados_finales['DPYD'].apply(fenotipo_dpyd)
    df_resultados_finales['Fenotipo_UGT1A1'] = df_resultados_finales['UGT1A1'].apply(fenotipo_ugt1a1)
    
    # Lógica de mapeo de fenotipo CYP2D6
    def map_cyp2d6_pheno(geno_str, pheno_map):
        # La GUI pasa 'cyp2d6_phenotype_map' (un dict)
        # 'geno_str' NO ESTÁ ORDENADO (ej. *17/*3) según tu lógica.
        phenotype_info = pheno_map.get(geno_str, "Indeterminate")
        
        phenotype_en = phenotype_info.split(';')[0].strip().replace(" Metabolizer", "")
        
        phenotype_map_es = {
            "Normal": "Metabolizador normal",
            "Intermediate": "Metabolizador intermedio",
            "Poor": "Metabolizador lento",
            "Ultrarapid": "Metabolizador ultrarrápido",
            "Indeterminate": "Indeterminado"
        }
        return phenotype_map_es.get(phenotype_en, "Indeterminado")

    df_resultados_finales['Fenotipo_CYP2D6'] = df_resultados_finales['CYP2D6'].apply(
        lambda geno: map_cyp2d6_pheno(geno, cyp2d6_phenotype_map)
    )

    return df_resultados_finales, None # Devuelve el DF y no-error


# === 3. Función de recomendaciones (la mantenemos del anterior) ===

def get_recommendations(phenotypes):
    """
    Devuelve recomendaciones terapéuticas.
    """
    recs = {}
    
    # Palabras clave para poner en negrita
    keywords = [
        "reducir dosis", "reducción de la dosis", "evitar el uso", "riesgo aumentado",
        "alternativa terapéutica", "fármaco alternativo", "terapia endocrina alternativa",
        "fracaso terapéutico", "toxicidad grave o mortal", "toxicidad aumentada"
    ]

    def highlight_keywords(text):
        """Envuelve las palabras clave en <b> para ReportLab."""
        for keyword in keywords:
            text = text.replace(keyword, f"<b>{keyword}</b>")
            text = text.replace(keyword.capitalize(), f"<b>{keyword.capitalize()}</b>")
        return text

    # DPYD
    dpyd_pheno = phenotypes.get('DPYD', 'Indeterminado')
    if dpyd_pheno == 'Metabolizador normal':
        rec_text = 'Dosis estándar según ficha técnica.'
    elif dpyd_pheno == 'Metabolizador intermedio':
        rec_text = 'Riesgo de toxicidad aumentada. Considerar una reducción de la dosis inicial del 50% seguida de titulación según tolerancia.'
    elif dpyd_pheno == 'Metabolizador lento':
        rec_text = 'Alto riesgo de toxicidad grave o mortal. Evitar el uso de fluoropirimidinas. Considerar un fármaco alternativo.'
    else:
        rec_text = 'Fenotipo indeterminado. No se pueden realizar recomendaciones.'
    recs['DPYD'] = highlight_keywords(rec_text)

    # CYP2D6
    cyp2d6_pheno = phenotypes.get('CYP2D6', 'Indeterminado')
    if cyp2d6_pheno in ['Metabolizador normal', 'Metabolizador ultrarrápido']:
        rec_text = 'Dosis estándar según ficha técnica.'
    elif cyp2d6_pheno == 'Metabolizador intermedio':
        rec_text = 'Riesgo de menor eficacia. Considerar una terapia endocrina alternativa (ej. inhibidor de la aromatasa).'
    elif cyp2d6_pheno == 'Metabolizador lento':
        rec_text = 'Alto riesgo de fracaso terapéutico. Se recomienda el uso de una terapia endocrina alternativa (ej. inhibidor de la aromatasa).'
    else:
        rec_text = 'Fenotipo indeterminado. No se pueden realizar recomendaciones.'
    recs['CYP2D6'] = highlight_keywords(rec_text)

    # UGT1A1
    ugt1a1_pheno = phenotypes.get('UGT1A1', 'Indeterminado')
    if ugt1a1_pheno == 'Metabolizador normal':
        rec_text = 'Dosis estándar según ficha técnica.'
    elif ugt1a1_pheno == 'Metabolizador intermedio':
        rec_text = 'Riesgo aumentado de neutropenia. Considerar iniciar con la dosis estándar; vigilar toxicidad hematológica.'
    elif ugt1a1_pheno == 'Metabolizador lento':
        rec_text = 'Alto riesgo de neutropenia grave. Se recomienda una reducción de la dosis inicial de al menos un 30%.'
    else:
        rec_text = 'Fenotipo indeterminado. No se pueden realizar recomendaciones.'
    recs['UGT1A1'] = highlight_keywords(rec_text)
    
    return recs