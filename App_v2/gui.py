# gui.py
# (ACTUALIZADO con Threading para evitar congelamiento y botón de guardar)

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import os
import json
import platform
import subprocess
import threading  # --- CAMBIO: Importar threading ---
import logging   # --- CAMBIO: Importar logging para errores en lote ---

from logic_engine import run_full_analysis, get_recommendations
from pdf_generator import create_pdf_report

# --- CAMBIO: Configurar un logging básico para errores ---
logging.basicConfig(filename='app_errors.log', 
                    level=logging.ERROR, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Generador de Informes Farmacogenéticos v3.1")
        self.root.geometry("850x750")
        
        self.style = ttk.Style(self.root)
        self.theme_var = tk.StringVar(value='clam')
        try:
            self.style.theme_use(self.theme_var.get())
        except tk.TclError:
            self.theme_var.set(self.style.theme_use())

        self.style.configure('TLabel', font=('Helvetica', 10))
        self.style.configure('Bold.TLabel', font=('Helvetica', 10, 'bold'))
        self.style.configure('TButton', font=('Helvetica', 10, 'bold'))
        self.style.configure('TLabelframe.Label', font=('Helvetica', 11, 'bold'))

        self.genotype_df_raw = None
        self.cyp2d6_phenotype_map = {}
        self.results_df = None
        
        self.current_genotypes, self.current_phenotypes = None, None
        
        self.patient_db = {}
        self.db_filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "patient_data.json")

        self._setup_ui()
        self._load_cyp2d6_map()
        self._load_patient_db()

    def _setup_ui(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        edit_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Editar", menu=edit_menu)
        edit_menu.add_command(label="Limpiar Formulario", command=self._clear_form)
        edit_menu.add_command(label="Abrir Carpeta de Informes", command=self._open_reports_folder)
        edit_menu.add_separator()
        theme_menu = tk.Menu(edit_menu, tearoff=0)
        edit_menu.add_cascade(label="Tema", menu=theme_menu)
        for theme in ['clam', 'alt', 'default', 'vista']:
            theme_menu.add_radiobutton(label=theme, variable=self.theme_var, command=self._change_theme)
            
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Ayuda", menu=help_menu)
        help_menu.add_command(label="Acerca de...", command=self._show_about)
        
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        file_frame = ttk.LabelFrame(main_frame, text="1. Cargar Archivo", padding="10")
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        self.file_path_var = tk.StringVar(value="No se ha seleccionado ningún archivo.")
        ttk.Label(file_frame, textvariable=self.file_path_var, wraplength=550).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        # --- CAMBIO: Guardar referencia al botón para desactivarlo ---
        self.load_button = ttk.Button(file_frame, text="Seleccionar Archivo (.csv)", command=self.load_csv)
        self.load_button.pack(side=tk.RIGHT, padx=5)

        patient_select_frame = ttk.LabelFrame(main_frame, text="2. Seleccionar Paciente", padding="10")
        patient_select_frame.pack(fill=tk.X, padx=5, pady=10)
        self.patient_var = tk.StringVar()
        self.patient_combobox = ttk.Combobox(patient_select_frame, textvariable=self.patient_var, state="disabled", font=('Helvetica', 10))
        self.patient_combobox.pack(fill=tk.X, expand=True)
        self.patient_combobox.bind("<<ComboboxSelected>>", self.on_patient_select)

        info_frame = ttk.LabelFrame(main_frame, text="3. Datos del Paciente", padding="10")
        info_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        info_frame.columnconfigure(1, weight=1)
        info_frame.columnconfigure(3, weight=1)
        
        self.entries = {}
        fields = ["Nombre", "Apellidos", "N° Historia", "Sexo", "Fecha de Nacimiento"]
        for i, field in enumerate(fields):
            row, col = divmod(i, 2)
            ttk.Label(info_frame, text=f"{field}:").grid(row=row, column=col*2, padx=5, pady=5, sticky="w")
            self.entries[field] = ttk.Entry(info_frame, font=('Helvetica', 10))
            self.entries[field].grid(row=row, column=col*2 + 1, padx=5, pady=5, sticky="ew")

        # --- CAMBIO: Añadir botón de Guardar Datos ---
        action_buttons_frame = ttk.Frame(info_frame)
        action_buttons_frame.grid(row=2, column=2, columnspan=2, padx=5, pady=5, sticky='e')
        
        self.save_button = ttk.Button(action_buttons_frame, text="Guardar Datos", command=self._save_current_patient, state="disabled")
        self.save_button.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(action_buttons_frame, text="Limpiar", command=self._clear_form).pack(side=tk.LEFT, padx=5)
        # --- FIN DEL CAMBIO ---

        ttk.Label(info_frame, text="Enfermedades Actuales:").grid(row=3, column=0, padx=5, pady=5, sticky="nw")
        self.entries['Enfermedades'] = tk.Text(info_frame, height=4, width=30, font=('Helvetica', 10))
        self.entries['Enfermedades'].grid(row=3, column=1, padx=5, pady=5, sticky="nsew")
        info_frame.rowconfigure(3, weight=1)
        ttk.Label(info_frame, text="Tratamiento Habitual:").grid(row=3, column=2, padx=5, pady=5, sticky="nw")
        self.entries['Tratamiento'] = tk.Text(info_frame, height=4, width=30, font=('Helvetica', 10))
        self.entries['Tratamiento'].grid(row=3, column=3, padx=5, pady=5, sticky="nsew")

        preview_frame = ttk.LabelFrame(main_frame, text="4. Previsualización", padding="15")
        preview_frame.pack(fill=tk.X, padx=5, pady=10)
        preview_frame.columnconfigure(1, weight=1)
        preview_frame.columnconfigure(3, weight=1)
        self.preview_vars = {
            'DPYD_geno': tk.StringVar(value="..."), 'DPYD_pheno': tk.StringVar(value="..."),
            'CYP2D6_geno': tk.StringVar(value="..."), 'CYP2D6_pheno': tk.StringVar(value="..."),
            'UGT1A1_geno': tk.StringVar(value="..."), 'UGT1A1_pheno': tk.StringVar(value="...")
        }
        self.pheno_labels = {}
        ttk.Label(preview_frame, text="Genotipo DPYD:", style='Bold.TLabel').grid(row=0, column=0, sticky='w', padx=5)
        ttk.Label(preview_frame, textvariable=self.preview_vars['DPYD_geno']).grid(row=0, column=1, sticky='w', padx=5)
        ttk.Label(preview_frame, text="Fenotipo DPYD:", style='Bold.TLabel').grid(row=0, column=2, sticky='w', padx=5)
        self.pheno_labels['DPYD'] = tk.Label(preview_frame, textvariable=self.preview_vars['DPYD_pheno'], font=('Helvetica', 10, 'bold'))
        self.pheno_labels['DPYD'].grid(row=0, column=3, sticky='w', padx=5)
        ttk.Label(preview_frame, text="Genotipo CYP2D6:", style='Bold.TLabel').grid(row=1, column=0, sticky='w', padx=5, pady=4)
        ttk.Label(preview_frame, textvariable=self.preview_vars['CYP2D6_geno']).grid(row=1, column=1, sticky='w', padx=5)
        ttk.Label(preview_frame, text="Fenotipo CYP2D6:", style='Bold.TLabel').grid(row=1, column=2, sticky='w', padx=5)
        self.pheno_labels['CYP2D6'] = tk.Label(preview_frame, textvariable=self.preview_vars['CYP2D6_pheno'], font=('Helvetica', 10, 'bold'))
        self.pheno_labels['CYP2D6'].grid(row=1, column=3, sticky='w', padx=5)
        ttk.Label(preview_frame, text="Genotipo UGT1A1:", style='Bold.TLabel').grid(row=2, column=0, sticky='w', padx=5)
        ttk.Label(preview_frame, textvariable=self.preview_vars['UGT1A1_geno']).grid(row=2, column=1, sticky='w', padx=5)
        ttk.Label(preview_frame, text="Fenotipo UGT1A1:", style='Bold.TLabel').grid(row=2, column=2, sticky='w', padx=5)
        self.pheno_labels['UGT1A1'] = tk.Label(preview_frame, textvariable=self.preview_vars['UGT1A1_pheno'], font=('Helvetica', 10, 'bold'))
        self.pheno_labels['UGT1A1'].grid(row=2, column=3, sticky='w', padx=5)
        
        action_frame = ttk.LabelFrame(main_frame, text="5. Acciones", padding="10")
        action_frame.pack(fill=tk.X, padx=5, pady=5)
        action_frame.columnconfigure(0, weight=1)
        action_frame.columnconfigure(1, weight=1)

        self.generate_button = ttk.Button(action_frame, text="Generar Informe para Paciente", command=self.generate_report, state="disabled")
        self.generate_button.grid(row=0, column=0, sticky='ew', padx=5, ipady=8)
        self.batch_button = ttk.Button(action_frame, text="Generar Todos los Informes", command=self.generate_batch_reports, state="disabled")
        self.batch_button.grid(row=0, column=1, sticky='ew', padx=5, ipady=8)
        
        self.progress_bar = ttk.Progressbar(action_frame, orient='horizontal', mode='determinate')
        self.progress_bar.grid(row=1, column=0, columnspan=2, sticky='ew', padx=5, pady=8)
        self.progress_bar.grid_remove()

        self.status_label = ttk.Label(action_frame, text="Procesando...", font=('Helvetica', 10, 'italic'))
        self.status_label.grid(row=2, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.status_label.grid_remove() 

    def _load_cyp2d6_map(self):
        try:
            filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CYP2D6_Diplotype_Phenotype_Table modificada.csv")
            df = pd.read_csv(filepath, sep=';')
            df.columns = df.columns.str.strip()
            
            # --- CAMBIO: CORREGIDO ---
            # Esta lógica AHORA coincide con la nueva lógica de 'logic_engine.py'
            # (ambos ordenan el diplotipo).
            df['CYP2D6 Diplotype'] = df['CYP2D6 Diplotype'].apply(
                lambda x: '/'.join(sorted(str(x).split('/')))
            )
            df = df.drop_duplicates(subset=['CYP2D6 Diplotype'])
            
            self.cyp2d6_phenotype_map = pd.Series(
                df['Coded Diplotype/Phenotype Summary'].values, 
                index=df['CYP2D6 Diplotype']
            ).to_dict()
        except Exception as e:
            messagebox.showerror("Error", f"No se encontró o no se pudo leer el archivo de fenotipos CYP2D6.\nError: {e}")
            self.root.quit()

    def _load_patient_db(self):
        if os.path.exists(self.db_filepath):
            try:
                with open(self.db_filepath, 'r') as f:
                    self.patient_db = json.load(f)
            except Exception as e:
                messagebox.showwarning("Advertencia", f"No se pudo cargar la base de datos de pacientes: {e}")
                self.patient_db = {}
    
    # --- CAMBIO: Función de ayuda para guardar ---
    def _save_patient_db(self, patient_info):
        patient_id = patient_info.get("N° Historia")
        if not patient_id: 
            return False
        self.patient_db[patient_id] = patient_info
        try:
            with open(self.db_filepath, 'w') as f:
                json.dump(self.patient_db, f, indent=4)
            return True
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar la información del paciente: {e}")
            return False

    def _open_folder(self, path):
        try:
            if platform.system() == "Windows":
                os.startfile(path)
            elif platform.system() == "Darwin":
                subprocess.run(["open", path], check=True)
            else:
                subprocess.run(["xdg-open", path], check=True)
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo abrir la carpeta:\n{e}")
            
    def _open_reports_folder(self):
        folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Informes_Lote")
        if not os.path.isdir(folder_path):
            folder_path = os.path.dirname(os.path.abspath(__file__))
        self._open_folder(folder_path)

    def _clear_form(self):
        self.patient_combobox.set('')
        for entry in self.entries.values():
            if isinstance(entry, tk.Text):
                entry.delete("1.0", tk.END)
            else:
                entry.delete(0, tk.END)
        for key in self.preview_vars:
            self.preview_vars[key].set("...")
        self._update_phenotype_colors()
        self.generate_button.config(state="disabled")
        self.save_button.config(state="disabled") # --- CAMBIO ---
        
    def _show_about(self):
        about_text = (
            "Generador de Informes Farmacogenéticos v3.1 (Integrado con motor Pandas)\n\n"
            "Esta aplicación ha sido desarrollada para la asignatura de 'Prevención y terapéutica de precisión' del Grado en Ingeniería de la Salud.\n\n"
            "Desarrolladores:\n"
            " • Pablo Fernández\n"
            " • Guillermo de la Fuente\n"
            " • Margarita Aguzarova\n"
            " • Daniel Marina\n"
            " • Darío Meneses\n"
            " • Isabel García\n"
            " • María García"
        )
        messagebox.showinfo("Acerca de", about_text)

    def _change_theme(self):
        try:
            self.style.theme_use(self.theme_var.get())
        except tk.TclError:
            messagebox.showwarning("Tema no disponible", f"El tema '{self.theme_var.get()}' no está disponible en este sistema operativo.")
            self.theme_var.set(self.style.theme_use())

    # --- CAMBIO: Función modificada para usar Threading ---
    def load_csv(self):
        """
        Abre el diálogo para cargar CSV y lanza el análisis EN UN HILO SEPARADO
        para no congelar la GUI.
        """
        filepath = filedialog.askopenfilename(filetypes=(("Archivos CSV", "*.csv"), ("Todos los archivos", "*.*")))
        if not filepath: return
        try:
            # La lectura inicial del CSV es rápida, se puede hacer aquí
            temp_df = pd.read_csv(filepath, sep=';', dtype={'Sample/Assay': str})
            
            if 'Sample/Assay' not in temp_df.columns:
                messagebox.showerror("Error de Formato", "El archivo CSV no contiene la columna 'Sample/Assay'.")
                return
            
            temp_df.dropna(subset=['Sample/Assay'], inplace=True)
            self.genotype_df_raw = temp_df.set_index('Sample/Assay')

            # --- Preparar la GUI para la carga ---
            self.status_label.config(text="Archivo cargado. Procesando todos los pacientes...")
            self.status_label.grid()
            self._set_ui_state("disabled") # Desactiva botones
            self.root.update_idletasks()
            
            # --- Lanzar la tarea pesada (run_full_analysis) en un hilo ---
            def task():
                # Esta función se ejecuta en el hilo secundario
                results_df, error = run_full_analysis(self.genotype_df_raw, self.cyp2d6_phenotype_map)
                # Cuando termina, llama a 'on_analysis_complete' en el hilo principal
                self.root.after(0, self.on_analysis_complete, results_df, error, os.path.basename(filepath))
            
            threading.Thread(target=task, daemon=True).start()

        except Exception as e:
            self.status_label.grid_remove() 
            self._set_ui_state("normal") # Reactiva botones si falla
            messagebox.showerror("Error de Carga", f"No se pudo leer o procesar el archivo.\nError: {e}")

    # --- CAMBIO: Nueva función callback para el hilo de 'load_csv' ---
    def on_analysis_complete(self, results_df, error, basename):
        """
        Se ejecuta en el hilo principal cuando 'run_full_analysis' termina.
        """
        self.status_label.grid_remove()
        self._set_ui_state("normal") # Reactiva botones

        if error:
            messagebox.showerror("Error de Análisis", error)
            return
        
        self.results_df = results_df
        
        patients = self.results_df.index.tolist()
        self.patient_combobox['values'] = patients
        self.patient_combobox.config(state="readonly")
        self.batch_button.config(state="normal")
        self.file_path_var.set(basename)
        
        messagebox.showinfo("Éxito", f"Proceso completado. Se encontraron y analizaron {len(patients)} pacientes.")

    def on_patient_select(self, event=None):
        """
        Busca los resultados pre-calculados.
        """
        selected_patient_id = self.patient_var.get()
        if not selected_patient_id or self.results_df is None: 
            return
            
        self._clear_form()
        
        if selected_patient_id in self.patient_db:
            for key, value in self.patient_db[selected_patient_id].items():
                if key in self.entries and isinstance(self.entries[key], tk.Text):
                    self.entries[key].insert("1.0", value)
                elif key in self.entries:
                    self.entries[key].insert(0, value)
        
        self.entries["N° Historia"].insert(0, selected_patient_id)
        
        try:
            patient_results = self.results_df.loc[selected_patient_id]
            
            self.current_genotypes = {
                'DPYD': patient_results['DPYD'],
                'CYP2D6': patient_results['CYP2D6'],
                'UGT1A1': patient_results['UGT1A1']
            }
            self.current_phenotypes = {
                'DPYD': patient_results['Fenotipo_DPYD'],
                'CYP2D6': patient_results['Fenotipo_CYP2D6'],
                'UGT1A1': patient_results['Fenotipo_UGT1A1']
            }
            
            for gene in ['DPYD', 'CYP2D6', 'UGT1A1']:
                self.preview_vars[f'{gene}_geno'].set(self.current_genotypes[gene])
                self.preview_vars[f'{gene}_pheno'].set(self.current_phenotypes[gene])
                
            self._update_phenotype_colors()
            self.generate_button.config(state="normal")
            self.save_button.config(state="normal") # --- CAMBIO ---
            
        except KeyError:
            messagebox.showerror("Error", f"No se encontraron resultados para el paciente {selected_patient_id}.")
            self.generate_button.config(state="disabled")
            self.save_button.config(state="disabled") # --- CAMBIO ---
        except Exception as e:
            messagebox.showerror("Error de Visualización", f"Error al mostrar al paciente {selected_patient_id}:\n{e}")
            self.generate_button.config(state="disabled")
            self.save_button.config(state="disabled") # --- CAMBIO ---

    def _update_phenotype_colors(self):
        color_map = {
            'Metabolizador normal': '#008000', 'Metabolizador ultrarrápido': '#008000',
            'Metabolizador intermedio': '#FFA500',
            'Metabolizador lento': '#FF0000',
            '...': 'black', 'Indeterminado': 'black'
        }
        for gene in ['DPYD', 'CYP2D6', 'UGT1A1']:
            pheno_text = self.preview_vars[f'{gene}_pheno'].get()
            color = color_map.get(pheno_text, 'black')
            self.pheno_labels[gene].config(fg=color)

    # --- CAMBIO: Nueva función de ayuda para obtener datos del form ---
    def _get_patient_info_from_form(self):
        """Recoge los datos de los campos de entrada en un diccionario."""
        return {key: entry.get("1.0", "end-1c").strip() if isinstance(entry, tk.Text) else entry.get().strip() for key, entry in self.entries.items()}

    # --- CAMBIO: Nueva función para el botón "Guardar Datos" ---
    def _save_current_patient(self):
        """Guarda la información del paciente actual en el JSON."""
        patient_info = self._get_patient_info_from_form()
        if not patient_info.get("N° Historia"):
            messagebox.showwarning("Faltan Datos", "El 'N° Historia' es necesario para guardar.")
            return
        
        if self._save_patient_db(patient_info):
            messagebox.showinfo("Guardado", f"Datos del paciente {patient_info['N° Historia']} guardados.")

    def generate_report(self):
        patient_info = self._get_patient_info_from_form()
        
        if not all([patient_info['Nombre'], patient_info['Apellidos'], patient_info['N° Historia']]):
            messagebox.showwarning("Faltan Datos", "Por favor, complete al menos Nombre, Apellidos y N° de Historia.")
            return
        if not self.current_genotypes:
            messagebox.showerror("Error", "No se han calculado los resultados. Seleccione un paciente.")
            return
        try:
            # --- CAMBIO: Guardar siempre los datos antes de generar ---
            self._save_patient_db(patient_info)
            
            recommendations = get_recommendations(self.current_phenotypes)
            pdf_filename, error = create_pdf_report(patient_info, self.current_genotypes, self.current_phenotypes, recommendations)
            
            if error:
                messagebox.showerror("Error al crear PDF", error)
            elif messagebox.askyesno("Éxito", f"Se ha guardado el informe:\n{pdf_filename}\n\n¿Desea abrir la carpeta contenedora?"):
                self._open_folder(os.path.dirname(os.path.abspath(pdf_filename)))
        except Exception as e:
            messagebox.showerror("Error de Generación", f"Error inesperado al generar el informe:\n{e}")

    # --- CAMBIO: Función modificada para usar Threading ---
    def generate_batch_reports(self):
        """
        Lanza la generación de informes en lote EN UN HILO SEPARADO
        para no congelar la GUI.
        """
        if self.results_df is None:
            messagebox.showwarning("Datos no procesados", "Por favor, cargue primero un archivo de genotipado.")
            return
        if not messagebox.askyesno("Confirmar", f"Se generarán informes para los {len(self.results_df)} pacientes del archivo. El proceso puede tardar.\n\n¿Desea continuar?"):
            return

        # --- Preparar la GUI para la carga ---
        self.progress_bar.grid()
        self.progress_bar.config(maximum=len(self.results_df), value=0)
        self._set_ui_state("disabled") # Desactiva botones
        
        # --- Lanzar la tarea pesada en un hilo ---
        def task():
            success_count, fail_count = 0, 0
            output_folder = "Informes_Lote"
        
            for patient_id, row in self.results_df.iterrows():
                try:
                    genotypes = {
                        'DPYD': row['DPYD'],
                        'CYP2D6': row['CYP2D6'],
                        'UGT1A1': row['UGT1A1']
                    }
                    phenotypes = {
                        'DPYD': row['Fenotipo_DPYD'],
                        'CYP2D6': row['Fenotipo_CYP2D6'],
                        'UGT1A1': row['Fenotipo_UGT1A1']
                    }
                    
                    recommendations = get_recommendations(phenotypes)
                    patient_info = self.patient_db.get(patient_id, {"N° Historia": patient_id})
                    
                    _, error = create_pdf_report(patient_info, genotypes, phenotypes, recommendations, folder=output_folder)
                    
                    if error: 
                        fail_count += 1
                        logging.error(f"Fallo al generar PDF para {patient_id}: {error}")
                    else: 
                        success_count += 1
                except Exception as e:
                    fail_count += 1
                    logging.error(f"Fallo crítico al procesar paciente {patient_id}: {e}")
                
                # Actualiza la barra de progreso desde el hilo principal
                self.root.after(0, self.progress_bar.step, 1)
            
            # Llama al callback final en el hilo principal
            self.root.after(0, self.on_batch_complete, success_count, fail_count, output_folder)
        
        threading.Thread(target=task, daemon=True).start()

    # --- CAMBIO: Nueva función callback para el hilo de 'generate_batch_reports' ---
    def on_batch_complete(self, success_count, fail_count, output_folder):
        """
        Se ejecuta en el hilo principal cuando 'generate_batch_reports' termina.
        """
        self.progress_bar.grid_remove()
        self._set_ui_state("normal") # Reactiva botones
        
        summary_message = f"Proceso completado.\n\nInformes generados: {success_count}\nInformes fallidos: {fail_count}"
        
        if fail_count > 0:
            summary_message += f"\n\nSe registraron {fail_count} errores en el archivo 'app_errors.log'."

        if messagebox.askyesno("Proceso en Lote Terminado", f"{summary_message}\n\nLos archivos se han guardado en la carpeta '{output_folder}'.\n¿Desea abrir esta carpeta?"):
            self._open_folder(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_folder))

    # --- CAMBIO: Nueva función de ayuda para activar/desactivar la GUI ---
    def _set_ui_state(self, state="normal"):
        """Activa o desactiva los botones principales."""
        self.load_button.config(state=state)
        self.batch_button.config(state=state)
        self.generate_button.config(state=state if state == "normal" and self.current_genotypes else "disabled")
        self.save_button.config(state=state if state == "normal" and self.current_genotypes else "disabled")
        self.patient_combobox.config(state=state if state == "normal" else "disabled")