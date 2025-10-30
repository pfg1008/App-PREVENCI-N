# gui.py
# (ACTUALIZADO para eliminar el messagebox de "Procesando...")

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import os
import json
import platform
import subprocess

from logic_engine import run_full_analysis, get_recommendations
from pdf_generator import create_pdf_report

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
        # ... (El menú no cambia) ...
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

        # ... (Los frames 1, 2, 3 y 4 no cambian) ...
        file_frame = ttk.LabelFrame(main_frame, text="1. Cargar Archivo", padding="10")
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        self.file_path_var = tk.StringVar(value="No se ha seleccionado ningún archivo.")
        ttk.Label(file_frame, textvariable=self.file_path_var, wraplength=550).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        ttk.Button(file_frame, text="Seleccionar Archivo (.csv)", command=self.load_csv).pack(side=tk.RIGHT, padx=5)

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

        ttk.Button(info_frame, text="Limpiar", command=self._clear_form).grid(row=2, column=2, columnspan=2, padx=5, pady=5, sticky='e')
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
        
        # --- CAMBIO 1 (Añadir self.status_label) ---
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

        # Esta es la nueva etiqueta que reemplaza al messagebox
        self.status_label = ttk.Label(action_frame, text="Procesando...", font=('Helvetica', 10, 'italic'))
        self.status_label.grid(row=2, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.status_label.grid_remove() # Empezará oculta
        # --- FIN CAMBIO 1 ---

    def _load_cyp2d6_map(self):
        try:
            filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CYP2D6_Diplotype_Phenotype_Table modificada.csv")
            df = pd.read_csv(filepath, sep=';')
            df.columns = df.columns.str.strip()
            
            # Aseguramos que los diplotipos estén ordenados para el mapeo
            # (Incluso si tu script de 'logic_engine' no los ordena,
            # 'gui.py' usaba una versión de logic_engine que NO ordenaba,
            # así que es más seguro estandarizar el mapa aquí).
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
    
    def _save_patient_db(self, patient_info):
        patient_id = patient_info.get("N° Historia")
        if not patient_id: return
        self.patient_db[patient_id] = patient_info
        try:
            with open(self.db_filepath, 'w') as f:
                json.dump(self.patient_db, f, indent=4)
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar la información del paciente: {e}")

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

    # --- CAMBIO 2 (Modificar load_csv) ---
    def load_csv(self):
        """
        Carga el CSV crudo y lanza el motor de análisis completo.
        Muestra un 'label' de estado en lugar de un 'messagebox'.
        """
        filepath = filedialog.askopenfilename(filetypes=(("Archivos CSV", "*.csv"), ("Todos los archivos", "*.*")))
        if not filepath: return
        try:
            temp_df = pd.read_csv(filepath, sep=';', dtype={'Sample/Assay': str})
            
            if 'Sample/Assay' not in temp_df.columns:
                messagebox.showerror("Error de Formato", "El archivo CSV no contiene la columna 'Sample/Assay'.")
                return
            
            temp_df.dropna(subset=['Sample/Assay'], inplace=True)
            self.genotype_df_raw = temp_df.set_index('Sample/Assay')

            # --- ESTE ES EL BLOQUE MODIFICADO ---
            # 1. Quitar el messagebox molesto
            # messagebox.showinfo("Procesando...", "Archivo cargado. Procesando todos los pacientes. Por favor, espere...")
            
            # 2. Mostrar nuestra nueva etiqueta de estado
            self.status_label.config(text="Archivo cargado. Procesando todos los pacientes...")
            self.status_label.grid()
            self.root.update_idletasks() # Forzar que la GUI se actualice AHORA
            # --- FIN BLOQUE MODIFICADO ---
            
            results_df, error = run_full_analysis(self.genotype_df_raw, self.cyp2d6_phenotype_map)
            
            # 3. Ocultar la etiqueta de estado cuando terminemos
            self.status_label.grid_remove()
            
            if error:
                messagebox.showerror("Error de Análisis", error)
                return
            
            self.results_df = results_df
            
            patients = self.results_df.index.tolist()
            self.patient_combobox['values'] = patients
            self.patient_combobox.config(state="readonly")
            self.batch_button.config(state="normal")
            self.file_path_var.set(os.path.basename(filepath))
            
            # El messagebox de "Éxito" sí lo dejamos, porque es útil
            messagebox.showinfo("Éxito", f"Proceso completado. Se encontraron y analizaron {len(patients)} pacientes.")
            
        except Exception as e:
            # 4. Ocultar la etiqueta si hay un error
            self.status_label.grid_remove() 
            messagebox.showerror("Error de Carga", f"No se pudo leer o procesar el archivo.\nError: {e}")

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
            
        except KeyError:
            messagebox.showerror("Error", f"No se encontraron resultados para el paciente {selected_patient_id}.")
            self.generate_button.config(state="disabled")
        except Exception as e:
            messagebox.showerror("Error de Visualización", f"Error al mostrar al paciente {selected_patient_id}:\n{e}")
            self.generate_button.config(state="disabled")

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

    def generate_report(self):
        patient_info = {key: entry.get("1.0", "end-1c").strip() if isinstance(entry, tk.Text) else entry.get().strip() for key, entry in self.entries.items()}
        
        if not all([patient_info['Nombre'], patient_info['Apellidos'], patient_info['N° Historia']]):
            messagebox.showwarning("Faltan Datos", "Por favor, complete al menos Nombre, Apellidos y N° de Historia.")
            return
        if not self.current_genotypes:
            messagebox.showerror("Error", "No se han calculado los resultados. Seleccione un paciente.")
            return
        try:
            self._save_patient_db(patient_info)
            recommendations = get_recommendations(self.current_phenotypes)
            pdf_filename, error = create_pdf_report(patient_info, self.current_genotypes, self.current_phenotypes, recommendations)
            if error:
                messagebox.showerror("Error al crear PDF", error)
            elif messagebox.askyesno("Éxito", f"Se ha guardado el informe:\n{pdf_filename}\n\n¿Desea abrir la carpeta contenedora?"):
                self._open_folder(os.path.dirname(os.path.abspath(pdf_filename)))
        except Exception as e:
            messagebox.showerror("Error de Generación", f"Error inesperado al generar el informe:\n{e}")

    def generate_batch_reports(self):
        """
        Itera sobre el DataFrame de resultados pre-calculados.
        """
        if self.results_df is None:
            messagebox.showwarning("Datos no procesados", "Por favor, cargue primero un archivo de genotipado.")
            return
        if not messagebox.askyesno("Confirmar", f"Se generarán informes para los {len(self.results_df)} pacientes del archivo. El proceso puede tardar.\n\n¿Desea continuar?"):
            return

        success_count, fail_count = 0, 0
        output_folder = "Informes_Lote"
        
        self.progress_bar.grid()
        self.progress_bar.config(maximum=len(self.results_df), value=0)
        
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
                
                if error: fail_count += 1
                else: success_count += 1
            except Exception:
                fail_count += 1
            
            self.progress_bar.step(1)
            self.root.update_idletasks()

        self.progress_bar.grid_remove()
        
        summary_message = f"Proceso completado.\n\nInformes generados: {success_count}\nInformes fallidos: {fail_count}"
        if messagebox.askyesno("Proceso en Lote Terminado", f"{summary_message}\n\nLos archivos se han guardado en la carpeta '{output_folder}'.\n¿Desea abrir esta carpeta?"):
            self._open_folder(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_folder))