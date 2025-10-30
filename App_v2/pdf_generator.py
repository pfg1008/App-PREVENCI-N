# pdf_generator.py

import os
from datetime import datetime
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from reportlab.platypus import Table, TableStyle, Paragraph
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

def create_pdf_report(patient_info, genotypes, phenotypes, recommendations, folder=""):
    """Generates the final PDF report with hyperlinks and bold keywords."""
    try:
        if folder:
            os.makedirs(folder, exist_ok=True)
            filename = os.path.join(folder, f"Informe_PGx_{patient_info['N° Historia']}_{datetime.now().strftime('%Y%m%d')}.pdf")
        else:
            filename = f"Informe_PGx_{patient_info['N° Historia']}_{datetime.now().strftime('%Y%m%d')}.pdf"

        c = canvas.Canvas(filename, pagesize=A4)
        width, height = A4

        c.setFont("Helvetica-Bold", 16)
        c.drawCentredString(width / 2.0, height - 3*cm, "INFORME FARMACOGENÉTICO - PERFIL ONCOLOGÍA")
        
        text_y = height - 4.5*cm
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, text_y, "PACIENTE:")
        c.drawString(11*cm, text_y, "N° HISTORIA:")
        c.setFont("Helvetica", 11)
        c.drawString(4.5*cm, text_y, f"{patient_info.get('Nombre', '')} {patient_info.get('Apellidos', '')}")
        c.drawString(14*cm, text_y, patient_info.get('N° Historia', ''))

        text_y -= 0.7*cm
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, text_y, "FECHA DE NACIMIENTO:")
        c.drawString(11*cm, text_y, "SEXO:")
        c.setFont("Helvetica", 11)
        c.drawString(6.5*cm, text_y, patient_info.get('Fecha de Nacimiento', ''))
        c.drawString(12.5*cm, text_y, patient_info.get('Sexo', ''))
        
        text_y -= 1*cm
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, text_y, "ENFERMEDADES ACTUALES:")
        text_y -= 0.6*cm
        text_object = c.beginText(2*cm, text_y)
        text_object.setFont("Helvetica", 10)
        for line in patient_info.get('Enfermedades', '').split('\n'):
            text_object.textLine(line)
        c.drawText(text_object)
        
        text_y -= (len(patient_info.get('Enfermedades', '').split('\n')) * 0.4 + 0.8) * cm
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, text_y, "TRATAMIENTO HABITUAL:")
        text_y -= 0.6*cm
        text_object = c.beginText(2*cm, text_y)
        text_object.setFont("Helvetica", 10)
        for line in patient_info.get('Tratamiento', '').split('\n'):
            text_object.textLine(line)
        c.drawText(text_object)
        
        table_top_y = text_y - (len(patient_info.get('Tratamiento', '').split('\n')) * 0.4 + 1.2) * cm

        c.line(2*cm, table_top_y, width - 2*cm, table_top_y)
        c.setFont("Helvetica-Bold", 14)
        c.drawString(2*cm, table_top_y - 0.7*cm, "RESULTADOS")

        # Define URLs for hyperlinks
        GUIDELINE_URLS = {
            "DPYD": "https://www.clinpgx.org/chemical/PA128406956/guidelineAnnotation/PA166122686",
            "CYP2D6": "https://www.clinpgx.org/chemical/PA451581/guidelineAnnotation/PA166176068",
            "UGT1A1": "https://www.clinpgx.org/chemical/PA450085/guidelineAnnotation/PA166104951"
        }

        styles = getSampleStyleSheet()
        cell_style = ParagraphStyle('cell_style', parent=styles['Normal'], fontSize=9, leading=12)
        link_style = ParagraphStyle('link_style', parent=cell_style, textColor=colors.blue, fontName='Helvetica-Bold')
        header_style = ParagraphStyle('header_style', parent=styles['Normal'], fontSize=10, textColor=colors.whitesmoke, fontName='Helvetica-Bold', alignment=1)
        
        data = [
            [Paragraph(col, header_style) for col in ['Gen', 'Genotipo', 'Fenotipo', 'Fármaco', 'Recomendación']],
            [
                Paragraph(f'<link href="{GUIDELINE_URLS["DPYD"]}" color="blue"><u>DPYD</u></link>', cell_style),
                Paragraph(genotypes.get('DPYD', 'N/A'), cell_style),
                Paragraph(phenotypes.get('DPYD', 'N/A'), cell_style),
                Paragraph('Fluorouracilo,<br/>Capecitabina,<br/>Tegafur', cell_style),
                Paragraph(recommendations.get('DPYD', ''), cell_style)
            ],
            [
                Paragraph(f'<link href="{GUIDELINE_URLS["CYP2D6"]}" color="blue"><u>CYP2D6</u></link>', cell_style),
                Paragraph(genotypes.get('CYP2D6', 'N/A'), cell_style),
                Paragraph(phenotypes.get('CYP2D6', 'N/A'), cell_style),
                Paragraph('Tamoxifeno', cell_style),
                Paragraph(recommendations.get('CYP2D6', ''), cell_style)
            ],
            [
                Paragraph(f'<link href="{GUIDELINE_URLS["UGT1A1"]}" color="blue"><u>UGT1A1</u></link>', cell_style),
                Paragraph(genotypes.get('UGT1A1', 'N/A'), cell_style),
                Paragraph(phenotypes.get('UGT1A1', 'N/A'), cell_style),
                Paragraph('Irinotecan', cell_style),
                Paragraph(recommendations.get('UGT1A1', ''), cell_style)
            ]
        ]

        table = Table(data, colWidths=[1.8*cm, 2.7*cm, 3.5*cm, 3.5*cm, 5.5*cm])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.grey),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('GRID', (0,0), (-1,-1), 1, colors.black),
            ('TOPPADDING', (0,0), (-1,-1), 6),
            ('BOTTOMPADDING', (0,0), (-1,-1), 6),
        ]))
        
        table_width, table_height = table.wrapOn(c, width, height)
        table_y = table_top_y - 1*cm - table_height
        table.drawOn(c, 2*cm, table_y)

        c.setFont("Helvetica-Oblique", 8)
        c.drawString(2*cm, 3*cm, "Este informe ha sido elaborado de acuerdo a las guías clínicas del Consorcio para la Implementación de la Farmacogenética Clínica (CPIC).")

        c.save()
        return filename, None
    except Exception as e:
        return None, str(e)