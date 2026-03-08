#!/usr/bin/env python3
"""
=============================================================================
DESeq2 & rMATS Pipeline — Interactive GUI Launcher
=============================================================================
Run this script to configure and execute the analysis pipeline without
editing any code.

Usage:
    python pipeline_launcher.py

Requires the pipeline script in the same directory:
    deseq2_rmats_filter_pipeline.py
=============================================================================
"""

import sys
import threading
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser, scrolledtext

# ---------------------------------------------------------------------------
# Locate the pipeline module (same directory as this launcher)
# ---------------------------------------------------------------------------
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

try:
    import deseq2_rmats_filter_pipeline as pipeline
except ImportError as e:
    import tkinter as _tk
    _r = _tk.Tk(); _r.withdraw()
    messagebox.showerror(
        "Import Error",
        f"Could not import deseq2_rmats_filter_pipeline.py:\n{e}\n\n"
        "Make sure it is in the same folder as this launcher."
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Thread-safe stdout/stderr redirect to a Tkinter Text widget
# ---------------------------------------------------------------------------

class _TextRedirector:
    def __init__(self, widget, app):
        self.widget = widget
        self.app = app

    def write(self, s):
        self.app.after(0, self._append, s)

    def _append(self, s):
        self.widget.configure(state="normal")
        self.widget.insert("end", s)
        self.widget.see("end")
        self.widget.configure(state="disabled")

    def flush(self):
        pass


# ---------------------------------------------------------------------------
# Default column names (mirrors pipeline defaults)
# ---------------------------------------------------------------------------

_DESEQ2_COL_DEFAULTS = {
    "gene_id":   "gene_id",
    "gene_name": "gene_name",
    "log2fc":    "log2FoldChange",
    "basemean":  "baseMean",
    "padj":      "padj",
    "pvalue":    "pvalue",
    "biotype":   "biotype",
}

_RMATS_COL_DEFAULTS = {
    "event_id":      "ID",
    "gene_id":       "GeneID",
    "gene_name":     "geneSymbol",
    "pvalue":        "PValue",
    "fdr":           "FDR",
    "inclevel_diff": "IncLevelDifference",
}

_DESEQ2_COL_LABELS = {
    "gene_id":   "Gene ID column",
    "gene_name": "Gene name column",
    "log2fc":    "log2FoldChange column",
    "basemean":  "baseMean column",
    "padj":      "Adjusted p-value column",
    "pvalue":    "Raw p-value column",
    "biotype":   "Biotype column",
}

_RMATS_COL_LABELS = {
    "event_id":      "Event ID column",
    "gene_id":       "Gene ID column",
    "gene_name":     "Gene name column",
    "pvalue":        "P-value column",
    "fdr":           "FDR column",
    "inclevel_diff": "IncLevelDifference column",
}

MAX_CONDITIONS = 5


# ---------------------------------------------------------------------------
# Main application window
# ---------------------------------------------------------------------------

class PipelineLauncherApp(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("DESeq2 & rMATS Analysis Pipeline Launcher")
        self.minsize(960, 780)
        self.resizable(True, True)

        # Internal state
        self.condition_rows: list[dict] = []   # one dict per condition row

        # Threshold variables
        self.log2fc_var    = tk.StringVar(value="0.4")
        self.basemean_var  = tk.StringVar(value="20")
        self.padj_var      = tk.StringVar(value="0.01")
        self.auto_bio_var  = tk.BooleanVar(value=True)
        self.rmats_fdr_var = tk.StringVar(value="0.01")
        self.rmats_pval_var= tk.StringVar(value="0.01")
        self.dpsi_var      = tk.StringVar(value="0.1")
        self.use_fdr_var   = tk.BooleanVar(value=True)

        # Output / figure variables
        default_output = ""
        self.output_dir_var   = tk.StringVar(value=default_output)
        self.fig_format_var   = tk.StringVar(value="png")
        self.dpi_var          = tk.StringVar(value="300")
        self.font_size_var    = tk.StringVar(value="12")
        self.color_up_var     = tk.StringVar(value="#E69F00")
        self.color_down_var   = tk.StringVar(value="#1B98E0")
        self.color_ns_var     = tk.StringVar(value="#BFBFBF")
        self.interactive_var  = tk.BooleanVar(value=True)

        # Gene name lookup
        self.gene_name_lookup_var = tk.BooleanVar(value=True)
        self.species_var = tk.StringVar(value="human")

        # Genes of interest for labeled volcano plots
        self.genes_of_interest_var = tk.StringVar(value="MIAT, QKI, QKI-5, QKI-6, QKI-7")

        # GSEA database selection
        self.gsea_db_vars = {}
        for db in ["GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                    "GO_Molecular_Function_2023", "KEGG_2021_Human",
                    "Reactome_2022", "MSigDB_Hallmark_2020", "WikiPathway_2021_Human"]:
            self.gsea_db_vars[db] = tk.BooleanVar(value=True)

        # GO ORA database selection
        self.ora_db_vars = {}
        for db in ["GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                    "GO_Molecular_Function_2023", "KEGG_2021_Human", "Reactome_2022"]:
            self.ora_db_vars[db] = tk.BooleanVar(value=True)

        # Column-name variables (keyed by "deseq2_<key>" and "rmats_<key>")
        self.col_vars: dict[str, tk.StringVar] = {}
        for k, v in _DESEQ2_COL_DEFAULTS.items():
            self.col_vars[f"deseq2_{k}"] = tk.StringVar(value=v)
        for k, v in _RMATS_COL_DEFAULTS.items():
            self.col_vars[f"rmats_{k}"] = tk.StringVar(value=v)

        self._build_ui()

        # Redirect stdout / stderr to the log widget
        redir = _TextRedirector(self.log_widget, self)
        sys.stdout = redir
        sys.stderr = redir

        # Start with one blank condition row
        self._add_condition_row()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        # Top: notebook tabs
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(side="top", fill="both", expand=True, padx=10, pady=(10, 4))

        self._build_conditions_tab()
        self._build_thresholds_tab()
        self._build_output_tab()
        self._build_enrichment_tab()
        self._build_columns_tab()

        # Middle: run controls
        ctrl_frame = ttk.Frame(self)
        ctrl_frame.pack(side="top", fill="x", padx=10, pady=4)

        self.run_btn = ttk.Button(ctrl_frame, text="▶  Run Pipeline",
                                  command=self._run_pipeline, width=22)
        self.run_btn.pack(side="left", padx=(0, 12))

        self.status_lbl = ttk.Label(ctrl_frame, text="Status: Ready", foreground="#444444")
        self.status_lbl.pack(side="left")

        # Bottom: log output
        log_frame = ttk.LabelFrame(self, text="Pipeline Log")
        log_frame.pack(side="top", fill="both", expand=False, padx=10, pady=(0, 10))

        self.log_widget = scrolledtext.ScrolledText(
            log_frame, height=14, state="disabled",
            wrap="word", font=("Consolas", 9), background="#1e1e1e", foreground="#d4d4d4")
        self.log_widget.pack(fill="both", expand=True, padx=4, pady=4)

    # ---- Tab 1: Conditions ------------------------------------------------

    def _build_conditions_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="  Conditions  ")

        btn_bar = ttk.Frame(tab)
        btn_bar.pack(side="top", fill="x", padx=8, pady=(8, 4))
        self.add_btn = ttk.Button(btn_bar, text="＋  Add Condition",
                                  command=self._add_condition_row, width=18)
        self.add_btn.pack(side="left", padx=(0, 8))
        self.remove_btn = ttk.Button(btn_bar, text="－  Remove Last",
                                     command=self._remove_last_condition, width=18)
        self.remove_btn.pack(side="left")
        ttk.Label(btn_bar, text=f"(maximum {MAX_CONDITIONS} conditions)",
                  foreground="#888888").pack(side="left", padx=12)

        # Scrollable canvas for condition rows
        canvas_frame = ttk.Frame(tab)
        canvas_frame.pack(fill="both", expand=True, padx=8, pady=4)

        self._cond_canvas = tk.Canvas(canvas_frame, borderwidth=0, highlightthickness=0)
        scrollbar = ttk.Scrollbar(canvas_frame, orient="vertical",
                                  command=self._cond_canvas.yview)
        self._cond_canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self._cond_canvas.pack(side="left", fill="both", expand=True)

        self._cond_inner = ttk.Frame(self._cond_canvas)
        self._cond_canvas_window = self._cond_canvas.create_window(
            (0, 0), window=self._cond_inner, anchor="nw")

        self._cond_inner.bind("<Configure>", self._on_cond_frame_configure)
        self._cond_canvas.bind("<Configure>", self._on_canvas_configure)

    def _on_cond_frame_configure(self, event=None):
        self._cond_canvas.configure(scrollregion=self._cond_canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self._cond_canvas.itemconfig(self._cond_canvas_window, width=event.width)

    def _add_condition_row(self):
        if len(self.condition_rows) >= MAX_CONDITIONS:
            return

        idx = len(self.condition_rows)
        frame = ttk.LabelFrame(self._cond_inner,
                               text=f"  Condition {idx + 1}  ",
                               padding=(10, 6))
        frame.pack(fill="x", padx=6, pady=6)

        # Row 0: name + label
        r0 = ttk.Frame(frame)
        r0.pack(fill="x", pady=2)
        ttk.Label(r0, text="Short Name:", width=14, anchor="e").pack(side="left")
        name_var = tk.StringVar()
        ttk.Entry(r0, textvariable=name_var, width=22).pack(side="left", padx=(4, 16))
        ttk.Label(r0, text="Display Label:", width=14, anchor="e").pack(side="left")
        label_var = tk.StringVar()
        ttk.Entry(r0, textvariable=label_var, width=30).pack(side="left", padx=4)

        # Row 1: DESeq2 file
        r1 = ttk.Frame(frame)
        r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="DESeq2 File:", width=14, anchor="e").pack(side="left")
        deseq2_var = tk.StringVar()
        ttk.Entry(r1, textvariable=deseq2_var, width=64).pack(side="left", padx=(4, 6))
        ttk.Button(r1, text="Browse…",
                   command=lambda i=idx: self._browse_deseq2_file(i),
                   width=9).pack(side="left")

        # Row 2: rMATS dir
        r2 = ttk.Frame(frame)
        r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="rMATS Dir:", width=14, anchor="e").pack(side="left")
        rmats_var = tk.StringVar()
        rmats_none_var = tk.BooleanVar(value=False)
        rmats_entry = ttk.Entry(r2, textvariable=rmats_var, width=64)
        rmats_entry.pack(side="left", padx=(4, 6))
        rmats_btn = ttk.Button(r2, text="Browse…",
                               command=lambda i=idx: self._browse_rmats_dir(i),
                               width=9)
        rmats_btn.pack(side="left", padx=(0, 8))
        ttk.Checkbutton(r2, text="None (no rMATS data)",
                        variable=rmats_none_var,
                        command=lambda i=idx: self._toggle_rmats_none(i)
                        ).pack(side="left")

        row_data = {
            "frame":          frame,
            "name_var":       name_var,
            "label_var":      label_var,
            "deseq2_var":     deseq2_var,
            "rmats_var":      rmats_var,
            "rmats_none_var": rmats_none_var,
            "rmats_entry":    rmats_entry,
            "rmats_btn":      rmats_btn,
        }
        self.condition_rows.append(row_data)
        self._update_add_btn_state()

    def _remove_last_condition(self):
        if not self.condition_rows:
            return
        row = self.condition_rows.pop()
        row["frame"].destroy()
        self._update_add_btn_state()

    def _update_add_btn_state(self):
        state = "disabled" if len(self.condition_rows) >= MAX_CONDITIONS else "normal"
        self.add_btn.configure(state=state)
        rem_state = "disabled" if len(self.condition_rows) <= 1 else "normal"
        self.remove_btn.configure(state=rem_state)

    def _browse_deseq2_file(self, idx):
        path = filedialog.askopenfilename(
            title="Select DESeq2 results file",
            filetypes=[("Excel / CSV / TSV", "*.xlsx *.xls *.csv *.tsv *.txt"),
                       ("All files", "*.*")])
        if path:
            self.condition_rows[idx]["deseq2_var"].set(path)

    def _browse_rmats_dir(self, idx):
        path = filedialog.askdirectory(title="Select rMATS output directory")
        if path:
            self.condition_rows[idx]["rmats_var"].set(path)

    def _toggle_rmats_none(self, idx):
        row = self.condition_rows[idx]
        is_none = row["rmats_none_var"].get()
        state = "disabled" if is_none else "normal"
        row["rmats_entry"].configure(state=state)
        row["rmats_btn"].configure(state=state)
        if is_none:
            row["rmats_var"].set("")

    # ---- Tab 2: Thresholds ------------------------------------------------

    def _build_thresholds_tab(self):
        tab = ttk.Frame(self.notebook, padding=12)
        self.notebook.add(tab, text="  Thresholds  ")

        # DESeq2
        de_frame = ttk.LabelFrame(tab, text="  DESeq2 Cutoffs  ", padding=(12, 8))
        de_frame.pack(fill="x", pady=(0, 12))
        self._labeled_entry(de_frame, "| log2FC | threshold:", self.log2fc_var)
        self._labeled_entry(de_frame, "Minimum baseMean:", self.basemean_var)
        self._labeled_entry(de_frame, "Adjusted p-value (padj):", self.padj_var)
        ttk.Checkbutton(de_frame, text="Auto biotype split (protein-coding vs non-coding)",
                        variable=self.auto_bio_var).grid(
            row=3, column=0, columnspan=3, sticky="w", pady=4)

        # rMATS
        rm_frame = ttk.LabelFrame(tab, text="  rMATS Cutoffs  ", padding=(12, 8))
        rm_frame.pack(fill="x")
        self._labeled_entry(rm_frame, "FDR threshold:", self.rmats_fdr_var)
        self._labeled_entry(rm_frame, "P-value threshold:", self.rmats_pval_var)
        self._labeled_entry(rm_frame, "| ΔΨ | (IncLevelDifference):", self.dpsi_var)

        fdr_row = ttk.Frame(rm_frame)
        fdr_row.grid(row=3, column=0, columnspan=3, sticky="w", pady=6)
        ttk.Label(fdr_row, text="Filter by:", width=26, anchor="e").pack(side="left")
        ttk.Radiobutton(fdr_row, text="FDR", variable=self.use_fdr_var,
                        value=True).pack(side="left", padx=(8, 4))
        ttk.Radiobutton(fdr_row, text="P-value", variable=self.use_fdr_var,
                        value=False).pack(side="left")

    def _labeled_entry(self, parent, label_text, var, row=None):
        if row is None:
            row = parent.grid_size()[1]
        ttk.Label(parent, text=label_text, width=30, anchor="e").grid(
            row=row, column=0, sticky="e", padx=(0, 6), pady=3)
        ttk.Entry(parent, textvariable=var, width=12).grid(
            row=row, column=1, sticky="w", pady=3)

    # ---- Tab 3: Output & Figures ------------------------------------------

    def _build_output_tab(self):
        tab = ttk.Frame(self.notebook, padding=12)
        self.notebook.add(tab, text="  Output & Figures  ")

        out_frame = ttk.LabelFrame(tab, text="  Output  ", padding=(12, 8))
        out_frame.pack(fill="x", pady=(0, 12))

        ttk.Label(out_frame, text="Output Directory:", width=20, anchor="e").grid(
            row=0, column=0, sticky="e", padx=(0, 6), pady=4)
        ttk.Entry(out_frame, textvariable=self.output_dir_var, width=60).grid(
            row=0, column=1, sticky="ew", pady=4)
        ttk.Button(out_frame, text="Browse…",
                   command=self._browse_output_dir, width=9).grid(
            row=0, column=2, padx=6, pady=4)
        out_frame.columnconfigure(1, weight=1)

        fig_frame = ttk.LabelFrame(tab, text="  Figure Settings  ", padding=(12, 8))
        fig_frame.pack(fill="x", pady=(0, 12))

        # Format
        ttk.Label(fig_frame, text="Format:", width=20, anchor="e").grid(
            row=0, column=0, sticky="e", padx=(0, 6), pady=4)
        ttk.OptionMenu(fig_frame, self.fig_format_var, "png", "png", "svg", "pdf").grid(
            row=0, column=1, sticky="w", pady=4)

        self._labeled_entry(fig_frame, "DPI:", self.dpi_var, row=1)
        self._labeled_entry(fig_frame, "Font size:", self.font_size_var, row=2)

        # Colors
        for row_idx, (label, var_name) in enumerate([
            ("Color — Up-regulated:",   "color_up_var"),
            ("Color — Down-regulated:", "color_down_var"),
            ("Color — Not significant:", "color_ns_var"),
        ], start=3):
            var = getattr(self, var_name)
            ttk.Label(fig_frame, text=label, width=24, anchor="e").grid(
                row=row_idx, column=0, sticky="e", padx=(0, 6), pady=4)
            entry = ttk.Entry(fig_frame, textvariable=var, width=10)
            entry.grid(row=row_idx, column=1, sticky="w", pady=4)
            swatch = tk.Button(fig_frame, bg=var.get(), width=3, relief="groove",
                               command=lambda v=var, s=None: self._pick_color(v))
            swatch.grid(row=row_idx, column=2, padx=6, pady=4)
            # Keep swatch in sync with the entry
            var.trace_add("write", lambda *_, v=var, s=swatch: self._sync_swatch(v, s))
            # Store swatch reference in the variable so trace can reach it
            swatch._color_var = var

        ttk.Checkbutton(fig_frame,
                        text="Generate interactive HTML plots (requires plotly)",
                        variable=self.interactive_var).grid(
            row=6, column=0, columnspan=3, sticky="w", pady=6)

    def _browse_output_dir(self):
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.output_dir_var.set(path)

    def _pick_color(self, var):
        current = var.get()
        result = colorchooser.askcolor(color=current, title="Choose colour")
        if result and result[1]:
            var.set(result[1].upper())

    def _sync_swatch(self, var, swatch):
        try:
            swatch.configure(bg=var.get())
        except tk.TclError:
            pass

    # ---- Tab 4: Enrichment ------------------------------------------------

    def _build_enrichment_tab(self):
        tab = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(tab, text="Enrichment")

        # --- GSEA Databases ---
        gsea_frame = ttk.LabelFrame(tab, text="GSEA Pre-ranked Databases", padding=8)
        gsea_frame.pack(fill="x", pady=(0, 8))

        _DB_DISPLAY = {
            "GO_Biological_Process_2023": "GO Biological Process",
            "GO_Cellular_Component_2023": "GO Cellular Component",
            "GO_Molecular_Function_2023": "GO Molecular Function",
            "KEGG_2021_Human": "KEGG (Human)",
            "Reactome_2022": "Reactome",
            "MSigDB_Hallmark_2020": "MSigDB Hallmark",
            "WikiPathway_2021_Human": "WikiPathway (Human)",
        }

        for db, var in self.gsea_db_vars.items():
            ttk.Checkbutton(gsea_frame, text=_DB_DISPLAY.get(db, db), variable=var).pack(anchor="w")

        btn_frame = ttk.Frame(gsea_frame)
        btn_frame.pack(fill="x", pady=(4, 0))
        ttk.Button(btn_frame, text="Select All",
                   command=lambda: [v.set(True) for v in self.gsea_db_vars.values()]).pack(side="left", padx=2)
        ttk.Button(btn_frame, text="Deselect All",
                   command=lambda: [v.set(False) for v in self.gsea_db_vars.values()]).pack(side="left", padx=2)

        # --- ORA Databases ---
        ora_frame = ttk.LabelFrame(tab, text="GO / Pathway ORA Databases", padding=8)
        ora_frame.pack(fill="x", pady=(0, 8))

        for db, var in self.ora_db_vars.items():
            ttk.Checkbutton(ora_frame, text=_DB_DISPLAY.get(db, db), variable=var).pack(anchor="w")

        btn_frame2 = ttk.Frame(ora_frame)
        btn_frame2.pack(fill="x", pady=(4, 0))
        ttk.Button(btn_frame2, text="Select All",
                   command=lambda: [v.set(True) for v in self.ora_db_vars.values()]).pack(side="left", padx=2)
        ttk.Button(btn_frame2, text="Deselect All",
                   command=lambda: [v.set(False) for v in self.ora_db_vars.values()]).pack(side="left", padx=2)

        # --- General Enrichment Settings ---
        gen_frame = ttk.LabelFrame(tab, text="General Settings", padding=8)
        gen_frame.pack(fill="x", pady=(0, 8))

        row = ttk.Frame(gen_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Species:").pack(side="left")
        ttk.Combobox(row, textvariable=self.species_var, width=15,
                     values=["human", "mouse", "rat", "zebrafish", "fly", "worm"],
                     state="readonly").pack(side="left", padx=(6, 0))

        ttk.Checkbutton(gen_frame, text="Auto-lookup gene names from Ensembl IDs (requires internet)",
                        variable=self.gene_name_lookup_var).pack(anchor="w", pady=2)

        row2 = ttk.Frame(gen_frame)
        row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Genes of interest:").pack(side="left")
        ttk.Entry(row2, textvariable=self.genes_of_interest_var, width=50).pack(side="left", padx=(6, 0), fill="x", expand=True)

        ttk.Label(gen_frame, text="Comma-separated gene symbols to highlight on volcano plots",
                  foreground="gray").pack(anchor="w")

    # ---- Tab 5: Column Names (Advanced) -----------------------------------

    def _build_columns_tab(self):
        tab = ttk.Frame(self.notebook, padding=12)
        self.notebook.add(tab, text="  Column Names (Advanced)  ")

        de_frame = ttk.LabelFrame(tab, text="  DESeq2 Column Mappings  ", padding=(12, 8))
        de_frame.pack(fill="x", pady=(0, 12))
        for row_idx, (key, lbl) in enumerate(_DESEQ2_COL_LABELS.items()):
            ttk.Label(de_frame, text=lbl + ":", width=30, anchor="e").grid(
                row=row_idx, column=0, sticky="e", padx=(0, 6), pady=3)
            ttk.Entry(de_frame, textvariable=self.col_vars[f"deseq2_{key}"],
                      width=30).grid(row=row_idx, column=1, sticky="w", pady=3)

        rm_frame = ttk.LabelFrame(tab, text="  rMATS Column Mappings  ", padding=(12, 8))
        rm_frame.pack(fill="x", pady=(0, 12))
        for row_idx, (key, lbl) in enumerate(_RMATS_COL_LABELS.items()):
            ttk.Label(rm_frame, text=lbl + ":", width=30, anchor="e").grid(
                row=row_idx, column=0, sticky="e", padx=(0, 6), pady=3)
            ttk.Entry(rm_frame, textvariable=self.col_vars[f"rmats_{key}"],
                      width=30).grid(row=row_idx, column=1, sticky="w", pady=3)

        ttk.Button(tab, text="Reset to Defaults",
                   command=self._reset_column_defaults, width=18).pack(anchor="w")

    def _reset_column_defaults(self):
        for k, v in _DESEQ2_COL_DEFAULTS.items():
            self.col_vars[f"deseq2_{k}"].set(v)
        for k, v in _RMATS_COL_DEFAULTS.items():
            self.col_vars[f"rmats_{k}"].set(v)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self) -> bool:
        errors = []

        if not self.condition_rows:
            errors.append("Add at least one condition.")

        seen_names = set()
        for i, row in enumerate(self.condition_rows, start=1):
            name  = row["name_var"].get().strip()
            label = row["label_var"].get().strip()
            deseq = row["deseq2_var"].get().strip()
            rmats = row["rmats_var"].get().strip()
            is_none = row["rmats_none_var"].get()

            if not name:
                errors.append(f"Condition {i}: Short Name cannot be empty.")
            elif name in seen_names:
                errors.append(f"Condition {i}: Short Name '{name}' is not unique.")
            else:
                seen_names.add(name)

            if not label:
                errors.append(f"Condition {i}: Display Label cannot be empty.")

            if not deseq:
                errors.append(f"Condition {i}: DESeq2 file path is empty.")
            elif not Path(deseq).is_file():
                errors.append(f"Condition {i}: DESeq2 file not found:\n  {deseq}")

            if not is_none and rmats and not Path(rmats).is_dir():
                errors.append(f"Condition {i}: rMATS directory not found:\n  {rmats}")

        numeric_fields = [
            (self.log2fc_var,    "| log2FC | threshold"),
            (self.basemean_var,  "Minimum baseMean"),
            (self.padj_var,      "Adjusted p-value"),
            (self.rmats_fdr_var, "rMATS FDR threshold"),
            (self.rmats_pval_var,"rMATS P-value threshold"),
            (self.dpsi_var,      "| ΔΨ | threshold"),
            (self.dpi_var,       "DPI"),
            (self.font_size_var, "Font size"),
        ]
        for var, lbl in numeric_fields:
            try:
                float(var.get())
            except ValueError:
                errors.append(f"'{lbl}' must be a valid number.")

        if not self.output_dir_var.get().strip():
            errors.append("Output directory cannot be empty.")

        if errors:
            messagebox.showerror("Validation Errors", "\n\n".join(errors))
            return False
        return True

    # ------------------------------------------------------------------
    # Config collection
    # ------------------------------------------------------------------

    def _collect_config(self) -> dict:
        conditions = []
        for row in self.condition_rows:
            is_none = row["rmats_none_var"].get()
            rmats   = None if is_none else (row["rmats_var"].get().strip() or None)
            conditions.append({
                "name":        row["name_var"].get().strip(),
                "label":       row["label_var"].get().strip(),
                "deseq2_file": row["deseq2_var"].get().strip(),
                "rmats_dir":   rmats,
            })

        return {
            "CONDITIONS":           conditions,
            "OUTPUT_DIR":           self.output_dir_var.get().strip(),
            "LOG2FC_CUTOFF":        float(self.log2fc_var.get()),
            "BASEMEAN_CUTOFF":      float(self.basemean_var.get()),
            "PADJ_CUTOFF":          float(self.padj_var.get()),
            "AUTO_BIOTYPE_SPLIT":   self.auto_bio_var.get(),
            "RMATS_FDR_CUTOFF":     float(self.rmats_fdr_var.get()),
            "RMATS_PVAL_CUTOFF":    float(self.rmats_pval_var.get()),
            "INCLEVEL_DIFF_CUTOFF": float(self.dpsi_var.get()),
            "USE_FDR":              self.use_fdr_var.get(),
            "FIG_DPI":              int(float(self.dpi_var.get())),
            "FIG_FORMAT":           self.fig_format_var.get(),
            "FONT_SIZE":            int(float(self.font_size_var.get())),
            "COLOR_UP":             self.color_up_var.get(),
            "COLOR_DOWN":           self.color_down_var.get(),
            "COLOR_NS":             self.color_ns_var.get(),
            "INTERACTIVE_PLOTS":    self.interactive_var.get(),
            "GENE_NAME_LOOKUP":     self.gene_name_lookup_var.get(),
            "SPECIES":              self.species_var.get(),
            "GENES_OF_INTEREST":    [g.strip() for g in self.genes_of_interest_var.get().split(",") if g.strip()],
            "GSEA_DATABASES":       [db for db, var in self.gsea_db_vars.items() if var.get()],
            "ORA_DATABASES":        [db for db, var in self.ora_db_vars.items() if var.get()],
            "DESEQ2_COLS": {k: self.col_vars[f"deseq2_{k}"].get()
                            for k in _DESEQ2_COL_DEFAULTS},
            "RMATS_COLS":  {k: self.col_vars[f"rmats_{k}"].get()
                            for k in _RMATS_COL_DEFAULTS},
        }

    # ------------------------------------------------------------------
    # Pipeline execution
    # ------------------------------------------------------------------

    def _run_pipeline(self):
        if not self._validate():
            return
        config = self._collect_config()
        self.run_btn.configure(state="disabled", text="⏳  Running…")
        self.status_lbl.configure(text="Status: Running…", foreground="#0066cc")
        self._log("\n" + "=" * 60)
        self._log("  Starting pipeline…")
        self._log("=" * 60 + "\n")
        t = threading.Thread(target=self._run_in_thread, args=(config,), daemon=True)
        t.start()

    def _run_in_thread(self, config):
        try:
            pipeline.run_pipeline(config)
            self.after(0, self._on_pipeline_success)
        except Exception:
            msg = traceback.format_exc()
            self.after(0, self._on_pipeline_error, msg)

    def _on_pipeline_success(self):
        self.status_lbl.configure(text="Status: Completed ✓", foreground="#008800")
        self.run_btn.configure(state="normal", text="▶  Run Pipeline")
        self._log("\n✓ Pipeline completed successfully.\n")

    def _on_pipeline_error(self, msg):
        self.status_lbl.configure(text="Status: Error ✗", foreground="#cc0000")
        self.run_btn.configure(state="normal", text="▶  Run Pipeline")
        self._log(f"\n✗ Pipeline error:\n{msg}\n")

    def _log(self, text):
        self.log_widget.configure(state="normal")
        self.log_widget.insert("end", text + "\n")
        self.log_widget.see("end")
        self.log_widget.configure(state="disabled")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app = PipelineLauncherApp()
    app.mainloop()
