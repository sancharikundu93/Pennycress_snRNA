# Pennycress snRNA-seq Explorer (beta)

### 👉 [Open the explorer](https://SHINYAPPS_ACCOUNT.shinyapps.io/pennycress_snrna/)

> **Beta.** Preliminary, unpublished data. Clusters, labels and layout may change.

Interactive viewer for *Thlaspi arvense* (pennycress) root single-nucleus RNA-seq,
Control vs Salt. Type a gene name and see:

- its expression on the **UMAP** (combined, or Control and Salt side by side)
- **violin plots** by cell type, split **Control vs Salt**
- the **3D UMAP**, coloured by cell type, condition or the selected gene

Genes can be searched by Arabidopsis symbol (e.g. `SUC2`), AGI (`AT1G22710`) or
pennycress locus (`TaMN106.1G…`). Symbols come from the Arabidopsis best-hit ortholog.

## Dataset

- 25,206 nuclei after QC (nCount_RNA > 500): Control 13,362 · Salt 11,844
- 23 clusters (Seurat resolution 0.6), annotated with marker genes
- Reference genome: *T. arvense* MN106 v4.1 (Phytozome)

## Repository layout

| Folder | Contents |
|---|---|
| `app/` | The Shiny app (`app.R`). Its data file is built by `data_prep/` and not stored in git. |
| `data_prep/` | `step1_export_app_data.R`: builds the app's data file from the integrated Seurat object |
| `analysis_scripts/` | Analysis code behind the UMAPs and gene expression figures (clustering, annotation, feature plots, violins, dot plots, 3D UMAP) |
| `docs/` | Stand-alone 3D UMAP pages (single HTML files, work offline) |

## Run the app locally

```r
# after building app/data/app_data.rds with data_prep/step1_export_app_data.R
shiny::runApp("app")
```

## Contact

Sanchari Kundu: [github.com/sancharikundu93](https://github.com/sancharikundu93)
