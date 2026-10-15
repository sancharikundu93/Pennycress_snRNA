# Pennycress snRNA-seq Explorer -- Shiny app
#
# Type a gene (TAIR symbol, Arabidopsis AGI, or pennycress locus) and see it on the
# UMAP and as violins by cell type, split Control vs Salt.
#
# Reads data/meta.rds + data/expr.bin, written by data_prep/step1 and step2.
#   Run locally:  shiny::runApp("app")
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(data.table)
})

APP_TITLE <- "Pennycress snRNA-seq Explorer"
D <- readRDS("data/meta.rds")
cells <- D$cells; genes <- D$genes; N_CELLS <- D$n_cells
stopifnot(nrow(cells) == N_CELLS)
CT <- D$ct_levels; CT_COL <- D$ct_colours[CT]
cells[, cell_type := factor(cell_type, levels = CT)]
cells[, condition := factor(condition, levels = c("Control", "Salt"))]
COND_COL <- c(Control = "#4c78a8", Salt = "#e45756")
EXPR_SCALE <- list(c(0, "#d9d9d9"), c(0.15, "#fdd49e"), c(0.5, "#ef6548"), c(1, "#7f0000"))

# TAIR lists several aliases per gene ("ATSUC2;SUC2;SUT1"), so every alias goes in the
# label -- the search box only matches what the label shows. An "AT"-prefixed alias is
# dropped when the bare name is also there (ATSUC2 next to SUC2 adds nothing).
tidy_syms <- function(s) {
  if (is.na(s)) return(NA_character_)
  a <- unique(strsplit(s, ";")[[1]])
  a <- a[!(grepl("^AT", a) & sub("^AT", "", a) %in% a)]
  paste(a, collapse = ", ")
}
genes[, names := vapply(symbols, tidy_syms, character(1))]
# One searchable label per gene: "SUC2, SUT1 | AT1G22710 | TaMN106.1G029100"
genes[, label := paste(c(if (!is.na(names)) names, if (!is.na(agi)) agi, locus), collapse = " | "),
      by = feature]
GENE_CHOICES <- setNames(genes$feature, genes$label)
DEFAULT_GENE <- genes[grepl("(^|, )SUC2(,|$)", names), feature][1]
if (is.na(DEFAULT_GENE)) DEFAULT_GENE <- genes$feature[1]

# expr.bin holds one block per gene: int32 cell indices, then int16 values x 1000.
# Seek to the gene's block so only that gene is ever in memory.
gene_vec <- function(f) {
  r <- genes[feature == f]
  v <- numeric(N_CELLS)
  if (r$nnz == 0) return(v)
  con <- file("data/expr.bin", "rb"); on.exit(close(con))
  seek(con, r$offset)
  i <- readBin(con, "integer", n = r$nnz, size = 4)
  x <- readBin(con, "integer", n = r$nnz, size = 2, signed = TRUE)
  v[i] <- x / 1000
  v
}
gene_name <- function(f) genes[feature == f, label]

ui <- page_navbar(
  title = tags$span(APP_TITLE, tags$small(" beta", style = "color:#e45756;font-weight:600")),
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 320,
    selectizeInput("gene", "Gene", choices = NULL,
                   options = list(placeholder = "e.g. SUC2, AT1G22710, TaMN106...")),
    helpText("Search by Arabidopsis symbol, AGI, or pennycress locus.",
             "Symbols come from the Arabidopsis best-hit ortholog."),
    radioButtons("split", "UMAP view", inline = TRUE,
                 choices = c("Combined" = "all", "Control | Salt" = "split")),
    sliderInput("pt", "Point size", min = 1, max = 6, value = 3, step = 0.5),
    hr(),
    tags$small(sprintf("%s nuclei · %d cell types · Control vs Salt",
                       format(nrow(cells), big.mark = ","), length(CT)))
  ),
  nav_panel("Gene expression",
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("Cell types"),
           plotlyOutput("umap_ct", height = "460px")),
      card(full_screen = TRUE, card_header(textOutput("gene_title", inline = TRUE)),
           plotlyOutput("umap_gene", height = "460px"))),
    card(full_screen = TRUE,
         card_header("Expression by cell type, Control vs Salt"),
         plotlyOutput("violin", height = "420px"),
         uiOutput("pct_note"))
  ),
  nav_panel("3D UMAP",
    card(full_screen = TRUE,
         card_header(radioButtons("col3d", NULL, inline = TRUE,
                                  choices = c("Cell type" = "ct", "Selected gene" = "gene",
                                              "Condition" = "cond"))),
         plotlyOutput("umap3d", height = "720px"))
  ),
  nav_panel("About",
    card(card_body(
      h4(em("Thlaspi arvense"), " (pennycress) root single-nucleus RNA-seq"),
      p("Beta release of preliminary, unpublished data. Labels and layout may change."),
      tags$ul(
        tags$li(sprintf("%s nuclei after QC (nCount_RNA > 500): Control %s, Salt %s.",
                        format(nrow(cells), big.mark = ","),
                        format(sum(cells$condition == "Control"), big.mark = ","),
                        format(sum(cells$condition == "Salt"), big.mark = ","))),
        tags$li("Clusters at Seurat resolution 0.6, annotated with marker genes."),
        tags$li("Expression values are log-normalised counts (Seurat LogNormalize)."),
        tags$li("Gene symbols are from the Arabidopsis best-hit ortholog (TAIR); one",
                "symbol can match several pennycress loci, each listed separately.")),
      p("Code: ", tags$a(href = "https://github.com/sancharikundu93/Pennycress_snRNA",
                         "github.com/sancharikundu93/Pennycress_snRNA")),
      p(tags$small("Data built ", D$made))
    ))
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = GENE_CHOICES,
                       selected = DEFAULT_GENE, server = TRUE)

  g <- reactive({ req(input$gene, input$gene %in% genes$feature); input$gene })
  gv <- reactive(gene_vec(g()))

  output$gene_title <- renderText(paste("Expression:", gene_name(g())))

  ax <- list(title = "", showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE)

  output$umap_ct <- renderPlotly({
    plot_ly(cells, x = ~u1, y = ~u2, color = ~cell_type, colors = CT_COL,
            type = "scattergl", mode = "markers",
            marker = list(size = input$pt), text = ~cell_type, hoverinfo = "text") |>
      layout(xaxis = ax, yaxis = c(ax, list(scaleanchor = "x")),
             legend = list(itemsizing = "constant", font = list(size = 10))) |>
      config(displaylogo = FALSE)
  })

  gene_umap <- function(d, showscale = TRUE) {
    d <- d[order(d$e), ]                       # expressing nuclei drawn on top
    plot_ly(d, x = ~u1, y = ~u2, type = "scattergl", mode = "markers",
            marker = list(size = input$pt, color = ~e, colorscale = EXPR_SCALE,
                          cmin = 0, cmax = max(1e-6, max(gv())),
                          showscale = showscale,
                          colorbar = list(title = "log-norm", len = 0.6)),
            text = ~sprintf("%s<br>%s<br>%.2f", cell_type, condition, e),
            hoverinfo = "text", showlegend = FALSE) |>
      layout(xaxis = ax, yaxis = c(ax, list(scaleanchor = "x")))
  }

  output$umap_gene <- renderPlotly({
    d <- data.frame(cells[, .(u1, u2, cell_type, condition)], e = gv())
    p <- if (input$split == "all") gene_umap(d) else
      subplot(gene_umap(d[d$condition == "Control", ], FALSE) |>
                layout(annotations = list(text = "Control", x = 0.5, y = 1.05,
                                          xref = "paper", yref = "paper", showarrow = FALSE)),
              gene_umap(d[d$condition == "Salt", ], TRUE) |>
                layout(annotations = list(text = "Salt", x = 0.5, y = 1.05,
                                          xref = "paper", yref = "paper", showarrow = FALSE)),
              shareX = TRUE, shareY = TRUE)
    config(p, displaylogo = FALSE)
  })

  output$violin <- renderPlotly({
    d <- data.frame(cell_type = cells$cell_type, condition = cells$condition, e = gv())
    plot_ly(d, x = ~cell_type, y = ~e, color = ~condition, colors = COND_COL,
            type = "violin", box = list(visible = TRUE), meanline = list(visible = FALSE),
            points = FALSE, spanmode = "hard", scalemode = "width") |>
      layout(violinmode = "group", xaxis = list(title = "", tickangle = -40),
             yaxis = list(title = "log-normalised expression", rangemode = "tozero"),
             legend = list(orientation = "h", y = 1.1)) |>
      config(displaylogo = FALSE)
  })

  output$pct_note <- renderUI({
    e <- gv()
    pc <- tapply(e > 0, cells$condition, mean) * 100
    tags$small(sprintf("Detected in %.1f%% of Control and %.1f%% of Salt nuclei. Violins include non-expressing nuclei (zeros).",
                       pc["Control"], pc["Salt"]), style = "color:#666;padding:0 12px")
  })

  output$umap3d <- renderPlotly({
    m <- list(size = max(1, input$pt - 1.5))
    p <- switch(input$col3d,
      ct   = plot_ly(cells, x = ~x3, y = ~y3, z = ~z3, color = ~cell_type, colors = CT_COL,
                     type = "scatter3d", mode = "markers", marker = m,
                     text = ~cell_type, hoverinfo = "text"),
      cond = plot_ly(cells, x = ~x3, y = ~y3, z = ~z3, color = ~condition, colors = COND_COL,
                     type = "scatter3d", mode = "markers", marker = m,
                     text = ~cell_type, hoverinfo = "text"),
      gene = {
        d <- data.frame(cells[, .(x3, y3, z3, cell_type)], e = gv())
        plot_ly(d, x = ~x3, y = ~y3, z = ~z3, type = "scatter3d", mode = "markers",
                marker = c(m, list(color = ~e, colorscale = EXPR_SCALE, showscale = TRUE,
                                   colorbar = list(title = gene_name(g())))),
                text = ~sprintf("%s<br>%.2f", cell_type, e), hoverinfo = "text")
      })
    p |> layout(scene = list(xaxis = list(title = "UMAP 1"), yaxis = list(title = "UMAP 2"),
                             zaxis = list(title = "UMAP 3"), aspectmode = "data"),
                legend = list(itemsizing = "constant")) |>
      config(displaylogo = FALSE)
  })
}

shinyApp(ui, server)
