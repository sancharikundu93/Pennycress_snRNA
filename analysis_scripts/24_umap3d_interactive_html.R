#!/usr/bin/env Rscript
# Step 24 -- ONE self-contained, auto-rotating HTML of the 3D UMAP.
#
# WHY THIS EXISTS. Step 23's htmlwidgets output came up blank when opened. Two
# reasons, both fixed here:
#   1. saveWidget(selfcontained = TRUE) needs pandoc, which is not installed on the
#      compute nodes, so the page was written with selfcontained = FALSE and depends
#      on a sibling umap3d_interactive_files/ directory. Move or copy the .html
#      alone -- as anyone would -- and the page loads with nothing in it.
#   2. It only rotated if you dragged it. "3D" here means it should turn on its own.
#
# This writes a single HTML file with plotly.js (3.5 MB) inlined from the local R
# package and the coordinates embedded as JS arrays. No sibling directory, no CDN, no
# network: it works from a laptop with no internet, emailed as one attachment.
# The camera orbits on a timer, with a button to pause; dragging still works and
# pauses the spin so you can inspect an angle.
#
# Reads results/23_umap3d/umap3d_embedding.csv, so it does NOT reload the 5 GB object
# or re-run UMAP -- it is seconds, and the coordinates are identical to the ones the
# static panels were drawn from.
#
#   Rscript 24_umap3d_interactive_html.R
#   POINT_SIZE=2.2 Rscript 24_umap3d_interactive_html.R
#
# Outputs -> ncount500_08132026/results/23_umap3d/umap3d_rotating.html
suppressPackageStartupMessages({ library(data.table) })

BASE <- "/projects/songli_lab/A_Sanchari/A_pennycress_snRNA"
ROOT <- file.path(BASE, "ncount500_08132026")
OUT  <- file.path(ROOT, "results/23_umap3d")
CSV  <- file.path(OUT, "umap3d_embedding.csv")
source(file.path(ROOT, "scripts/palette_celltypes.R"))
stopifnot(file.exists(CSV))

PLOTLY_JS <- "/home/sancharikundu/R/owl-genoa/4.4.2/plotly/htmlwidgets/lib/plotlyjs/plotly-latest.min.js"
if (!file.exists(PLOTLY_JS)) stop("plotly.js not found at ", PLOTLY_JS)

PT   <- as.numeric(Sys.getenv("POINT_SIZE", "1.8"))
SPIN <- as.numeric(Sys.getenv("SPIN_MS", "60"))     # ms between camera steps

# THEME=dark is the default, chosen 2026-08-13. Points keep their PAL colours in both
# themes rather than being brightened for the dark one: a cell type has to be the same
# colour here as in the UMAP and the proportion bars, and a palette that shifts with
# the background stops being a shared palette. The dark greys below are for the page,
# axes and grid only.
THEME <- Sys.getenv("THEME", "dark")
if (!THEME %in% c("dark", "light")) stop("THEME must be 'dark' or 'light', got: ", THEME)
DARK <- identical(THEME, "dark")

BG      <- if (DARK) "#000000" else "#ffffff"   # page and plot ground
FG      <- if (DARK) "#e8e8e8" else "#111111"   # body text, axis titles
DIM_FG  <- if (DARK) "#9a9a9a" else "#555555"   # caption text
GRID    <- if (DARK) "#2e2e2e" else "#e0e0e0"   # scene gridlines
RULE    <- if (DARK) "#333333" else "#dddddd"   # header rule, button border
BTN_BG  <- if (DARK) "#1c1c1c" else "#f6f6f6"
# On black, 0.8 opacity over a dark ground reads dimmer than the same value on white.
OPACITY <- if (DARK) 0.9 else 0.8

d <- fread(CSV)
d[, cell_type := factor(cell_type, levels = CT_LEVELS[CT_LEVELS %in% unique(cell_type)])]
cat("nuclei:", nrow(d), " cell types:", nlevels(d$cell_type), "\n")

# 3 decimals is well past what a screen can resolve and keeps the file small
r3 <- function(v) formatC(round(v, 3), format = "f", digits = 3, drop0trailing = TRUE)
# JSON string literal. Cell-type names here are plain ASCII (the underscore in
# pericycle_like is the only oddity), but quote defensively rather than assume.
jq <- function(s) paste0('"', gsub('"', '\\\\"', s), '"')
PAL_CT <- ct_palette(levels(d$cell_type))

# One trace per cell type, so the legend toggles cell types on and off -- clicking a
# legend entry to isolate a population is the main thing this view is good for.
traces <- vapply(levels(d$cell_type), function(ct) {
  s <- d[cell_type == ct]
  sprintf('{"type":"scatter3d","mode":"markers","name":%s,
  "x":[%s],"y":[%s],"z":[%s],
  "marker":{"size":%s,"color":%s,"opacity":%s,"line":{"width":0}},
  "hoverinfo":"name"}',
          jq(ct), paste(r3(s$x), collapse = ","), paste(r3(s$y), collapse = ","),
          paste(r3(s$z), collapse = ","), PT, jq(unname(PAL_CT[ct])), OPACITY)
}, character(1))

html <- sprintf('<!doctype html>
<html><head><meta charset="utf-8"><title>3D UMAP - pennycress root, res 0.6</title>
<style>
  body{margin:0;font-family:Georgia,"Times New Roman",serif;background:%s;color:%s}
  #bar{padding:10px 14px;border-bottom:1px solid %s}
  #bar h1{font-size:16px;margin:0 0 4px}
  #bar p{font-size:12px;margin:0;color:%s}
  button{font-family:inherit;font-size:13px;padding:4px 12px;margin-right:8px;
         border:1px solid %s;background:%s;color:%s;border-radius:3px;cursor:pointer}
  #plot{width:100vw;height:calc(100vh - 96px)}
</style>
<script>%s</script></head>
<body>
<div id="bar">
  <h1>3D UMAP &mdash; <i>Thlaspi arvense</i> root, resolution 0.6, nCount_RNA &gt; 500</h1>
  <p>%s nuclei. Drag to rotate, scroll to zoom, click a legend entry to hide a cell type.
     A separate n.components = 3 embedding &mdash; not a rotation of the 2D UMAP.</p>
  <div style="margin-top:8px">
    <button id="toggle">Pause rotation</button>
    <button id="reset">Reset view</button>
  </div>
</div>
<div id="plot"></div>
<script>
var DATA = [%s];
// showbackground:false removes plotly default grey wall panels, which on black read
// as three lit rectangles floating behind the cloud.
var AX = {color:"%s", gridcolor:"%s", zerolinecolor:"%s", showbackground:false};
var LAYOUT = {
  margin:{l:0,r:0,t:0,b:0}, showlegend:true,
  paper_bgcolor:"%s", plot_bgcolor:"%s",
  legend:{font:{size:11,color:"%s"},itemsizing:"constant"},
  scene:{ xaxis:Object.assign({title:"UMAP 1"}, AX),
          yaxis:Object.assign({title:"UMAP 2"}, AX),
          zaxis:Object.assign({title:"UMAP 3"}, AX),
          aspectmode:"data",
          camera:{eye:{x:1.8,y:0,z:0.6}} }
};
var gd = document.getElementById("plot");
Plotly.newPlot(gd, DATA, LAYOUT, {responsive:true, displaylogo:false});

// --- turntable ------------------------------------------------------------
// Orbit the camera on a timer. Radius and height are taken from the starting eye so
// a manual zoom is preserved when the spin resumes.
var angle = 0, spinning = true, R = 1.8, H = 0.6;
setInterval(function(){
  if(!spinning) return;
  angle += 0.02;
  Plotly.relayout(gd, {"scene.camera.eye":{x:R*Math.cos(angle), y:R*Math.sin(angle), z:H}});
}, %d);

// Dragging while it spins fights the timer, so any manual interaction stops it.
gd.on("plotly_relayouting", function(){ if(spinning) setSpin(false); });

function setSpin(on){
  spinning = on;
  document.getElementById("toggle").textContent = on ? "Pause rotation" : "Resume rotation";
  if(on){
    var e = (gd.layout.scene.camera||{}).eye || {x:1.8,y:0,z:0.6};
    R = Math.sqrt(e.x*e.x + e.y*e.y); H = e.z; angle = Math.atan2(e.y, e.x);
  }
}
document.getElementById("toggle").onclick = function(){ setSpin(!spinning); };
document.getElementById("reset").onclick = function(){
  R = 1.8; H = 0.6; angle = 0;
  Plotly.relayout(gd, {"scene.camera.eye":{x:R,y:0,z:H}});
};
</script></body></html>',
  # order matters: these fill the %s/%d above, top to bottom
  BG, FG,                                    # body background, body text
  RULE, DIM_FG,                              # header rule, caption text
  RULE, BTN_BG, FG,                          # button border, fill, text
  paste(readLines(PLOTLY_JS, warn = FALSE), collapse = "\n"),
  format(nrow(d), big.mark = ","),
  paste(traces, collapse = ",\n"),
  FG, GRID, GRID,                            # axis text, gridlines, zero lines
  BG, BG, FG,                                # paper, plot ground, legend text
  SPIN)

# The dark page is the one to open, so it keeps the plain name; light goes to its own
# file rather than overwriting it.
dest <- file.path(OUT, if (DARK) "umap3d_rotating.html" else "umap3d_rotating_light.html")
writeLines(html, dest)
cat("wrote", dest, "--", round(file.info(dest)$size / 1024^2, 1), "MB\n")
cat("single file: no sibling directory, no network needed\n")
cat("DONE\n")
