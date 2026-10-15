#!/bin/bash
# Upload app/ to shinyapps.io.
#   SHINY_ACCOUNT=... SHINY_TOKEN=... SHINY_SECRET=... ./deploy.sh
# The token and secret come from shinyapps.io -> Account -> Tokens -> Show.
# rsconnect stores them under ~/.config/R/rsconnect after the first run, so later
# deploys only need SHINY_ACCOUNT.
set -eo pipefail
: "${SHINY_ACCOUNT:?set SHINY_ACCOUNT to your shinyapps.io account name}"
module reset >/dev/null 2>&1
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 >/dev/null 2>&1
export R_LIBS_USER=/home/sancharikundu/R/owl-genoa/4.4.2
cd "$(dirname "$0")"
[ -s app/data/meta.rds ] && [ -s app/data/expr.bin ] || { echo "app/data missing: run data_prep step1 + step2 first"; exit 1; }
Rscript -e '
a <- Sys.getenv("SHINY_ACCOUNT"); t <- Sys.getenv("SHINY_TOKEN"); s <- Sys.getenv("SHINY_SECRET")
if (nzchar(t) && nzchar(s)) rsconnect::setAccountInfo(name = a, token = t, secret = s)
rsconnect::deployApp("app", appName = "pennycress_snrna", account = a,
                     appFiles = c("app.R", "data/meta.rds", "data/expr.bin"),
                     forceUpdate = TRUE, launch.browser = FALSE)'
