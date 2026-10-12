## final_analysis/R_figures/palette.R
##
## Shared figure conventions for every manuscript and supplementary figure
## (PK, 2026-10-12). Colour always means HOW markers are tested -- grey
## marker-wise, blue Simes, orange consensus (near-black for the
## floor-filtered marker-wise arm of the floor decomposition). The
## association engine is never a colour: EMMAX is solid lines / filled
## points, LFMM dashed lines / open points.
COL_MARKER    <- "#7a7a7a"
COL_FILTERED  <- "#252525"
COL_SIMES     <- "#2a78d6"
COL_CONSENSUS <- "#eb6834"
APPROACH_LEVELS  <- c("Marker-wise", "Floor-filtered marker-wise", "Simes", "Consensus")
APPROACH_COLOURS <- setNames(c(COL_MARKER, COL_FILTERED, COL_SIMES, COL_CONSENSUS), APPROACH_LEVELS)
ENGINE_LINETYPES <- c(EMMAX = "solid", LFMM = "22")
ENGINE_SHAPES    <- c(EMMAX = 16, LFMM = 1)

## method ids such as "emmax_simes_region" or "lfmm_snp_region"
method_approach <- function(m) data.table::fcase(grepl("simes", m), "Simes",
                                                 grepl("consensus", m), "Consensus",
                                                 default = "Marker-wise")
method_engine <- function(m) data.table::fifelse(grepl("^lfmm", m), "LFMM", "EMMAX")

## Readable cell labels. V is inversely related to selection strength and c
## to gene flow (PK, 2026-09-18); the design has 7 of the 9 (V, c) cells.
CELL_ORDER <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
SELECTION_LABEL <- c(V0.5 = "Strong selection", V1 = "Medium selection", V2 = "Weak selection")
FLOW_LABEL <- c(c1 = "High gene flow", c1.5 = "Medium gene flow", c2 = "Low gene flow")
CELL_LABELS <- setNames(paste0(SELECTION_LABEL[sub("_c.*$", "", CELL_ORDER)], "\n",
                               FLOW_LABEL[sub("^V[0-9.]+_", "", CELL_ORDER)]), CELL_ORDER)
