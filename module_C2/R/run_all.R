## =====================================================================
## module_C2 / R/run_all.R -- driver. NB 02 depends on 03's cache.
## Rscript module_C2/R/run_all.R
## =====================================================================
for (s in c("01_reproduce_current_C2", "03_anchor_loci", "02_compare_denominators",
            "04_grid_sensitivity", "05_null_resolution", "06_make_figures", "07_sim_validation")) {
  cat(sprintf("\n########## %s ##########\n", s))
  source(file.path("module_C2/R", paste0(s, ".R")), local = new.env(), echo = FALSE)
}
