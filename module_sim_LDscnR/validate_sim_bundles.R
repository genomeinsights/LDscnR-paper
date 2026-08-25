suppressMessages({ library(data.table) })
D <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs"
fs <- sort(list.files(D, pattern = "[.]rds$"))
need <- c("GTs","map","env","LD_decay","ld_ws","GRM","grm_markers","grm_method","emx_gif")

rows <- rbindlist(lapply(fs, function(f) {
  x <- tryCatch(readRDS(file.path(D, f)), error = function(e) NULL)
  if (is.null(x)) return(data.table(file = f, ok = FALSE, err = "UNREADABLE"))
  m <- as.data.table(x$map)
  q <- p.adjust(m$emx_p, "BH")
  data.table(file = f, ok = TRUE, err = NA_character_,
    chr = as.integer(sub(".*_chr([0-9]+)_.*","\\1",f)), env = as.integer(sub(".*_env([0-9]+)[.]rds","\\1",f)),
    missing = paste(setdiff(need, names(x)), collapse=","),
    n_ind = nrow(x$GTs), n_snp = ncol(x$GTs),
    grm_ok = !is.null(x$GRM) && identical(dim(x$GRM), c(nrow(x$GTs), nrow(x$GTs))),
    grm_n = length(x$grm_markers), grm_method = x$grm_method %||% NA_character_,
    decay_ok = !is.null(x$LD_decay$decay_sum) && nrow(x$LD_decay$decay_sum) == 2L,
    ldw_ok = identical(nrow(x$ld_ws), ncol(x$GTs)),
    map_ok = nrow(m) == ncol(x$GTs),
    n_qtn = sum(m$true_QTN %in% TRUE),
    qtn_chr2 = sum(m$true_QTN %in% TRUE & m$Chr == "Chr2"),
    gif = round(x$emx_gif, 3), n_na_p = sum(is.na(m$emx_p)),
    n_bh05 = sum(q < 0.05, na.rm = TRUE),
    envsum = digest_env <- paste0(round(sum(x$env$env), 6), "_", nrow(x$env)),
    mk1 = colnames(x$GTs)[1])
}), fill = TRUE)
`%||%` <- function(a,b) if (is.null(a)) b else a

cat("=== integrity ===\n")
cat("  unreadable:", sum(!rows$ok), "| missing fields:", sum(nzchar(rows$missing), na.rm=TRUE), "\n")
cat("  GRM ok:", sum(rows$grm_ok), "/", nrow(rows), " | decay ok:", sum(rows$decay_ok),
    " | ld_ws aligned:", sum(rows$ldw_ok), " | map aligned:", sum(rows$map_ok), "\n")
cat("  grm_method:", paste(unique(rows$grm_method), collapse=","), "| grm_markers range:",
    paste(range(rows$grm_n), collapse="-"), "\n")
cat("  n_ind:", paste(unique(rows$n_ind), collapse=","), "| n_snp range:", paste(range(rows$n_snp), collapse="-"), "\n")
cat("  NA emx_p:", sum(rows$n_na_p), "\n")

cat("\n=== pooling invariant: within an env, all 10 chr must share the SAME individuals/env ===\n")
bad <- rows[, .(distinct_envsum = uniqueN(envsum), n = .N), by = env][distinct_envsum != 1]
if (nrow(bad)) { cat("  VIOLATED for env:", paste(bad$env, collapse=","), "\n"); print(bad) } else
  cat("  OK: all 10 chr share one env vector within each of the 10 env cells\n")
cat("  markers disjoint across chr within env1:",
    { s <- rows[env==1]; uniqueN(s$mk1) == nrow(s) }, "\n")

cat("\n=== truth ===\n")
print(rows[, .(n_files=.N, total_QTN=sum(n_qtn), QTN_on_neutral_Chr2=sum(qtn_chr2),
               files_with_0_QTN=sum(n_qtn==0)), by=env][order(env)])

cat("\n=== calibration + detectability (single-SNP BH q<0.05, per file) ===\n")
print(rows[, .(median_gif=median(gif), min_gif=min(gif), max_gif=max(gif),
               files_with_signal=sum(n_bh05>0), total_hits=sum(n_bh05)), by=env][order(env)])
cat("\noverall: ", sum(rows$n_bh05>0), "of", nrow(rows), "files have single-SNP signal\n")
saveRDS(rows, "/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/f5c2953d-f16b-4266-bda5-08c843e9b161/scratchpad/newsim_validation.rds")
