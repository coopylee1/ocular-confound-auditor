# =============================================================================
# fig_revision.R  --  regenerate the leakage figure for the camera-ready.
#
# The submitted figure was titled "Neural leakage is BOUNDED by the
# cognitive-vs-OCULAR subspace overlap". Reviewer 1 objected to both words: the
# score is a geometric diagnostic rather than a derived bound, and the subspace
# is gaze-predictive rather than demonstrably ocular. This regenerates it from
# the E9 sweep with corrected terminology and 95% confidence intervals.
# =============================================================================
OUT <- "camera_ready"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

E9 <- read.csv("results/tables/E9_haufe_overlap.csv")
E9 <- E9[order(E9$s_weights), ]
col_pos <- "#c0392b"

png(file.path(OUT, "fig8_leakage_angle.png"), width = 1200, height = 820, res = 150)
par(mar = c(4.8, 4.6, 3.0, 1))
plot(E9$s_weights, E9$opi_leak, type = "b", pch = 19, col = col_pos, lwd = 2,
     xlab = expression("subspace overlap  s = ||P t|| / ||t|| = cos(principal angle)"),
     ylab = "pure neural-leakage OPI (no ocular coupling)",
     main = "Neural leakage tracks the cognitive-vs-gaze-predictive subspace overlap",
     cex.main = 0.95,
     ylim = range(c(E9$opi_leak - E9$opi_ci, E9$opi_leak + E9$opi_ci)))
arrows(E9$s_weights, E9$opi_leak - E9$opi_ci,
       E9$s_weights, E9$opi_leak + E9$opi_ci,
       angle = 90, code = 3, length = 0.03, col = col_pos)
abline(h = 0, col = "gray60", lty = 3)
text(E9$s_weights[1], E9$opi_leak[1] + 0.012,
     "posterior source\n(near-orthogonal → leakage ≈ 0)",
     pos = 4, cex = 0.72, col = "#1e8449")
text(max(E9$s_weights), E9$opi_leak[which.max(E9$s_weights)],
     "frontal source\n(aligned → leakage)", pos = 2, cex = 0.72, col = "#a93226")
legend("topleft",
       sprintf("Spearman rho = %.2f  (95%% CI bars over 8 studies)",
               cor(E9$s_weights, E9$opi_leak, method = "spearman")),
       bty = "n", cex = 0.85)
dev.off()

cat(sprintf("wrote %s/fig8_leakage_angle.png\n", OUT))
cat(sprintf("  s range %.3f-%.3f   leakage %.3f-%.3f   rho=%.3f\n",
            min(E9$s_weights), max(E9$s_weights),
            min(E9$opi_leak), max(E9$opi_leak),
            cor(E9$s_weights, E9$opi_leak, method = "spearman")))
