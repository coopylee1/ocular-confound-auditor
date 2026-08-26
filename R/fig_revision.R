# =============================================================================
# fig_revision.R  --  leakage figure for the camera-ready.
#
# The submitted figure was titled "Neural leakage is BOUNDED by the
# cognitive-vs-OCULAR subspace overlap". Reviewer 1 objected to both words: the
# score is a geometric diagnostic rather than a derived bound, and the subspace
# is gaze-predictive rather than demonstrably ocular. Regenerated from the E9
# sweep with corrected terminology and 95% confidence intervals.
#
# Annotations are placed in the empty regions above the left tail and below the
# right shoulder, with leader arrows, so no label overlaps the data.
# =============================================================================
OUT <- "camera_ready"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

E9 <- read.csv("results/tables/E9_haufe_overlap.csv")
E9 <- E9[order(E9$s_weights), ]
col_pos <- "#c0392b"

lo <- E9$opi_leak - E9$opi_ci
hi <- E9$opi_leak + E9$opi_ci
ylim <- c(min(lo) - 0.005, max(hi) + 0.030)   # headroom for the legend

png(file.path(OUT, "fig8_leakage_angle.png"), width = 1200, height = 820, res = 150)
par(mar = c(4.8, 4.6, 3.0, 1.2))
plot(E9$s_weights, E9$opi_leak, type = "b", pch = 19, col = col_pos, lwd = 2,
     xlab = expression("subspace overlap  s = ||P t|| / ||t|| = cos(principal angle)"),
     ylab = "pure neural-leakage OPI (no ocular coupling)",
     main = "Neural leakage tracks the cognitive-vs-gaze-predictive subspace overlap",
     cex.main = 0.95, ylim = ylim)
arrows(E9$s_weights, lo, E9$s_weights, hi,
       angle = 90, code = 3, length = 0.03, col = col_pos)
abline(h = 0, col = "gray60", lty = 3)

# --- annotation 1: posterior / near-orthogonal end -------------------------
# label sits in the empty space ABOVE the low left tail; arrow points down to it
px <- E9$s_weights[1]; py <- E9$opi_leak[1]
lx <- px + 0.035; ly <- py + 0.072
text(lx, ly, "posterior source\n(near-orthogonal, leakage ~ 0)",
     pos = 4, cex = 0.72, col = "#1e8449")
arrows(lx + 0.012, ly - 0.011, px + 0.004, py + 0.016,
       length = 0.06, col = "#1e8449", lwd = 1.3)

# --- annotation 2: frontal / aligned end -----------------------------------
# label sits in the empty space BELOW the right shoulder; arrow points up to it
i  <- which.max(E9$opi_leak)
qx <- E9$s_weights[i]; qy <- E9$opi_leak[i]
mx <- qx - 0.105; my <- qy - 0.085
text(mx, my, "frontal source\n(aligned, leakage grows)",
     pos = 4, cex = 0.72, col = "#a93226")
arrows(mx + 0.055, my + 0.013, qx - 0.004, qy - 0.017,
       length = 0.06, col = "#a93226", lwd = 1.3)

legend("topleft",
       sprintf("Spearman rho = %.2f   (bars: 95%% CI over 8 studies)",
               cor(E9$s_weights, E9$opi_leak, method = "spearman")),
       bty = "n", cex = 0.82)
dev.off()
cat(sprintf("wrote %s/fig8_leakage_angle.png\n", OUT))
