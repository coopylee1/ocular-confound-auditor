# =============================================================================
# figures.R  --  publication figures from the saved simulation results.
# Reads results/tables/sim_results.rds; writes PNGs to results/figures/.
# Base-R graphics only (deterministic, no theme dependencies).
# =============================================================================
source("R/simulate.R")
R <- readRDS("results/tables/sim_results.rds")
OUT <- "results/figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
col_pos <- "#c0392b"; col_neg <- "#2c3e50"; col_dense <- "#1f77b4"; col_erp <- "#d62728"

mont <- make_montage()
draw_head <- function(xlim = c(-1.25,1.25), ylim = c(-1.15,1.45), main = "") {
  plot(NA, xlim = xlim, ylim = ylim, asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)
  th <- seq(0, 2*pi, length = 200)
  lines(cos(th), sin(th)*1.15 + 0.1, lwd = 2, col = "gray40")      # head
  polygon(c(-0.13,0,0.13), c(1.22,1.42,1.22), border = "gray40", lwd = 2)  # nose
  lines(c(-1.0,-1.08,-1.0), c(0.25,0.1,-0.05), col="gray40", lwd=2)        # ears
  lines(c( 1.0, 1.08, 1.0), c(0.25,0.1,-0.05), col="gray40", lwd=2)
}

# ---- Figure 1: montage (dense auditor vs ERP-CORE-30 target) -----------------
png(file.path(OUT, "fig1_montage.png"), width = 1500, height = 780, res = 150)
par(mfrow = c(1,2), mar = c(1,1,3,1))
C <- mont$coords
for (which in c("full","erpcore")) {
  ch <- if (which=="full") mont$full else mont$erpcore
  draw_head(main = if (which=="full") "Full auditor montage (64 ch; ERP-CORE 30 in white)"
                   else "ERP-CORE target montage (30 ch)")
  peri <- setdiff(mont$full, mont$erpcore)
  cc <- C[ch,,drop=FALSE]
  is_peri <- rownames(cc) %in% peri
  points(cc[!is_peri,"x"], cc[!is_peri,"y"], pch = 21, bg = "white", cex = 2.4, col = "gray30")
  if (any(is_peri)) points(cc[is_peri,"x"], cc[is_peri,"y"], pch = 21, bg = "#f1c40f", cex = 2.4, col = "gray30")
  text(cc[,"x"], cc[,"y"], rownames(cc), cex = 0.38)
  points(mont$eye_left[1], mont$eye_left[2], pch = 13, cex = 2.2, col = "#2980b9", lwd = 2)
  points(mont$eye_right[1], mont$eye_right[2], pch = 13, cex = 2.2, col = "#2980b9", lwd = 2)
  if (which=="full") legend("bottomright", c("ERP-CORE 30","added in full 64","eye"),
         pch = c(21,21,13), pt.bg = c("white","#f1c40f",NA), col = c("gray30","gray30","#2980b9"),
         bty = "n", cex = 0.8)
}
dev.off()

# ---- Figure 2: source topographies (forward model) ---------------------------
grid_field <- function(f, n = 120) {
  g <- seq(-1.2, 1.2, length = n); G <- expand.grid(x = g, y = g)
  inside <- (G$x^2 + ((G$y-0.1)/1.15)^2) <= 1.02
  z <- rep(NA, nrow(G)); z[inside] <- f(G$x[inside], G$y[inside])
  list(x = g, y = g, z = matrix(z, n, n))
}
pr <- function(x,y,px,py,s) exp(-((x-px)^2+(y-py)^2)/(2*s^2))
el <- mont$eye_left; er <- mont$eye_right
fields <- list(
  "HEOG / corneo-retinal (horizontal gaze)" = function(x,y) pr(x,y,er[1],er[2],0.40)-pr(x,y,el[1],el[2],0.40),
  "VEOG / blink (vertical, frontopolar)"     = function(x,y) pr(x,y,el[1],el[2],0.40)+pr(x,y,er[1],er[2],0.40),
  "Saccadic spike (focal periocular)"        = function(x,y) pr(x,y,er[1],er[2],0.22)-pr(x,y,el[1],el[2],0.22),
  "N2pc (lateralised posterior)"             = function(x,y) pr(x,y,0.55,-0.72,0.28)-pr(x,y,-0.55,-0.72,0.28)
)
png(file.path(OUT, "fig2_topographies.png"), width = 1500, height = 430, res = 150)
par(mfrow = c(1,4), mar = c(1,1,3,1))
pal <- colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(64)
for (nm in names(fields)) {
  fl <- grid_field(fields[[nm]])
  zr <- max(abs(fl$z), na.rm = TRUE)
  image(fl$x, fl$y, fl$z, col = pal, zlim = c(-zr,zr), asp = 1, axes = FALSE,
        xlab = "", ylab = "", main = nm, cex.main = 0.85)
  th <- seq(0,2*pi,length=200); lines(cos(th), sin(th)*1.15+0.1, lwd=2, col="gray40")
  polygon(c(-0.13,0,0.13), c(1.22,1.42,1.22), border="gray40", lwd=2)
}
dev.off()

# ---- Figure 3: H3 sensitivity battery ----------------------------------------
E1 <- R$E1
ord <- c("ocular_corr","blink_corr","label_shuffle","covert","motor","neutral")
lab <- c(ocular_corr="saccade-coupled\nocular (POS)", blink_corr="blink-coupled\n(POS)",
         label_shuffle="label\nshuffle", covert="covert\nattention",
         motor="motor\nimagery", neutral="neutral")
ord <- ord[ord %in% E1$scenario]                       # tolerate runs without blink_corr
m  <- sapply(ord, function(s) mean(E1$opi[E1$scenario==s]))
sd <- sapply(ord, function(s) sd(E1$opi[E1$scenario==s]))
is_pos <- ord %in% c("ocular_corr","blink_corr")
margin <- max(sapply(c("covert","motor","neutral","label_shuffle"),
                     function(s) mean(E1$opi[E1$scenario==s]))) + 0.02
png(file.path(OUT, "fig3_h3_battery.png"), width = 1400, height = 820, res = 150)
par(mar = c(4.5,4.5,3,1))
bp <- barplot(m, names.arg = lab[ord], col = ifelse(is_pos, col_pos, col_neg), border = NA,
              ylim = c(min(0,min(m-sd))-0.02, 0.5), ylab = "Ocular Predictability Index (OPI)",
              main = "H3: sensitivity battery  (saccade & blink positive controls vs negatives)",
              cex.names = 0.76)
arrows(bp, m-sd, bp, m+sd, angle = 90, code = 3, length = 0.04, col = "gray20", lwd = 1.5)
abline(h = 0, col = "gray60"); abline(h = margin, lty = 2, col = "#27ae60", lwd = 2)
text(bp[1], margin+0.03, sprintf("calibrated neural-leakage margin = %.3f", margin),
     col = "#1e8449", pos = 4, cex = 0.8, offset = -0.2)
text(bp, pmax(m+sd,0)+0.018, sprintf("%.3f", m), cex = 0.8)
dev.off()

# ---- Figure 4: graded coupling sweep -----------------------------------------
E2 <- R$E2; rho <- sort(unique(E2$coupling))
om <- sapply(rho, function(r) mean(E2$opi[E2$coupling==r]))
os <- sapply(rho, function(r) sd(E2$opi[E2$coupling==r]))
pm <- sapply(rho, function(r) mean(E2$auc_proj[E2$coupling==r]))
png(file.path(OUT, "fig4_coupling.png"), width = 1200, height = 820, res = 150)
par(mar = c(4.5,4.5,3,4.5))
plot(rho, om, type = "b", pch = 19, col = col_pos, lwd = 2, ylim = c(-0.05,0.55),
     xlab = "gaze ↔ label coupling  ρ", ylab = "OPI", main = "H3: graded sensitivity to ocular coupling")
arrows(rho, om-os, rho, om+os, angle = 90, code = 3, length = 0.03, col = col_pos)
abline(h = 0, col = "gray70", lty = 3)
par(new = TRUE); plot(rho, pm, type = "b", pch = 17, col = col_dense, lwd = 1.5,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0.45,1.02))
axis(4, col = col_dense, col.axis = col_dense); mtext("projected-EEG AUC", side = 4, line = 3, col = col_dense)
legend("topleft", c("OPI (left axis)","ocular-projection AUC (right axis)"),
       col = c(col_pos,col_dense), pch = c(19,17), lwd = 2, bty = "n", cex = 0.85)
dev.off()

# ---- Figure 5: Method-5 bias envelope ----------------------------------------
agg <- function(df, mtg) {
  s <- df[df$montage==mtg,]; a <- sort(unique(s$sacc_deg))
  data.frame(amp=a, bias=sapply(a,function(v)mean(s$bias[s$sacc_deg==v])),
                     bsd =sapply(a,function(v)sd(s$bias[s$sacc_deg==v])))
}
E3 <- R$E3; D <- agg(E3,"full"); Ee <- agg(E3,"erpcore")
png(file.path(OUT, "fig5_bias_envelope.png"), width = 1250, height = 820, res = 150)
par(mar = c(4.5,4.8,3,1))
xr <- range(D$amp)
plot(NA, xlim = xr, ylim = c(-0.09,0.05), log = "x", xlab = "saccade amplitude (deg, log scale)",
     ylab = "bias  =  auditor OPI  −  ground-truth-gaze OPI",
     main = "Method 5: measured two-sided bias envelope")
polygon(c(xr,rev(xr)), c(-0.05,-0.05,0.05,0.05), col = "#eafaf1", border = NA)  # +-0.05 band
abline(h = 0, col = "gray50")
for (dd in list(list(D,col_dense,"full (64 ch)"), list(Ee,col_erp,"ERP-CORE (30 ch)"))) {
  g <- dd[[1]]; lines(g$amp, g$bias, col = dd[[2]], lwd = 2, type = "b", pch = 19)
  arrows(g$amp, g$bias-g$bsd, g$amp, g$bias+g$bsd, angle=90, code=3, length=0.03, col=dd[[2]])
}
text(xr[1], -0.045, "conservative (auditor under-states)", pos = 4, cex = 0.75, col = "#1e8449")
legend("bottomright", c("full (64 ch)","ERP-CORE (30 ch)","±0.05 acceptable band"),
       col = c(col_dense,col_erp,"#eafaf1"), lwd = c(2,2,8), pch = c(19,19,NA), bty = "n", cex = 0.85)
dev.off()

# ---- Figure 6: reliability vs amplitude with calibrated tier thresholds -------
relagg <- function(mtg){ s<-E3[E3$montage==mtg,]; a<-sort(unique(s$sacc_deg))
  data.frame(amp=a, r2=sapply(a,function(v)mean(s$gaze_r2[s$sacc_deg==v]))) }
RD <- relagg("full"); RE <- relagg("erpcore")
tc <- R$calib$t_conf; te <- R$calib$t_expl
png(file.path(OUT, "fig6_reliability_gate.png"), width = 1250, height = 820, res = 150)
par(mar = c(4.5,4.8,3,1))
plot(NA, xlim = range(RD$amp), ylim = c(0.35,0.95), log = "x",
     xlab = "saccade amplitude (deg, log scale)", ylab = "auditor reliability  (mean gaze R²)",
     main = "Method 3: reliability gate (empirically calibrated)")
rect(xleft=10^par("usr")[1], xright=10^par("usr")[2], ybottom=tc, ytop=0.95, col="#eafaf1", border=NA)
rect(xleft=10^par("usr")[1], xright=10^par("usr")[2], ybottom=te, ytop=tc, col="#fef9e7", border=NA)
rect(xleft=10^par("usr")[1], xright=10^par("usr")[2], ybottom=0.35, ytop=te, col="#fdedec", border=NA)
abline(h = c(tc,te), lty = 2, col = c("#27ae60","#e67e22"))
lines(RD$amp, RD$r2, col = col_dense, lwd = 2.5, type = "b", pch = 19)
lines(RE$amp, RE$r2, col = col_erp, lwd = 2.5, type = "b", pch = 17)
text(RD$amp[1], 0.92, "CONFIRMATORY", col = "#1e8449", pos = 4, cex = 0.8)
text(RD$amp[1], (tc+te)/2, "EXPLORATORY", col = "#b9770e", pos = 4, cex = 0.8)
text(RD$amp[1], (te+0.36)/2, "INCONCLUSIVE", col = "#a93226", pos = 4, cex = 0.8)
legend("right", c("full (64 ch)","ERP-CORE (30 ch)"), col = c(col_dense,col_erp),
       lwd = 2.5, pch = c(19,17), bty = "n", cex = 0.85)
dev.off()

# ---- Figure 7: real-data audit (eegmmidb) ------------------------------------
rf <- "results/tables/real_eegmmidb_audit.csv"
if (file.exists(rf)) {
  rd <- read.csv(rf)
  rd$short <- ifelse(grepl("left_vs_right", rd$contrast), "L vs R imagery\n(no ocular confound)",
                                                          "rest vs move\n(frontal-proxy leak)")
  png(file.path(OUT, "fig7_real_eegmmidb.png"), width = 1150, height = 800, res = 150)
  par(mar = c(4.5,4.5,3.5,1))
  M <- rbind(ceiling = rd$auc_full, ocular = rd$auc_proj)
  bp <- barplot(M, beside = TRUE, names.arg = rd$short, col = c("#7f8c8d", col_pos),
                ylim = c(0.45,0.75), xpd = FALSE, ylab = "decoding AUC (subject-wise)",
                main = "Real EEG (PhysioNet eegmmidb, 64-ch): ceiling vs ocular subspace")
  abline(h = 0.5, lty = 3, col = "gray50")
  text(bp[1,], M[1,]+0.012, sprintf("%.2f", M[1,]), cex = 0.85)
  text(bp[2,], M[2,]+0.012, sprintf("%.2f", M[2,]), cex = 0.85)
  for (i in seq_len(nrow(rd)))
    text(mean(bp[,i]), 0.47, sprintf("OPI %+.3f\np=%.3f", rd$opi[i], rd$p_value[i]), cex = 0.72)
  legend("topright", c("full 64-ch EEG (ceiling)","ocular subspace"),
         fill = c("#7f8c8d", col_pos), bty = "n", cex = 0.85)
  dev.off()
}

# ---- Figure 8: neural-leakage term vs principal-angle overlap (E6) -----------
ef <- "results/tables/E6_leakage.rds"
if (file.exists(ef)) {
  E6 <- readRDS(ef)
  png(file.path(OUT, "fig8_leakage_angle.png"), width = 1200, height = 820, res = 150)
  par(mar = c(4.6,4.6,3.2,1))
  plot(E6$overlap, E6$opi_leak, type = "b", pch = 19, col = col_pos, lwd = 2,
       xlab = expression("ocular-subspace overlap  s = ||P t|| / ||t||  = cos(principal angle)"),
       ylab = "pure neural-leakage OPI (no ocular coupling)",
       main = "Neural leakage is bounded by the cognitive-vs-ocular subspace overlap")
  arrows(E6$overlap, E6$opi_leak-E6$opi_sd, E6$overlap, E6$opi_leak+E6$opi_sd,
         angle = 90, code = 3, length = 0.03, col = col_pos)
  abline(h = 0, col = "gray60", lty = 3)
  text(E6$overlap[1], 0.012, "posterior N2pc\n(orthogonal → leakage ≈ 0)", pos = 4, cex = 0.72, col = "#1e8449")
  text(max(E6$overlap), E6$opi_leak[which.max(E6$overlap)], "frontal source\n(aligned → leakage)",
       pos = 2, cex = 0.72, col = "#a93226")
  legend("topleft", sprintf("Spearman rho = %.2f", cor(E6$overlap, E6$opi_leak, method="spearman")),
         bty = "n", cex = 0.9)
  dev.off()
}

# ---- Figure 9: electric field of eye movement vs electrode position ----------
# The field of an eye movement depends on electrode POSITION (dipolar: opposite
# sign left vs right of the eyes) and DISTANCE from the eye (decay). A blink (a
# different artifact) has a distinct, symmetric, one-signed profile.
C64 <- mont$coords[mont$full, , drop = FALSE]
el <- mont$eye_left; er <- mont$eye_right
.pe <- function(M,p,s) exp(-((M[,"x"]-p[1])^2+(M[,"y"]-p[2])^2)/(2*s^2))
sacc_ref <- 10                                          # reference saccade size (deg)
heog_uV  <- sacc_ref * 16 * (.pe(C64,er,0.40) - .pe(C64,el,0.40))   # HEOG field / electrode
blink_uV <- 80 * .pe(C64, c(0,1.15), 0.30)              # blink field / electrode
dist_eye <- pmin(sqrt((C64[,"x"]-el[1])^2+(C64[,"y"]-el[2])^2),
                 sqrt((C64[,"x"]-er[1])^2+(C64[,"y"]-er[2])^2))
eside <- ifelse(C64[,"x"] < -0.04, "left", ifelse(C64[,"x"] > 0.04, "right", "midline"))
png(file.path(OUT, "fig9_efield.png"), width = 1550, height = 660, res = 150)
par(mfrow = c(1,2), mar = c(1,1,3,1))
fl <- grid_field(function(x,y) sacc_ref*16*(pr(x,y,er[1],er[2],0.40)-pr(x,y,el[1],el[2],0.40)))
zr <- max(abs(fl$z), na.rm = TRUE)
palf <- colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(64)
image(fl$x, fl$y, fl$z, col = palf, zlim = c(-zr,zr), asp = 1, axes = FALSE, xlab="", ylab="",
      main = sprintf("Electric field of a %d-deg eye movement (uV)", sacc_ref))
th <- seq(0,2*pi,length=200); lines(cos(th), sin(th)*1.15+0.1, lwd=2, col="gray40")
polygon(c(-0.13,0,0.13), c(1.22,1.42,1.22), border="gray40", lwd=2)
points(C64[,"x"], C64[,"y"], pch=21, bg="white", cex=0.8, col="gray30")
points(rbind(el,er), pch=13, cex=2, col="#1a5276", lwd=2)
mtext(sprintf("peak ~ %.0f uV at frontopolar; sign flips L vs R (dipole)", max(abs(heog_uV))),
      side=1, cex=0.72, col="gray25", line=-0.5)
par(mar = c(4.6,4.6,3,1))
cols <- c(left="#2166ac", right="#b2182b", midline="gray55")
plot(dist_eye, heog_uV, pch=19, col=cols[eside], cex=1.1, ylim=range(c(heog_uV,blink_uV)),
     xlab="electrode distance from nearest eye", ylab="field amplitude (uV)",
     main="Field vs distance and position", cex.main=0.98)
abline(h=0, col="gray70", lty=3)
points(dist_eye, blink_uV, pch=4, col="#1e8449", cex=1.0)
legend("topright", c("eye movement, left electrodes","eye movement, right electrodes",
                     "eye blink (all electrodes)"),
       pch=c(19,19,4), col=c(cols[["left"]],cols[["right"]],"#1e8449"), bty="n", cex=0.74)
text(max(dist_eye)*0.60, max(heog_uV)*0.55,
     "same distance ->\nopposite sign\n(left vs right of eye)", cex=0.7, col="gray25")
dev.off()

cat("Figures written:\n"); print(list.files(OUT, pattern="png$"))
