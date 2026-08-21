# =============================================================================
# Closed capture-mark-recapture, fit separately to each year
# =============================================================================
# Produces Pop[t], the input to the state-space component of the IPM.
#
# Following Creel et al. 2024: the same time bins and the same logit-normal
# individual random effect on p as the CJS model, fit OUTSIDE the IPM.  The
# posterior median of N is used as a point estimate, and its uncertainty is
# NOT propagated forward -- the state-space model's own sampling variance
# (sigma.c) absorbs it.  The point of doing this rather than using raw counts
# is that detection changed over the study, so uncorrected counts would bias
# the trend (Schaub & Kery 2022, section 4.3).
#
# Model Mh with data augmentation (Royle & Dorazio; Kery & Schaub):
#   z[i] ~ Bern(omega)                     is individual i in the population?
#   logit(p[i]) <- mu + eps[i]             individual heterogeneity in capture
#   eps[i] ~ N(0, sigma^2)
#   N <- sum(z[])
#
# Input:  lion_ipm_data.RDS  (chc, one matrix per year, from the create_ch script)
# Output: closed_capture_results.RDS, Pop_estimates.csv, and diagnostic plots

rm(list = ls())

required_packages <- c("jagsUI", "dplyr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

source("R/00_paths.R")

d <- readRDS("data/lion_ipm_data.RDS")
study_years <- d$settings$STUDY_YEARS

# --- settings ---------------------------------------------------------------

# Occasion scale within each year.
#   "monthly" - one occasion per month (8 per year for a May-Dec season)
#   "bins"    - the same 2-month bins the CJS uses (4 per year)
#
# "bins" matches Creel et al. literally, but four occasions leaves Mh badly
# under-identified in these data: sigma and N trade off along a ridge, so
# sigma correlated 0.78 with the SD of N, four of sixteen years failed to
# converge, and the resulting series swung -63% then +103% between adjacent
# years.  Closure is a within-year assumption either way, so splitting the
# same months into monthly occasions buys identifiability for free.
OCCASION_SCALE <- "monthly"

# Fit all years in one model with a single shared individual-heterogeneity SD
# (year-specific N, omega and p as before).
#
# TRIED AND NOT RECOMMENDED -- kept so the comparison can be reproduced.
# Pooling looks like an improvement on every convergence metric: all 16 years
# reach Rhat <= 1.016 against a worst case of 1.322, and the mean CV of N
# halves from 0.142 to 0.079.  But it buys that precision by making N a near
# deterministic function of survey intensity.  With sigma fixed across years,
# the only way the model can explain a year in which each lion was seen fewer
# times is to posit more undetected lions, so the correlation between
# detections-per-known-lion and N/known strengthens from -0.63 to -0.97 among
# well-sampled years.  That is precision without accuracy: tighter intervals
# around estimates that track effort rather than abundance.
POOL_SIGMA <- FALSE

chc <- switch(
  OCCASION_SCALE,
  monthly = d$chc_monthly,
  bins = d$chc,
  stop("OCCASION_SCALE must be 'monthly' or 'bins'")
)
n_bins <- ncol(chc[[1]])
cat("occasion scale:", OCCASION_SCALE, "->", n_bins, "occasions per year\n")

# Configuration tag: keeps runs of different configurations side by side
# instead of overwriting, and names the fit cache.
tag <- paste0(OCCASION_SCALE, if (POOL_SIGMA) "_pooled" else "_independent")
CACHE_FILE <- sprintf("output/closed_capture_fits_%s.RDS", tag)

# The full run is 16 models and takes about half an hour, so completed fits are
# cached.  Set REFIT to TRUE after changing anything that affects the fitting
# itself (the capture histories, augmentation, or MCMC settings); changes to
# the effort restriction or the summaries below do not need a refit.
REFIT <- FALSE
# Data augmentation: M = n observed + n augmented.  Too small and the N
# posterior piles up against the ceiling; the check below catches that.
AUG_FACTOR <- 5 # augment to AUG_FACTOR x the number of known individuals
AUG_MIN <- 300

N_CHAINS <- 3
N_ADAPT <- 1000
N_BURNIN <- 2000
N_ITER <- 12000
N_THIN <- 2

set.seed(20260820)

# --- model ------------------------------------------------------------------
cat(
  file = "models/closed_capture_Mh.txt",
  "
model {
  # Priors
  omega ~ dunif(0, 1)              # inclusion probability
  mean.p ~ dunif(0, 1)             # mean detection on the probability scale
  mu <- log(mean.p / (1 - mean.p))
  sigma ~ dunif(0, 5)              # SD of the individual random effect
  tau <- pow(sigma, -2)

  for (i in 1:M){
    z[i] ~ dbern(omega)
    eps[i] ~ dnorm(0, tau)T(-16, 16)
    logit(p[i]) <- mu + eps[i]
    p.eff[i] <- z[i] * p[i]

    for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i])
      yrep[i,j] ~ dbern(p.eff[i])   # replicate data for goodness of fit
    }
    freq[i] <- sum(yaug[i,])        # times individual i was detected
    freq.rep[i] <- sum(yrep[i,])
  }

  # Capture-frequency distribution, observed and replicated.  fobs[] is fixed
  # (yaug is data); frep[] carries the posterior, so comparing them is a
  # direct check that the Mh detection model reproduces the observed pattern
  # of how often individuals were seen.
  for (k in 1:T){
    for (i in 1:M){
      ind.obs[i,k] <- equals(freq[i], k)
      ind.rep[i,k] <- equals(freq.rep[i], k)
    }
    fobs[k] <- sum(ind.obs[,k])
    frep[k] <- sum(ind.rep[,k])
  }

  N <- sum(z[])
  sigma2 <- pow(sigma, 2)
}
"
)

params <- c("N", "omega", "mean.p", "sigma", "sigma2", "fobs", "frep", "eps")

# --- fit one year -----------------------------------------------------------
fit_year <- function(yr) {
  ch <- chc[[as.character(yr)]]
  n_known <- nrow(ch)
  if (n_known == 0) return(NULL)

  nz <- max(AUG_MIN, AUG_FACTOR * n_known) - n_known
  yaug <- rbind(ch, matrix(0, nrow = nz, ncol = ncol(ch)))
  M <- nrow(yaug)

  jd <- list(yaug = yaug, M = M, T = ncol(yaug))

  inits <- function() {
    list(
      z = rep(1, M),
      sigma = runif(1, 0.5, 2),
      mean.p = runif(1, 0.2, 0.6)
    )
  }

  cat(sprintf(
    "\n=== %d : %d known, M = %d ===\n",
    yr,
    n_known,
    M
  ))

  out <- jags(
    jd,
    inits = inits,
    parameters.to.save = params,
    model.file = "models/closed_capture_Mh.txt",
    n.chains = N_CHAINS,
    n.adapt = N_ADAPT,
    n.burnin = N_BURNIN,
    n.iter = N_ITER,
    n.thin = N_THIN,
    parallel = TRUE,
    verbose = FALSE
  )

  # Keep only the summary matrix.  The full jagsUI object carries every
  # posterior draw of eps -- one node per augmented individual -- which runs to
  # gigabytes across 16 years, and nothing downstream reads the raw samples.
  list(
    year = yr,
    n_known = n_known,
    M = M,
    out = list(summary = out$summary)
  )
}

# --- fit all years at once, sharing sigma -----------------------------------
fit_pooled <- function() {
  blocks <- lapply(study_years, function(yr) {
    ch <- chc[[as.character(yr)]]
    nz <- max(AUG_MIN, AUG_FACTOR * nrow(ch)) - nrow(ch)
    rbind(unname(ch), matrix(0L, nz, ncol(ch)))
  })
  sizes <- vapply(blocks, nrow, 0L)
  end <- cumsum(sizes)
  start <- end - sizes + 1L

  yaug <- do.call(rbind, blocks)
  M <- nrow(yaug)

  jd <- list(
    yaug = yaug,
    M = M,
    T = ncol(yaug),
    nyears = length(study_years),
    yr = rep(seq_along(study_years), sizes),
    start = start,
    end = end
  )

  cat(sprintf(
    "\n=== pooled fit: %d years, M = %d, %d occasions ===\n",
    length(study_years),
    M,
    ncol(yaug)
  ))

  out <- jags(
    jd,
    inits = function() {
      list(
        z = rep(1, M),
        sigma = runif(1, 0.5, 2),
        mean.p = runif(length(study_years), 0.2, 0.6)
      )
    },
    parameters.to.save = c("N", "omega", "mean.p", "sigma", "sigma2", "frep"),
    model.file = "models/closed_capture_Mh_pooled.txt",
    n.chains = N_CHAINS,
    n.adapt = N_ADAPT,
    n.burnin = N_BURNIN,
    n.iter = N_ITER,
    n.thin = N_THIN,
    parallel = TRUE,
    verbose = FALSE
  )

  slim <- list(summary = out$summary)
  lapply(seq_along(study_years), function(t) {
    list(
      year = study_years[t],
      n_known = nrow(chc[[as.character(study_years[t])]]),
      M = sizes[t],
      out = slim,
      idx = t,
      pooled = TRUE
    )
  })
}

if (!REFIT && file.exists(CACHE_FILE)) {
  cat("loading cached fits from", CACHE_FILE, "\n")
  cat("(set REFIT <- TRUE to refit)\n")
  fits <- readRDS(CACHE_FILE)
} else {
  fits <- if (POOL_SIGMA) {
    fit_pooled()
  } else {
    Filter(Negate(is.null), lapply(study_years, fit_year))
  }
  saveRDS(fits, CACHE_FILE)
}
names(fits) <- vapply(fits, function(f) as.character(f$year), "")

# --- collect N --------------------------------------------------------------
# In the pooled fit the parameters are indexed by year, so the row names in
# the summary differ; pull the right ones either way.
Pop_table <- bind_rows(lapply(fits, function(f) {
  s <- f$out$summary
  nm <- if (isTRUE(f$pooled)) {
    c(N = paste0("N[", f$idx, "]"), p = paste0("mean.p[", f$idx, "]"))
  } else {
    c(N = "N", p = "mean.p")
  }
  data.frame(
    year = f$year,
    n_known = f$n_known,
    M = f$M,
    N_median = s[nm["N"], "50%"],
    N_mean = s[nm["N"], "mean"],
    N_sd = s[nm["N"], "sd"],
    N_lower = s[nm["N"], "2.5%"],
    N_upper = s[nm["N"], "97.5%"],
    mean_p = s[nm["p"], "mean"],
    sigma = s["sigma", "mean"],
    Rhat_N = s[nm["N"], "Rhat"],
    stringsAsFactors = FALSE
  )
}))

Pop_table <- Pop_table |>
  mutate(
    # the estimate must not be pinned against the augmentation ceiling, and
    # can never be below the number of individuals actually identified
    ceiling_hit = N_upper > 0.95 * M,
    below_known = N_median < n_known,
    converged = Rhat_N < 1.1
  )

cat("\n\n=== population estimates ===\n")
print(
  Pop_table |>
    mutate(across(where(is.numeric), ~ round(.x, 3))) |>
    as.data.frame(),
  row.names = FALSE
)

# --- checks -----------------------------------------------------------------
cat("\n=== checks ===\n")
chk <- function(lbl, bad) {
  cat(sprintf(
    "%-46s %s\n",
    lbl,
    if (any(bad)) paste("** ", paste(Pop_table$year[bad], collapse = ", ")) else "ok"
  ))
}
chk("Rhat >= 1.1 (not converged)", !Pop_table$converged)
chk("N posterior near augmentation ceiling", Pop_table$ceiling_hit)
chk("N median below known individuals", Pop_table$below_known)
chk("sigma pressing its prior bound (> 4 of 5)", Pop_table$sigma > 4)

if (any(Pop_table$ceiling_hit)) {
  cat("\nRaise AUG_FACTOR and re-run the flagged years: the augmentation\n")
  cat("ceiling is truncating the posterior, so N is underestimated.\n")
}

if (any(Pop_table$sigma > 4)) {
  cat(
    "\nsigma near the top of its dunif(0, 5) prior means the individual\n",
    "random effect is absorbing very strong heterogeneity in detection.\n",
    "With only ", n_bins, " occasions per year the likelihood for N is then\n",
    "close to flat, so N is weakly identified and the credible intervals are\n",
    "wide -- the well-known difficulty with Mh.  The medians are still the\n",
    "best point estimates available and are what the IPM consumes, but treat\n",
    "N for the flagged years as soft, and check that the state-space model's\n",
    "sigma.c is large enough to absorb the resulting error.\n",
    sep = ""
  )
}

# --- goodness of fit --------------------------------------------------------
# 1. Capture-frequency distribution: observed vs posterior of replicated.
# 2. Q-Q plot of the individual random effects for DETECTED individuals only
#    (eps for augmented all-zero rows is drawn from the prior and says nothing
#    about fit).  This is the check the paper reports.
gof <- bind_rows(lapply(fits, function(f) {
  s <- f$out$summary
  ch <- chc[[as.character(f$year)]]
  k <- seq_len(ncol(ch))
  # observed capture-frequency distribution is fixed data either way
  obs <- vapply(k, function(kk) sum(rowSums(ch) == kk), 0L)
  rep_rows <- if (isTRUE(f$pooled)) {
    paste0("frep[", f$idx, ",", k, "]")
  } else {
    paste0("frep[", k, "]")
  }
  data.frame(
    year = f$year,
    times_detected = k,
    observed = obs,
    rep_mean = s[rep_rows, "mean"],
    rep_lower = s[rep_rows, "2.5%"],
    rep_upper = s[rep_rows, "97.5%"],
    stringsAsFactors = FALSE
  )
}))

gof <- gof |>
  mutate(inside = observed >= rep_lower & observed <= rep_upper)

cat("\n=== capture-frequency fit ===\n")
cat(
  "observed counts inside the 95% predictive interval:",
  sum(gof$inside),
  "of",
  nrow(gof),
  "\n"
)
print(
  gof |>
    filter(!inside) |>
    mutate(across(where(is.numeric), ~ round(.x, 1))) |>
    as.data.frame(),
  row.names = FALSE
)

pdf(
  sprintf(
    "output/closed_capture_diagnostics_%s.pdf",
    paste0(OCCASION_SCALE, if (POOL_SIGMA) "_pooled" else "_independent")
  ),
  width = 10,
  height = 7
)

# N through time
with(
  Pop_table,
  {
    plot(
      year,
      N_median,
      type = "b",
      pch = 16,
      ylim = range(c(N_lower, N_upper, n_known)),
      xlab = "Year",
      ylab = "Population size",
      main = "Closed-capture estimates of N"
    )
    arrows(year, N_lower, year, N_upper, code = 3, angle = 90, length = 0.03)
    points(year, n_known, pch = 4, col = "red")
    legend(
      "topleft",
      c("N (median, 95% CRI)", "individuals identified"),
      pch = c(16, 4),
      col = c("black", "red"),
      bty = "n"
    )
  }
)

# capture-frequency fit, one panel per year
op <- par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))
for (f in fits) {
  g <- gof[gof$year == f$year, ]
  plot(
    g$times_detected,
    g$observed,
    pch = 16,
    ylim = range(c(g$rep_lower, g$rep_upper, g$observed)),
    xlab = "",
    ylab = "",
    main = f$year
  )
  arrows(
    g$times_detected,
    g$rep_lower,
    g$times_detected,
    g$rep_upper,
    code = 3,
    angle = 90,
    length = 0.02,
    col = "grey40"
  )
}
par(op)
mtext("Capture frequency: observed (points) vs replicated (bars)", outer = TRUE)

# Q-Q plots of the individual random effects, detected individuals only.
# The pooled fit does not monitor eps (it would be one node per stacked row),
# so this panel is only produced for independent per-year fits.
if (!POOL_SIGMA) {
  op <- par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))
  for (f in fits) {
    s <- f$out$summary
    e <- s[grep("^eps\\[", rownames(s)), "mean"]
    e <- e[seq_len(f$n_known)]
    qqnorm(e, main = f$year, xlab = "", ylab = "", pch = 16, cex = 0.5)
    qqline(e, col = "red")
  }
  par(op)
  mtext("Q-Q plots of individual random effects on logit(p)", outer = TRUE)
}

dev.off()
cat("
wrote diagnostics PDF
")

# --- save -------------------------------------------------------------------
# --- restrict Pop[t] to years of comparable survey effort --------------------
# N from an Mh model tracks how hard you looked as well as how many lions were
# there: with no effort covariate available (tblDaySheets covers only Jan 2012
# to Apr 2013), a year in which each lion was seen fewer times is explained by
# positing more undetected lions.  Across these data the correlation between
# detections per known lion and N/known is -0.73.
#
# Creel et al. handled this by analysing density only over 2018-2021, "over
# which the intensity of monitoring was constant".  The same logic here keeps
# years at or above MIN_DET_PER_LION and sets the rest to NA.  NA entries drop
# out of the state-space likelihood in JAGS and NIMBLE alike, so the excluded
# years cost nothing elsewhere -- the CJS and fecundity components still use
# every year of data.
MIN_DET_PER_LION <- 2.5

det_per_lion <- vapply(
  Pop_table$year,
  function(yr) {
    m <- chc[[as.character(yr)]]
    sum(m) / nrow(m)
  },
  0
)

Pop_table$det_per_lion <- det_per_lion
Pop_table$effort_ok <- det_per_lion >= MIN_DET_PER_LION &
  !Pop_table$ceiling_hit &
  Pop_table$converged

Pop <- Pop_table$N_median
Pop[!Pop_table$effort_ok] <- NA_real_
names(Pop) <- Pop_table$year

cat("\n--- effort restriction on Pop[t] ---\n")
cat("threshold: detections per known lion >=", MIN_DET_PER_LION, "\n")
print(
  Pop_table |>
    transmute(
      year,
      n_known,
      det_per_lion = round(det_per_lion, 2),
      N_median,
      ceiling_hit,
      converged,
      used = effort_ok
    ) |>
    as.data.frame(),
  row.names = FALSE
)
cat(
  "years retained:",
  sum(Pop_table$effort_ok),
  "of",
  nrow(Pop_table),
  "->",
  paste(Pop_table$year[Pop_table$effort_ok], collapse = ", "),
  "\n"
)

closed_capture_results <- list(
  Pop = Pop,
  Pop_table = Pop_table,
  gof = gof,
  fits = fits,
  settings = list(
    OCCASION_SCALE = OCCASION_SCALE,
    POOL_SIGMA = POOL_SIGMA,
    AUG_FACTOR = AUG_FACTOR,
    AUG_MIN = AUG_MIN,
    n_bins = n_bins,
    study_years = study_years
  )
)

# Tag the outputs so a pooled run and an independent run can sit side by side
# for comparison instead of overwriting each other.
saveRDS(
  closed_capture_results,
  sprintf("output/closed_capture_results_%s.RDS", tag)
)
write.csv(
  Pop_table,
  sprintf("output/Pop_estimates_%s.csv", tag),
  row.names = FALSE
)

cat("\nPop[t] for the IPM state-space component:\n")
print(round(Pop, 1))
cat(
  sprintf("
Saved output/*_%s.{RDS,csv,pdf}
", tag),
  "Pass Pop (and n.occ = ", length(Pop), ") into the IPM.  Creel et al. also\n",
  "divide by the area occupied to get density -- supply that separately.\n",
  sep = ""
)
