# =============================================================================
# Integrated population model -- assemble inputs and run (JAGS)
# =============================================================================
# Inputs:
#   data/lion_ipm_data.RDS                          capture histories, covariates
#   output/Pop_estimates_monthly_independent.csv    closed-capture N per year
# Output:
#   output/ipm_samples<TAG>.RDS, ipm_summary<TAG>.csv, ipm_diagnostics<TAG>.pdf
#
# WHY JAGS AND NOT NIMBLE
# The CJS component has nind x n.occasions = 37,760 discrete latent alive-state
# nodes.  JAGS has specialised machinery for that structure; NIMBLE gives each
# node its own binary sampler and spends essentially all of its time there.
# Three formulations were tried in NIMBLE and none was usable:
#   * latent states as written here   -- ran, but ~70x slower than JAGS
#   * marginalised with dCJS_vs       -- configureMCMC never finished
#   * marginalised with a custom CJS  -- same, because the variable-length
#     declarations dCJS forces, y[i, f[i]:n.occasions], make NIMBLE's graph
#     analysis blow up
# models/lion_ipm_nimble.R is kept for reference; models/lion_ipm_jags.txt is
# what this script fits.
#
# DEPARTURES from the JAGS code in docs/R_JAGS_IPM_CODE.docx, all marked in the
# model file:
#
#   1. Detection p is indexed by individual only.  The published code builds
#      n.occasions identical copies of the same node per lion.
#
#   2. Stratum-specific fecundity is the analytic ZIP mean
#      psi * exp(alpha + gamma[r]) rather than mean(mu.f[1:74]) over a
#      hard-coded row range.  Equal in expectation, but it does not depend on
#      the fecundity rows staying sorted by stratum, and it drops the Monte
#      Carlo noise from the latent zero-inflation indicators.
#
#   3. Cubs are produced by adult females and split between the sexes.  The
#      published code derives female cubs from adult females and male cubs from
#      adult MALES using the same rate; since C[i] counts cubs of both sexes per
#      female, that doubles total recruitment.
#
#   4. Survival is mapped onto the age classes it was estimated for.  The CJS
#      classes are two years wide, so cub survival covers ages 0->1 and 1->2,
#      subadult 2->3 and 3->4, prime 4->5 and 5->6.  The published code applies
#      cub survival once and shifts every later class one year early.
#
#   5. The Leslie projection runs through PROJECTION_BURN_IN years before the
#      first study year, with only overall scale free per stratum.  Flat priors
#      on 28 starting abundances start the model far from the stable age
#      distribution; without the burn-in the fit showed a spurious decline from
#      599 to 207 animals and pushed sigma.c to 73 against a bound of 100.
#
# Departures 3, 4 and 5 change lambda, so results are not directly comparable
# with the published Kafue estimates.  Departures 1 and 2 leave inference
# unchanged.

rm(list = ls())

required_packages <- c("jagsUI", "coda", "dplyr", "IPMbook")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

source("R/00_paths.R")

MODEL_FILE <- "models/lion_ipm_jags.txt"

# Suffix on the output filenames.  output/ is gitignored, so it does NOT switch
# with the branch -- without this the branch run would silently overwrite the
# results main's report renders from.
OUTPUT_TAG <- "_movement"

# --- settings ---------------------------------------------------------------

# Area occupied by the monitored population, km^2.  Used only for the derived
# density series.  Replace with the kernel-utilisation estimate for these
# prides before reporting any density.
AREA_KM2 <- 4000

CUB_SEX_RATIO <- 0.5 # proportion of cubs that are female
N_INIT_LOWER <- 0.1 # bounds on the free scale parameter per stratum
N_INIT_UPPER <- 100

# Years of Leslie projection run before the first study year, so the age
# structure converges to the stable distribution rather than starting at the
# arbitrary structure implied by flat priors.  See DEPARTURE 5.
PROJECTION_BURN_IN <- 25

N_CHAINS <- 3
N_ADAPT <- 1000
N_BURNIN <- 20000
N_ITER <- 60000
N_THIN <- 10

# jagsUI runs one chain per core itself -- no cluster to wire up, and failures
# come back as ordinary R errors rather than opaque cluster messages.
PARALLEL_CHAINS <- TRUE

set.seed(20260820)

# --- inputs -----------------------------------------------------------------

d <- readRDS("data/lion_ipm_data.RDS")
jd <- d$jags.data
fecd <- d$jags.data.fec
study_years <- d$settings$STUDY_YEARS
n_occ <- length(study_years)

pop_file <- "output/Pop_estimates_monthly_independent.csv"
if (!file.exists(pop_file)) {
  stop("Run R/03_closed_capture.R first -- ", pop_file, " is missing.")
}
pop_tab <- read.csv(pop_file)

Pop <- rep(NA_real_, n_occ)
Pop[match(pop_tab$year, study_years)] <- ifelse(
  pop_tab$effort_ok,
  pop_tab$N_median,
  NA_real_
)

cat("--- Pop[t] ---\n")
print(data.frame(year = study_years, Pop = Pop), row.names = FALSE)
cat(
  "years contributing to the state-space likelihood:",
  sum(!is.na(Pop)),
  "of",
  n_occ,
  "\n"
)
if (sum(!is.na(Pop)) < 4) {
  stop("Too few usable population estimates to fit the state-space component.")
}

# --- movement between strata ------------------------------------------------
# Annual transition counts from the capture histories: of lion-years where the
# stratum was observed in BOTH year t and year t+1, how many moved each way.
# The model estimates psi12 and psi21 from these, and uses them to move
# survivors between strata in the projection.
#
# Observed annual rates: 11.6% inside -> outside (133 of 1147 lion-years) and
# 31.3% outside -> inside (101 of 323).  The equilibrium composition those
# imply is 0.729 inside, against 0.72 actually observed -- which is the whole
# argument for this branch: composition is set by dispersal, not by
# differential demography.

stratum_by_year <- matrix(NA_integer_, nrow(jd$y), n_occ)
okey <- d$occasion_key
for (t in seq_len(n_occ)) {
  cols <- okey$occasion[okey$year == study_years[t]]
  seen <- which(rowSums(jd$y[, cols, drop = FALSE]) > 0)
  for (i in seen) {
    det <- cols[jd$y[i, cols] == 1]
    stratum_by_year[i, t] <- if (mean(jd$area[i, det] == 1) > 0.5) 1L else 2L
  }
}

n12 <- n1 <- n21 <- n2 <- 0L
for (t in seq_len(n_occ - 1)) {
  a <- stratum_by_year[, t]
  b <- stratum_by_year[, t + 1]
  ok <- !is.na(a) & !is.na(b)
  n1 <- n1 + sum(a[ok] == 1)
  n2 <- n2 + sum(a[ok] == 2)
  n12 <- n12 + sum(a[ok] == 1 & b[ok] == 2)
  n21 <- n21 + sum(a[ok] == 2 & b[ok] == 1)
}

cat('\n--- annual movement between strata ---\n')
cat('inside  -> outside:', n12, 'of', n1,
    sprintf('(%.1f%%)', 100 * n12 / n1), '\n')
cat('outside -> inside :', n21, 'of', n2,
    sprintf('(%.1f%%)', 100 * n21 / n2), '\n')
cat('implied equilibrium proportion inside:',
    round((n21 / n2) / (n12 / n1 + n21 / n2), 3), '\n')

# --- composition of the detected sample -------------------------------------
# Second observation process.  On main this destroyed the fit because the
# projection had no movement term; with movement it should inform the split
# without fighting the demography.  Restricted to the window of stable spatial
# coverage: over the full series the observed proportion inside falls from
# 0.87 to 0.72 (-0.045 logit/yr, p = 0.0004), which is the study expanding
# outward rather than the population moving.  Within 2016 onward there is no
# trend (+0.039 logit/yr, p = 0.16).
PROP_YEARS <- 2016:2023
PROP_MIN_PER_STRATUM <- 5

n_seen <- colSums(!is.na(stratum_by_year))
n_in <- colSums(stratum_by_year == 1, na.rm = TRUE)
n_out <- n_seen - n_in
prop_ok <- n_in >= PROP_MIN_PER_STRATUM &
  n_out >= PROP_MIN_PER_STRATUM &
  study_years %in% PROP_YEARS
n_in_obs <- ifelse(prop_ok, n_in, NA_integer_)

cat('\n--- sample composition by stratum ---\n')
print(
  data.frame(
    year = study_years,
    detected = n_seen,
    inside = n_in,
    outside = n_out,
    prop_inside = round(n_in / n_seen, 3),
    used = prop_ok
  ),
  row.names = FALSE
)
cat('years contributing composition data:', sum(prop_ok), 'of', n_occ, '\n')

n_region <- length(unique(as.vector(jd$area)))
if (n_region < 2) {
  stop(
    "Only one protection stratum present, so beta.prot and gamma[2] are not\n",
    "identifiable.  Re-run R/02_sighting_protection.R then R/01_capture_histories.R."
  )
}

# --- data bundle ------------------------------------------------------------
# idx.obs maps study year t onto its column in the projection, which starts
# PROJECTION_BURN_IN years before the study.  Passing it as data keeps the
# arithmetic out of the array indices inside the model.

jags.data <- list(
  y = jd$y,
  f = jd$f,
  nind = jd$nind,
  n.occasions = jd$n.occasions,
  nageclass = jd$nageclass,
  age = jd$age,
  sex = jd$sex,
  area = jd$area,
  surv.exp = jd$surv.exp,
  surv.exp.cub = jd$surv.exp.cub,
  C = fecd$C,
  region = fecd$region,
  n.fec = fecd$n,
  n.region = n_region,
  Pop = Pop,
  # composition of the detected sample (second observation process)
  n.in = n_in_obs,
  n.seen = n_seen,
  # annual movement counts, which identify psi12 and psi21
  n12 = n12,
  n1 = n1,
  n21 = n21,
  n2 = n2,
  n.occ = n_occ,
  n.proj = PROJECTION_BURN_IN + n_occ,
  idx.obs = PROJECTION_BURN_IN + seq_len(n_occ),
  Area = AREA_KM2,
  cub.sex.ratio = CUB_SEX_RATIO,
  N.lower = N_INIT_LOWER,
  N.upper = N_INIT_UPPER
)

# --- inits ------------------------------------------------------------------
# zInit gives NA up to and including first capture (those nodes are
# deterministic or data) and 1 afterwards.
z_init <- IPMbook::zInit(jd$y)

start_scale <- mean(Pop, na.rm = TRUE) / (7 * 2 * n_region)

inits <- function() {
  # mean.ageclass is a PER-INTERVAL logit, so a plausible annual survival has
  # to be back-transformed through the number of occasions per year.  A naive
  # init of 1.0 implies 0.29 annual survival and sends the deterministic
  # projection into an immediate crash.
  phi_interval <- 0.8^(1 / jd$surv.exp)
  psi_init <- runif(1, 0.4, 0.6)

  list(
    z = z_init,
    mean.ageclass = qlogis(phi_interval) + runif(jd$nageclass, -0.2, 0.2),
    beta.sex = rnorm(1, 0, 0.3),
    beta.prot = rnorm(1, 0, 0.3),
    mean.p = runif(1, 0.2, 0.5),
    sigma.p = runif(1, 0.5, 1.5),
    epsilon = rnorm(jd$nind, 0, 0.5),
    # fec = psi * exp(alpha), so alpha must absorb psi for initial fecundity to
    # match the observed mean cubs per female-year
    alpha = log(max(mean(fecd$C), 0.05) / psi_init) + rnorm(1, 0, 0.1),
    psi = psi_init,
    gamma = c(NA, rnorm(n_region - 1, 0, 0.3)),
    k = as.numeric(fecd$C > 0),
    sigma.c = runif(1, 5, 30),
    N1 = rep(start_scale, n_region)
  )
}

params <- c(
  "mean.ageclass", "beta.sex", "beta.prot", "phi", "surv.annual",
  "mean.p", "sigma.p",
  "alpha", "psi", "gamma", "fec",
  "sigma.c", "Ntot", "Nprot", "prop.in", "Density",
  "lambda.tot", "lambda.prot",
  "psi12", "psi21", "prop.equil"
)

# --- run --------------------------------------------------------------------

cat("\nrunning", N_CHAINS, "chains x", N_ITER, "iterations")
cat(if (PARALLEL_CHAINS) " (parallel)\n" else " (sequential)\n")

t0 <- Sys.time()
out <- jags(
  data = jags.data,
  inits = inits,
  parameters.to.save = params,
  model.file = MODEL_FILE,
  n.chains = N_CHAINS,
  n.adapt = N_ADAPT,
  n.burnin = N_BURNIN,
  n.iter = N_ITER,
  n.thin = N_THIN,
  parallel = PARALLEL_CHAINS,
  n.cores = min(N_CHAINS, parallel::detectCores()),
  verbose = TRUE
)
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# --- summarise --------------------------------------------------------------
# Column names match what the earlier NIMBLE version produced, so
# reports/ipm_report.qmd needs no changes.

s <- out$summary
ipm_summary <- data.frame(
  parameter = rownames(s),
  mean = s[, "mean"],
  sd = s[, "sd"],
  q2.5 = s[, "2.5%"],
  q50 = s[, "50%"],
  q97.5 = s[, "97.5%"],
  Rhat = s[, "Rhat"],
  n.eff = s[, "n.eff"],
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("\n=== convergence ===\n")
stochastic <- !is.na(ipm_summary$sd) & ipm_summary$sd > 0
bad <- ipm_summary[
  stochastic & !is.na(ipm_summary$Rhat) & ipm_summary$Rhat > 1.1,
]
if (nrow(bad) == 0) {
  cat("all monitored parameters Rhat <= 1.1\n")
} else {
  cat("** Rhat > 1.1 for", nrow(bad), "parameters:\n")
  print(utils::head(bad[order(-bad$Rhat), ], 15), row.names = FALSE)
}
cat(
  "minimum effective sample size (stochastic nodes):",
  round(min(ipm_summary$n.eff[stochastic], na.rm = TRUE)),
  "\n"
)

show <- function(pattern, label) {
  rows <- grep(pattern, ipm_summary$parameter)
  if (length(rows) == 0) return(invisible(NULL))
  cat("\n---", label, "---\n")
  print(
    ipm_summary[
      rows,
      c("parameter", "mean", "sd", "q2.5", "q50", "q97.5", "Rhat")
    ] |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))),
    row.names = FALSE
  )
}

show("^surv.annual", "annual apparent survival [ageclass, sex, protection]")
show("^beta", "effects on logit(phi)")
show("^fec|^psi$|^alpha$|^gamma", "fecundity")
show("^Ntot", "total population size")
show("^lambda.tot", "population growth rate")
show("^sigma", "variance components")

# --- save -------------------------------------------------------------------

saveRDS(
  list(
    samples = out$samples,
    summary = ipm_summary,
    Pop = Pop,
    study_years = study_years,
    settings = list(
      ENGINE = "JAGS",
      AREA_KM2 = AREA_KM2,
      CUB_SEX_RATIO = CUB_SEX_RATIO,
      PROJECTION_BURN_IN = PROJECTION_BURN_IN,
      PARALLEL_CHAINS = PARALLEL_CHAINS,
      N_CHAINS = N_CHAINS,
      N_ITER = N_ITER,
      N_BURNIN = N_BURNIN,
      N_THIN = N_THIN
    )
  ),
  sprintf("output/ipm_samples%s.RDS", OUTPUT_TAG)
)
write.csv(
  ipm_summary,
  sprintf("output/ipm_summary%s.csv", OUTPUT_TAG),
  row.names = FALSE
)

pdf(
  sprintf("output/ipm_diagnostics%s.pdf", OUTPUT_TAG),
  width = 10,
  height = 7
)
key <- intersect(
  c("beta.sex", "beta.prot", "psi", "alpha", "sigma.c", "sigma.p", "mean.p"),
  ipm_summary$parameter
)
plot(out$samples[, key])

nt <- ipm_summary[grep("^Ntot\\[", ipm_summary$parameter), ]
plot(
  study_years,
  nt$q50,
  type = "n",
  ylim = range(c(nt$q2.5, nt$q97.5, Pop), na.rm = TRUE),
  xlab = "Year",
  ylab = "Population size",
  main = "IPM population size (line, 95% CRI) vs closed-capture input (points)"
)
polygon(
  c(study_years, rev(study_years)),
  c(nt$q2.5, rev(nt$q97.5)),
  col = adjustcolor("grey", 0.4),
  border = NA
)
lines(study_years, nt$q50, lwd = 2)
points(study_years, Pop, pch = 16, col = "red")
dev.off()

cat(sprintf(
  "\nSaved output/ipm_{samples.RDS, summary.csv, diagnostics.pdf} tagged '%s'\n",
  OUTPUT_TAG
))
