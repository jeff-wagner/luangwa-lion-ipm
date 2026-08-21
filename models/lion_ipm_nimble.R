# =============================================================================
# Integrated population model for Luangwa lions -- NIMBLE model code
# =============================================================================
# Structure follows Creel et al. 2024 (Conserv Sci Pract): a CJS survival
# model, a zero-inflated Poisson fecundity model, a Leslie projection linking
# the two, and a Gaussian state-space observation model for population size,
# all in one joint likelihood.
#
# Sourced by R/04_ipm.R, which supplies the data, constants and inits.
#
# Four places where this departs from the JAGS code in docs/R_JAGS_IPM_CODE.docx
# are marked DEPARTURE below and explained in R/04_ipm.R.

library(nimble)
library(nimbleEcology) # supplies dCJS_vs, used to marginalise the latent states

# Look up a whole individual's survival vector in one node.
#
# phi has only nageclass x 2 sexes x 2 strata = 16 distinct values, and each
# lion's survival trajectory is a sequence of picks from those 16.  Written the
# obvious way -- phi.ind[i, t] <- phi[age[i, t], sex[i], area[i, t]] inside a
# loop -- that is nind x (n.occasions - 1) = 37,170 separate scalar nodes, and
# NIMBLE has to generate and compile C++ for every one of them.  That is what
# makes compilation appear to hang: it is not the sampling, and marginalising
# the latent states does not help, because the node count is unchanged.
#
# BUGS code cannot index a vector with a vector of indices, so this does it in
# a nimbleFunction.  The result is one vector node per lion -- 590 instead of
# 37,170.
# Marginalised CJS likelihood with a usable length argument.
#
# nimbleEcology's dCJS_vs insists that `len` equals the length of the data, so
# every lion needs its own variable-length declaration:
#     y[i, f[i]:n.occasions] ~ dCJS_vs(probSurvive = phi.ind[i, f[i]:63], ...)
# NIMBLE's graph analysis handles that badly.  configureMCMC did not finish in
# ten CPU-minutes on a 60-lion subset, and never finished at all on the full
# 590 -- which is what looks like a hang during compilation.
#
# This version takes histories left-aligned in a fixed-width matrix and stops
# at len[i], so every declaration in the model has identical dimensions and the
# graph is trivial.  The likelihood is the standard two-state forward
# algorithm, conditioning on first capture, renormalised at each step to avoid
# underflow over long histories.
dCJS_len <- nimbleFunction(
  run = function(
    x = double(1),
    probSurvive = double(1),
    probCapture = double(0),
    len = double(0, default = 0),
    log = integer(0, default = 0)
  ) {
    returnType(double(0))
    ll <- 0
    alphaAlive <- 1
    alphaDead <- 0

    if (len > 1) {
      for (t in 2:len) {
        # survive the interval
        aA <- alphaAlive * probSurvive[t - 1]
        aD <- alphaAlive * (1 - probSurvive[t - 1]) + alphaDead
        # observe (a dead animal is never detected)
        if (x[t] > 0) {
          aA <- aA * probCapture
          aD <- 0
        } else {
          aA <- aA * (1 - probCapture)
        }
        s <- aA + aD
        if (s <= 0) {
          if (log) return(-Inf) else return(0)
        }
        ll <- ll + log(s)
        alphaAlive <- aA / s
        alphaDead <- aD / s
      }
    }

    if (log) return(ll) else return(exp(ll))
  }
)

# NIMBLE requires a matching r function to register the distribution.  Nothing
# in this model simulates capture histories, so it is a stub.
rCJS_len <- nimbleFunction(
  run = function(
    n = integer(0),
    probSurvive = double(1),
    probCapture = double(0),
    len = double(0, default = 0)
  ) {
    returnType(double(1))
    out <- numeric(length = len)
    out[1] <- 1
    return(out)
  }
)

registerDistributions(list(
  dCJS_len = list(
    BUGSdist = "dCJS_len(probSurvive, probCapture, len)",
    types = c("value = double(1)", "probSurvive = double(1)"),
    discrete = TRUE
  )
), verbose = FALSE)


phiLookup <- nimbleFunction(
  run = function(phiFlat = double(1), idx = double(1)) {
    returnType(double(1))
    n <- length(idx)
    out <- numeric(length = n)
    for (k in 1:n) {
      out[k] <- phiFlat[idx[k]]
    }
    return(out)
  }
)

lion_ipm_code <- nimbleCode({
  # ==========================================================================
  # SURVIVAL COMPONENT
  # CJS with age-class, sex and protection effects on phi, and a logit-normal
  # individual random effect on detection.
  # ==========================================================================

  beta.sex ~ dlogis(0, 1) # male effect (protection level 1 is the reference)
  beta.prot ~ dlogis(0, 1) # effect of low protection

  for (a in 1:nageclass) {
    mean.ageclass[a] ~ dlogis(0, 1)
    logit(phi[a, 1, 1]) <- mean.ageclass[a] # female, high protection
    logit(phi[a, 2, 1]) <- mean.ageclass[a] + beta.sex # male,   high
    logit(phi[a, 1, 2]) <- mean.ageclass[a] + beta.prot # female, low
    logit(phi[a, 2, 2]) <- mean.ageclass[a] + beta.sex + beta.prot # male, low
  }

  mean.p ~ dunif(0, 1)
  mu.p <- logit(mean.p)
  sigma.p ~ dunif(0, 5)

  # DEPARTURE 1: p carries an individual effect only, so it is indexed by
  # individual alone.  The published code writes the same quantity inside a
  # loop over occasions, which creates n.occasions identical nodes per lion.
  # DEPARTURE 6: the latent alive-states are marginalised out rather than
  # sampled.  Written the textbook way -- z[i,t] ~ dbern(phi * z[i,t-1]) with
  # y[i,t] ~ dbern(p * z[i,t]) -- this model carries nind x n.occasions
  # discrete latent nodes (37,760 here), and NIMBLE gives each one its own
  # binary sampler.  JAGS copes with that; NIMBLE does not.  Fitting it that
  # way ran roughly 70x slower than the same model with the states summed out,
  # and compilation alone became the bottleneck.
  #
  # dCJS_vs from nimbleEcology computes the same likelihood analytically with
  # the forward algorithm: identical inference, no z nodes, no discrete
  # samplers.  "_vs" is survival varying over occasions (through age and area)
  # with a capture probability that is constant within an individual.
  #
  # Timing convention matches the latent-state version: element 1 of
  # probSurvive is survival from occasion f[i] to f[i]+1, and so uses the age
  # class and stratum the lion occupied at f[i].  CJS conditions on first
  # capture, so detection applies from f[i]+1 onward.
  # Flatten phi so a single integer indexes age class, sex and stratum.
  # phi.idx[i, t] is precomputed in R with the matching arithmetic.
  for (a in 1:nageclass) {
    for (s in 1:2) {
      for (r in 1:n.region) {
        phiFlat[a + (s - 1) * nageclass + (r - 1) * nageclass * 2] <-
          phi[a, s, r]
      }
    }
  }

  for (i in 1:nind) {
    epsilon[i] ~ T(dnorm(0, sd = sigma.p), -16, 16)
    logit(p[i]) <- mu.p + epsilon[i]

    # Fixed-width throughout: y and phi.idx arrive left-aligned at each lion's
    # first capture and zero-padded, and len[i] says how much of the row is
    # real.  Nothing here depends on f[i], so every declaration has identical
    # dimensions.
    phi.ind[i, 1:(maxlen - 1)] <- phiLookup(
      phiFlat[1:n.phi],
      phi.idx[i, 1:(maxlen - 1)]
    )

    y[i, 1:maxlen] ~ dCJS_len(
      probSurvive = phi.ind[i, 1:(maxlen - 1)],
      probCapture = p[i],
      len = len[i]
    )
  }

  # --- annualise per-interval survival --------------------------------------
  # surv.exp and surv.exp.cub come in as data so they follow the number of
  # occasions per year set in R/01_capture_histories.R.
  for (s in 1:2) {
    for (r in 1:n.region) {
      surv[1, s, r] <- pow(phi[1, s, r], surv.exp.cub)
      for (a in 2:nageclass) {
        surv[a, s, r] <- pow(phi[a, s, r], surv.exp)
      }
    }
  }

  # ==========================================================================
  # FECUNDITY COMPONENT
  # Zero-inflated Poisson for first-year cubs per adult female per year.
  # Zero inflation is real biology here: a female with surviving cubs from the
  # previous year does not conceive (lactational infertility).
  # ==========================================================================

  alpha ~ dnorm(0, sd = 10)
  psi ~ dunif(0, 1)

  gamma[1] <- 0 # corner constraint: protection level 1 is the intercept
  for (r in 2:n.region) {
    gamma[r] ~ T(dnorm(0, sd = 3.162), -10, 10)
  }

  for (i in 1:n.fec) {
    k[i] ~ dbern(psi)
    log(lambda[i]) <- alpha + gamma[region[i]]
    mu.f[i] <- lambda[i] * k[i] + 0.00001
    C[i] ~ dpois(mu.f[i])
  }

  # DEPARTURE 2: fecundity per protection level is the analytic ZIP mean
  # rather than the average of the realised mu.f over a hard-coded range of
  # rows (fec.prot1 <- mean(mu.f[1:74]) in the published code).  Same quantity
  # in expectation, without the Monte Carlo noise from the latent k[i] and
  # without depending on the data staying sorted by region.
  for (r in 1:n.region) {
    fec[r] <- psi * exp(alpha + gamma[r])
  }

  # ==========================================================================
  # POPULATION COMPONENT
  # Leslie projection ("growth in place": no immigration or emigration), with
  # seven one-year age classes -- ages 0, 1, 2, 3, 4, 5 and 6+ -- tracked
  # separately by sex and protection level.
  # ==========================================================================

  # DEPARTURE 5: the projection is run through a burn-in of n.burn years before
  # the first study year, and only one free parameter per stratum (overall
  # scale) is estimated rather than a free abundance for each of the 7 age
  # classes x 2 sexes.
  #
  # With free per-class starting values the model begins far from the stable
  # age distribution -- flat priors put roughly 1/7 of animals in each class,
  # against a stable structure that holds ~41% in the 6+ class -- and a
  # deterministic Leslie projection then spends many years working off that
  # transient.  Over the 9 years of the published analysis that matters less
  # than it does over 16.  Fit without the burn-in, this model produced a
  # spurious decline from 599 to 207 animals, overshooting the earliest
  # observed years by 187 and forcing sigma.c to 73 against a prior bound of
  # 100.  Running the projection forward from an arbitrary structure lets it
  # converge to the stable distribution before year 1, so what the data have
  # to determine is population size, not a starting age structure they carry
  # no information about.
  for (r in 1:n.region) {
    N1[r] ~ dunif(N.lower, N.upper)
    for (a in 1:7) {
      N.f[a, 1, r] <- N1[r]
      N.m[a, 1, r] <- N1[r]
    }

    for (j in 1:(n.proj - 1)) {
      # DEPARTURE 3: recruitment.  Cubs are produced by adult FEMALES aged 3+
      # and split between the sexes, rather than female cubs coming from adult
      # females and male cubs from adult males as in the published code.  C[i]
      # counts cubs of both sexes per female, so applying the whole rate to
      # each sex separately would double total recruitment and make male
      # recruitment a function of the number of adult males.
      N.f[1, j + 1, r] <- (N.f[4, j, r] +
        N.f[5, j, r] +
        N.f[6, j, r] +
        N.f[7, j, r]) *
        fec[r] *
        cub.sex.ratio
      N.m[1, j + 1, r] <- (N.f[4, j, r] +
        N.f[5, j, r] +
        N.f[6, j, r] +
        N.f[7, j, r]) *
        fec[r] *
        (1 - cub.sex.ratio)

      # DEPARTURE 4: survival is mapped onto the age classes it was estimated
      # for.  The CJS classes are two years wide (cub 0-1.99, subadult 2-3.99,
      # prime 4-5.99, old 6+), so cub survival applies to ages 0->1 AND 1->2,
      # subadult to 2->3 and 3->4, prime to 4->5 and 5->6.  The published code
      # applies cub survival once and shifts every later class a year early.
      N.f[2, j + 1, r] <- N.f[1, j, r] * surv[1, 1, r] # age 0 -> 1
      N.f[3, j + 1, r] <- N.f[2, j, r] * surv[1, 1, r] # age 1 -> 2
      N.f[4, j + 1, r] <- N.f[3, j, r] * surv[2, 1, r] # age 2 -> 3
      N.f[5, j + 1, r] <- N.f[4, j, r] * surv[2, 1, r] # age 3 -> 4
      N.f[6, j + 1, r] <- N.f[5, j, r] * surv[3, 1, r] # age 4 -> 5
      N.f[7, j + 1, r] <- N.f[6, j, r] * surv[3, 1, r] + # age 5 -> 6+
        N.f[7, j, r] * surv[4, 1, r] # 6+ stays 6+

      N.m[2, j + 1, r] <- N.m[1, j, r] * surv[1, 2, r]
      N.m[3, j + 1, r] <- N.m[2, j, r] * surv[1, 2, r]
      N.m[4, j + 1, r] <- N.m[3, j, r] * surv[2, 2, r]
      N.m[5, j + 1, r] <- N.m[4, j, r] * surv[2, 2, r]
      N.m[6, j + 1, r] <- N.m[5, j, r] * surv[3, 2, r]
      N.m[7, j + 1, r] <- N.m[6, j, r] * surv[3, 2, r] +
        N.m[7, j, r] * surv[4, 2, r]
    }
  }

  # --- observation model ----------------------------------------------------
  # Pop[t] holds the closed-capture estimates.  Years excluded by the effort
  # restriction in R/03_closed_capture.R arrive as NA and simply contribute no
  # likelihood term.  Uncertainty in the closed-capture estimates is not
  # propagated; sigma.c absorbs it, following Creel et al.
  sigma.c ~ dunif(0.5, 100)

  for (t in 1:n.occ) {
    Ntot[t] <- sum(N.f[1:7, n.burn + t, 1:n.region]) +
      sum(N.m[1:7, n.burn + t, 1:n.region])
    Pop[t] ~ dnorm(Ntot[t], sd = sigma.c)
  }

  # ==========================================================================
  # DERIVED QUANTITIES
  # ==========================================================================

  for (t in 1:n.occ) {
    for (r in 1:n.region) {
      Nprot[r, t] <- sum(N.f[1:7, n.burn + t, r]) + sum(N.m[1:7, n.burn + t, r])
    }
    Density[t] <- Ntot[t] / Area
  }

  for (t in 1:(n.occ - 1)) {
    for (r in 1:n.region) {
      lambda.prot[r, t] <- Nprot[r, t + 1] / Nprot[r, t]
    }
    lambda.tot[t] <- Ntot[t + 1] / Ntot[t]
  }

  # annual apparent survival by age class, sex and protection, for reporting
  for (a in 1:nageclass) {
    for (s in 1:2) {
      for (r in 1:n.region) {
        surv.annual[a, s, r] <- surv[a, s, r]
      }
    }
  }
})
