# ─────────────────────────────────────────────────────────────────────────────
#  IES 2022/23 — Household Income Explorer
#  Mid-Century Modern · Browns & Warm Tones
#  Stats SA · Fact_IES2023_Households.csv
# ─────────────────────────────────────────────────────────────────────────────

library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(scales)
library(pracma)   # provides erf(), used in the Gini coefficient calculation

# Data file locations change to correct location
HOUSEHOLDS_CSV <- "C:/Users/arabe/Documents/Research_Project/Fact_IES2023_Households.csv"
GEOGRAPHY_CSV  <- "C:/Users/arabe/Documents/Research_Project/Fact_IES2023_Geography.csv"

# ── Label look-ups ─────────────────────────────────────────────────────────────
sex_lbl     <- c("1"="Male",          "2"="Female")
pop_lbl     <- c("1"="Black African", "2"="Coloured","3"="Indian/Asian","4"="White")
status_lbl  <- c("1"="Wealthy",       "2"="Very comfortable",
                 "3"="Reasonably comfortable",
                 "4"="Just getting along","5"="Poor","6"="Very poor")
dwell_lbl   <- c("1"="Formal house/brick","2"="Traditional hut",
                 "3"="Flat/apartment",    "4"="Cluster house",
                 "5"="Town house",        "6"="Semi-detached",
                 "7"="House/flat in backyard","8"="Informal shack (backyard)",
                 "9"="Informal shack (other)","10"="Room/granny flat",
                 "11"="Caravan/tent","12"="Other")
elec_lbl    <- c("1"="Yes","2"="No")
province_lbl    <- c("1"="Western Cape","2"="Eastern Cape","3"="Northern Cape",
                     "4"="Free State","5"="KwaZulu-Natal","6"="North West",
                     "7"="Gauteng","8"="Mpumalanga","9"="Limpopo")
settlement_lbl  <- c("1"="Urban","2"="Traditional","3"="Farms")

# Financial & lifestyle indicator columns shown on the Financial & Lifestyle tab
ind_vars <- c(
  "Recreation Equipment"     = "RES_RECREATION",
  "Recreation Services"      = "RES_RECSERVICES",
  "Acquired Pets"            = "RES_ACQPETS",
  "Overnight Trips Away"     = "AWA_AWAY",
  "Timeshare/Holiday Accom." = "AWA_TSHARE",
  "Mortgage Bond"            = "FAB_MORT_BOND",
  "Credit Card Debt"         = "FAB_CRED_CARD",
  "Municipal Arrears"        = "FAB_ARREAR_MUN",
  "Worried About Food"       = "LCF_ANOMONEY"
)
ind_icons <- c(
  "Recreation Equipment"     = "🎮",
  "Recreation Services"      = "🎟️",
  "Acquired Pets"            = "🐾",
  "Overnight Trips Away"     = "🧳",
  "Timeshare/Holiday Accom." = "🏖️",
  "Mortgage Bond"            = "🏠",
  "Credit Card Debt"         = "💳",
  "Municipal Arrears"        = "⚠️",
  "Worried About Food"       = "🍽️"
)

# MCM palette
PAL <- c("#A0522D","#D4783E","#C4956A","#8B7D6B",
         "#6B4226","#3D2B1F","#EDD9C0","#D4A76A","#B8865C","#7A5C3E")
SETTLE_PAL <- c("Urban"="#A0522D","Traditional"="#D4783E","Farms"="#8B7D6B")

#estimates the most common income level for a population group by smoothing the data into a 
#continuous curve and identifying the point where that curve reaches its highest peak.
density_mode <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  d <- density(x)
  d$x[which.max(d$y)]
}

# ── Contaminated lognormal (mode-parameterized) fit + Gini coefficient ─────────
# Same model/method as RP_HI_Fit_and_Gini.R, condensed to a single fit used
# on-demand for whatever subset of households the current filters select.

# n_steps kept lower here (15 vs. the 50 used in the offline RP_HI_Fit_and_Gini.R
# script) so a click of "Calculate" stays responsive; the step sizes below match
# the 15-step multi-start "nudge path" already validated for cmLN fits in
# RP_Final_Simulation.R, so 15 steps still covers a comparable lambda/epsilon
# range as the offline script's 50 finer steps.
cmLN_eps_lo <- 0.5
cmLN_lam_hi <- 50
cmLN_inv_lam <- function(l) log((l - 1) / (cmLN_lam_hi - l))
cmLN_inv_eps <- function(e) log((e - cmLN_eps_lo) / (1 - e))

dlnL_cmLN <- function(par, x) {
  m       <- exp(par[1])
  sigma   <- exp(par[2])
  mu      <- log(m) + sigma^2
  lambda  <- 1 + (cmLN_lam_hi - 1) / (1 + exp(-par[3]))
  epsilon <- cmLN_eps_lo + (1 - cmLN_eps_lo) / (1 + exp(-par[4]))
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

## Gini coefficient for the contaminated mode-parametrized log-normal.
## Identical method to rp_gini_coeff_Balidis.R / RP_HI_Fit_and_Gini.R.
cmLN_H <- function(x, m, s2) {
  0.5 + 0.5 * erf((log(x) - log(m) - s2) / sqrt(2 * s2))
}
cmLN_F_contam <- function(x, eps, m, s2, lam) {
  eps * cmLN_H(x, m, s2) + (1 - eps) * cmLN_H(x, m, lam * s2)
}
cmLN_EX_contam <- function(eps, m, s2, lam) {
  m * (eps * exp(1.5 * s2) + (1 - eps) * exp(1.5 * lam * s2))
}

# Currency formatter, shared by the server and the animated-log helpers below.
fmt_r   <- function(x) paste0("R ", format(round(x), big.mark = ","))
fmt_num <- function(x, digits = 4) format(round(x, digits), big.mark = ",", nsmall = digits)

# ── Interactive calculation log (equations + values for the "loading screen") ──
# LaTeX source, transcribed from RP_Balidis_u21495361.pdf (Section 2 & Appendix).
# NOTE: the lambda/epsilon transforms below are the app's own BOUNDED logistic
# forms (lambda in (1,50), epsilon in (0.5,1)) rather than the simpler
# unbounded forms shown in the written report -- this matches what the code
# actually optimizes over (see the "Bounded so lambda/epsilon can't run off..."
# comment elsewhere in this codebase), so the popup stays faithful to what's
# really being computed rather than the report's original toy transform.
eq_fmLN    <- r"(f_{mLN}(x;m,\sigma^2)=\dfrac{1}{x\sigma\sqrt{2\pi}}\exp\!\left[-\dfrac{(\ln(x/m)-\sigma^2)^2}{2\sigma^2}\right])"
eq_fcmLN   <- r"(f_{cmLN}(x;\varepsilon,m,\sigma^2,\lambda)=\varepsilon\,f_{mLN}(x;m,\sigma^2)+(1-\varepsilon)\,f_{mLN}(x;m,\lambda\sigma^2))"
eq_loglik  <- r"(\ln L(\varepsilon,m,\sigma^2,\lambda)=\sum_{i=1}^{n}\ln\Big[\varepsilon\,f_{mLN}(x_i;m,\sigma^2)+(1-\varepsilon)\,f_{mLN}(x_i;m,\lambda\sigma^2)\Big])"
eq_transf  <- r"(\tilde m=\ln m,\quad\tilde\sigma=\ln\sigma,\quad\tilde\lambda=\ln\!\left(\dfrac{\lambda-1}{\lambda_{hi}-\lambda}\right),\quad\tilde\varepsilon=\ln\!\left(\dfrac{\varepsilon-\varepsilon_{lo}}{1-\varepsilon}\right))"

# One "step card" of the calculation log: a title, one or more (Display-mode)
# LaTeX equations typeset by MathJax, and optional HTML detail lines
# underneath. eq_latex may be a single string or a character vector of
# several equations -- each gets its OWN \[...\] block (and its own line),
# since a bare "\\" line break inside \[...\] is only valid within an
# align/gather environment, not plain display math.
calc_step_card <- function(step_no, title, eq_latex, detail_html = "", status = "info", size = "normal") {
  cls <- switch(status,
    best = "calc-step calc-step-best",
    fail = "calc-step calc-step-fail",
    done = "calc-step calc-step-done",
    "calc-step calc-step-info"
  )
  if (size == "live") cls <- paste(cls, "calc-step-live")
  eq_html <- paste0(
    sprintf(r"(<div class="calc-step-eq-line">\[%s\]</div>)", eq_latex),
    collapse = ""
  )
  sprintf(
    r"(<div class="%s"><div class="calc-step-head"><span class="calc-step-num">Step %d</span><span class="calc-step-title">%s</span></div><div class="calc-step-eq">%s</div>%s</div>)",
    cls, step_no, title, eq_html,
    if (nzchar(detail_html)) sprintf(r"(<div class="calc-step-detail">%s</div>)", detail_html) else ""
  )
}

# Runs the multi-start search exactly as fit_cmLN_mode() used to, but instead
# of a plain progress bar it calls on_step(title, eq_latex, detail_html,
# status) once per optim() attempt (converged or not) so the caller can
# render a live, worked-math log of every step -- including which attempt is
# the current running-maximum log-likelihood.
fit_cmLN_mode_animated <- function(x, on_step, n_steps = 15,
                                   base_lambda = 1, base_eps = 0.99,
                                   lambda_step = 0.5, eps_step = -0.032) {
  start_m     <- density_mode(x)
  start_sigma <- sd(log(x[x < median(x)]))

  on_step(
    title = "Model and objective function",
    eq_latex = c(eq_fmLN, eq_fcmLN, eq_loglik, eq_transf),
    detail_html = sprintf(
      r"(<div class="calc-line">Fitting to <b>n = %s</b> filtered households.</div><div class="calc-line">Starting values (held fixed across all %d multi-starts): \(\tilde m_0=\ln(%s)=%s\), \(\tilde\sigma_0=\ln(%s)=%s\)</div>)",
      format(length(x), big.mark = ","), n_steps,
      fmt_num(start_m, 2), fmt_num(log(start_m)), fmt_num(start_sigma, 4), fmt_num(log(start_sigma))
    ),
    status = "info"
  )

  results         <- vector("list", n_steps)
  converged_flags <- logical(n_steps)
  logliks         <- rep(NA_real_, n_steps)
  best_so_far     <- -Inf

  for (i in seq_len(n_steps)) {
    curr_L <- max(1.001, min(cmLN_lam_hi - 0.001, base_lambda + i * lambda_step))
    curr_E <- max(cmLN_eps_lo + 0.001, min(0.99, base_eps + i * eps_step))
    init_pars <- c(log(start_m), log(start_sigma),
                   cmLN_inv_lam(curr_L), cmLN_inv_eps(curr_E))

    est <- try(optim(par = init_pars, fn = dlnL_cmLN, x = x,
                     control = list(fnscale = -1, maxit = 2000)),
               silent = TRUE)

    ok <- !(inherits(est, "try-error") || !is.finite(est$value) || est$convergence != 0)

    init_detail <- sprintf(
      r"(<div class="calc-line">Initial guess: \(\lambda_0=%s,\ \varepsilon_0=%s\ \Rightarrow\ \tilde\lambda_0=%s,\ \tilde\varepsilon_0=%s\)</div>)",
      fmt_num(curr_L, 3), fmt_num(curr_E, 3), fmt_num(cmLN_inv_lam(curr_L)), fmt_num(cmLN_inv_eps(curr_E))
    )

    if (!ok) {
      on_step(
        title = sprintf("Optim call %d of %d — did not converge", i, n_steps),
        eq_latex = eq_loglik,
        detail_html = paste0(init_detail, r"(<div class="calc-line calc-fail">✗ Optimizer failed to converge — skipped.</div>)"),
        status = "fail"
      )
      next
    }

    converged_flags[i] <- TRUE
    logliks[i]          <- est$value
    results[[i]]         <- est

    par       <- est$par
    m_i       <- exp(par[1]); sigma_i <- exp(par[2])
    lambda_i  <- 1 + (cmLN_lam_hi - 1) / (1 + exp(-par[3]))
    epsilon_i <- cmLN_eps_lo + (1 - cmLN_eps_lo) / (1 + exp(-par[4]))

    is_new_best <- est$value > best_so_far
    prev_best_txt <- if (i == 1 || all(!converged_flags[seq_len(i - 1)])) "none yet" else fmt_num(max(logliks[seq_len(i - 1)], na.rm = TRUE), 2)
    if (is_new_best) best_so_far <- est$value

    result_detail <- sprintf(
      r"(<div class="calc-line">Output: \(\hat m=%s,\ \hat\sigma^2=%s,\ \hat\lambda=%s,\ \hat\varepsilon=%s\)</div><div class="calc-line">\(\ln\hat L=%s\)</div><div class="calc-line %s">%s</div>)",
      fmt_r(m_i), fmt_num(sigma_i^2), fmt_num(lambda_i), fmt_num(epsilon_i), fmt_num(est$value, 2),
      if (is_new_best) "calc-best" else "calc-notbest",
      if (is_new_best) sprintf("\U0001F3C6 New maximum log-likelihood (previous best: %s)", prev_best_txt)
      else sprintf("Below current maximum of %s", fmt_num(best_so_far, 2))
    )

    on_step(
      title = sprintf("Optim call %d of %d%s", i, n_steps, if (is_new_best) " — new best" else ""),
      eq_latex = eq_loglik,
      detail_html = paste0(init_detail, result_detail),
      status = if (is_new_best) "best" else "info"
    )
  }

  n_converged <- sum(converged_flags)
  if (n_converged == 0) {
    on_step(
      title = "No converged fit",
      eq_latex = eq_loglik,
      detail_html = r"(<div class="calc-line calc-fail">✗ None of the multi-starts converged.</div>)",
      status = "fail"
    )
    return(list(m = NA, sigma = NA, lambda = NA, epsilon = NA,
               logLik = NA, n_converged = 0, n_tries = n_steps))
  }

  best_idx <- which.max(logliks)
  best_est <- results[[best_idx]]
  par <- best_est$par
  fit <- list(
    m           = exp(par[1]),
    sigma       = exp(par[2]),
    lambda      = 1 + (cmLN_lam_hi - 1) / (1 + exp(-par[3])),
    epsilon     = cmLN_eps_lo + (1 - cmLN_eps_lo) / (1 + exp(-par[4])),
    logLik      = best_est$value,
    n_converged = n_converged,
    n_tries     = n_steps
  )

  on_step(
    title = sprintf("Best fit selected (from %d/%d converged attempts)", n_converged, n_steps),
    eq_latex = eq_loglik,
    detail_html = sprintf(
      r"(<div class="calc-line">\(\hat m=%s,\ \hat\sigma^2=%s,\ \hat\lambda=%s,\ \hat\varepsilon=%s,\ \ln\hat L=%s\)</div>)",
      fmt_r(fit$m), fmt_num(fit$sigma^2), fmt_num(fit$lambda), fmt_num(fit$epsilon), fmt_num(fit$logLik, 2)
    ),
    status = "done"
  )

  pct_typical  <- round(100 * fit$epsilon, 1)
  pct_outlying <- round(100 * (1 - fit$epsilon), 1)
  on_step(
    title = "What do these results mean?",
    eq_latex = r"(\hat m,\quad \hat\varepsilon,\quad \hat\lambda)",
    detail_html = paste0(
      sprintf(
        r"(<div class="calc-line">The modal household income is <b>%s</b> — the single most common income level among these households.</div>)",
        fmt_r(fit$m)
      ),
      sprintf(
        r"(<div class="calc-line">\(\hat\varepsilon=%s\) means about <b>%s%%</b> of households have "typical" incomes, while the remaining <b>%s%%</b> are "outlying" — more dispersed incomes that pull the distribution's right tail.</div>)",
        fmt_num(fit$epsilon), fmt_num(pct_typical, 1), fmt_num(pct_outlying, 1)
      ),
      sprintf(
        r"(<div class="calc-line">\(\hat\lambda=%s\) means those outlying incomes are about <b>%s×</b> more variable than typical ones — the higher this is, the more extreme the gap between ordinary and outlying households.</div>)",
        fmt_num(fit$lambda), fmt_num(fit$lambda, 2)
      )
    ),
    status = "info"
  )

  fit
}

# Computes the Gini coefficient exactly as gini_coeff() used to, but calls
# on_step() once per equation (2)-(4) from the report's appendix, so the
# derivation from the fitted parameters to the final G is shown worked out.
gini_coeff_animated <- function(eps, m, s2, lam, on_step, tol = 1e-12) {
  EX_ref <- m * exp(1.5 * s2)
  on_step(
    title = "Expected value of the reference component",
    eq_latex = r"(E(X;m,\sigma^2)=m\,e^{1.5\sigma^2})",
    detail_html = sprintf(
      r"(<div class="calc-line">\(E(X;m,\sigma^2)=%s\times e^{1.5\times\,%s}=%s\)</div>)",
      fmt_r(m), fmt_num(s2), fmt_r(EX_ref)
    ),
    status = "info"
  )

  EX_out <- m * exp(1.5 * lam * s2)
  mu <- cmLN_EX_contam(eps, m, s2, lam)
  on_step(
    title = "Expected value of the contaminated model",
    eq_latex = r"(E(X;\varepsilon,m,\sigma^2,\lambda)=\varepsilon\,E(X;m,\sigma^2)+(1-\varepsilon)\,E(X;m,\lambda\sigma^2))",
    detail_html = sprintf(
      r"(<div class="calc-line">\(E(X)=%s\times\,%s+%s\times\,%s=%s\)</div>)",
      fmt_num(eps, 4), fmt_r(EX_ref), fmt_num(1 - eps, 4), fmt_r(EX_out), fmt_r(mu)
    ),
    status = "info"
  )

  target <- function(x) (1 - cmLN_F_contam(x, eps, m, s2, lam)) - tol
  X_hi <- uniroot(target, lower = m, upper = m * 1e8)$root
  on_step(
    title = "Contaminated CDF and integration cutoff",
    eq_latex = c(
      r"(H(x;m,\sigma^2)=0.5+0.5\,\mathrm{erf}\!\left(\dfrac{\ln x-\ln m-\sigma^2}{\sqrt{2\sigma^2}}\right))",
      r"(F(x;\varepsilon,m,\sigma^2,\lambda)=\varepsilon\,H(x;m,\sigma^2)+(1-\varepsilon)\,H(x;m,\lambda\sigma^2))"
    ),
    detail_html = sprintf(
      r"(<div class="calc-line">Solved numerically for \(X_{hi}\) where \(1-F(X_{hi})=10^{-12}\): \(X_{hi}=%s\)</div>)",
      fmt_r(X_hi)
    ),
    status = "info"
  )

  survival_sq <- function(x) (1 - cmLN_F_contam(x, eps, m, s2, lam))^2
  integrand_t <- function(t) { xx <- exp(t); survival_sq(xx) * xx }
  res <- integrate(integrand_t, lower = log(1e-6), upper = log(X_hi),
                   rel.tol = 1e-10, subdivisions = 1000)
  on_step(
    title = "Numerical integral",
    eq_latex = r"(\int_0^{\infty}\left[1-F(x;\varepsilon,m,\sigma^2,\lambda)\right]^2\,dx)",
    detail_html = sprintf(
      r"(<div class="calc-line">Evaluated via R's <code>integrate()</code> on a log-scale substitution \(x=e^t\) for numerical stability.</div><div class="calc-line">Result \(=%s\) (estimated error \(%s\))</div>)",
      format(res$value, scientific = TRUE, digits = 6), format(res$abs.error, scientific = TRUE, digits = 3)
    ),
    status = "info"
  )

  G <- 1 - res$value / mu
  on_step(
    title = "Gini coefficient",
    eq_latex = r"(G(\varepsilon,m,\sigma^2,\lambda)=1-\dfrac{1}{E(X;\varepsilon,m,\sigma^2,\lambda)}\int_0^{\infty}\left[1-F(x;\varepsilon,m,\sigma^2,\lambda)\right]^2dx)",
    detail_html = sprintf(
      r"(<div class="calc-line">\(G=1-\dfrac{%s}{%s}=%s\)</div>)",
      format(res$value, scientific = TRUE, digits = 6), fmt_r(mu), fmt_num(G, 4)
    ),
    status = "done"
  )

  list(G = G, abs.error = res$abs.error, X_hi = X_hi, mean = mu)
}

# ── ggplot theme ───────────────────────────────────────────────────────────────
theme_mcm <- function() {
  theme_minimal(base_family = "serif") +
    theme(
      plot.background   = element_rect(fill="#F5E6D3", colour=NA),
      panel.background  = element_rect(fill="#EDD9C0", colour=NA),
      panel.grid.major  = element_line(colour="#C4956A50"),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(colour="#3D2B1F", size=10),
      axis.title        = element_text(colour="#3D2B1F", size=11, face="bold"),
      plot.title        = element_text(colour="#3D2B1F", size=13, face="bold"),
      plot.subtitle     = element_text(colour="#6B4226", size=9),
      legend.background = element_rect(fill="#F5E6D3", colour=NA),
      legend.text       = element_text(colour="#3D2B1F", size=9),
      legend.title      = element_text(colour="#3D2B1F", face="bold", size=9),
      strip.background  = element_rect(fill="#6B4226", colour=NA),
      strip.text        = element_text(colour="#FAF3EA", face="bold")
    )
}

# ── CSS ────────────────────────────────────────────────────────────────────────
mcm_css <- "
@import url('https://fonts.googleapis.com/css2?family=Playfair+Display:wght@400;700&family=Lato:wght@300;400;700&display=swap');

* { box-sizing: border-box; }

html, body {
  height: 100%;
  margin: 0;
  overflow: hidden;
}

body {
  background-color: #F5E6D3;
  font-family: 'Lato', sans-serif;
  color: #3D2B1F;
  padding: 0;
}

/* ── Layout shell: header is fixed, sidebar & main scroll independently ── */
.container-fluid {
  height: 100vh;
  padding: 0 !important;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}
.app-body-row {
  flex: 1 1 auto;
  overflow: hidden;
  margin: 0 !important;
}
.app-body-row > [class*='col-'] {
  height: 100%;
}

/* ── Header ── */
.app-header {
  flex: 0 0 auto;
  background: linear-gradient(135deg, #3D2B1F 0%, #6B4226 60%, #A0522D 100%);
  color: #FAF3EA;
  padding: 22px 32px 18px;
  border-bottom: 4px solid #D4783E;
  margin-bottom: 0;
}
.app-header h1 {
  font-family: 'Playfair Display', Georgia, serif;
  font-size: 26px;
  margin: 0 0 4px 0;
  letter-spacing: 0.5px;
}
.app-header p {
  font-size: 11px;
  margin: 0;
  color: #C4956A;
  letter-spacing: 1.5px;
  text-transform: uppercase;
}

/* ── Sidebar ── */
.sidebar-wrap {
  background-color: #EDD9C0;
  border-right: 3px solid #C4956A;
  padding: 20px 16px 40px;
  height: 100%;
  overflow-y: auto;
}
.sidebar-wrap h4 {
  font-family: 'Playfair Display', serif;
  color: #3D2B1F;
  border-bottom: 2px solid #A0522D;
  padding-bottom: 5px;
  margin: 20px 0 10px;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 1.2px;
}
.sidebar-wrap h4:first-child { margin-top: 0; }

/* ── Inputs ── */
.selectize-input, .selectize-dropdown {
  background-color: #FAF3EA !important;
  border: 1px solid #C4956A !important;
  color: #3D2B1F !important;
}
.selectize-input.focus {
  border-color: #A0522D !important;
  box-shadow: 0 0 0 2px #D4783E30 !important;
}
.form-control {
  background-color: #FAF3EA;
  border: 1px solid #C4956A;
  color: #3D2B1F;
  font-size: 13px;
}
.form-control:focus {
  border-color: #A0522D;
  box-shadow: 0 0 0 2px #D4783E30;
}
.irs--shiny .irs-bar {
  background: #2E5A3E !important;
  border-top: 1px solid #2E5A3E !important;
  border-bottom: 1px solid #2E5A3E !important;
}
.irs--shiny .irs-handle {
  background: #4C8058 !important;
  border: 2px solid #2E5A3E !important;
  box-shadow: none !important;
}
.irs--shiny .irs-handle > i:first-child { background: #FAF3EA !important; }
.irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
  background: #1F3D28 !important;
}
.irs--shiny .irs-from:before, .irs--shiny .irs-to:before, .irs--shiny .irs-single:before {
  border-top-color: #1F3D28 !important;
}
.irs--shiny .irs-min, .irs--shiny .irs-max { background: #EDD9C0 !important; color: #6B4226 !important; }
.irs--shiny .irs-line { background: #D9C4A5 !important; border-color: #C4956A !important; }

/* ── Buttons ── */
.btn-mcm {
  background-color: #A0522D;
  color: #FAF3EA;
  border: none;
  font-family: 'Lato', sans-serif;
  font-weight: 700;
  letter-spacing: 0.5px;
  padding: 8px 16px;
  border-radius: 4px;
  width: 100%;
  margin-top: 8px;
  cursor: pointer;
}
.btn-mcm:hover { background-color: #6B4226; color: #FAF3EA; }
.btn-reset {
  background-color: #8B7D6B !important;
}
.btn-reset:hover { background-color: #6B4226 !important; }

/* ── Tabs ── */
.nav-tabs { border-bottom: 3px solid #A0522D; }
.nav-tabs > li > a {
  font-family: 'Lato', sans-serif;
  font-weight: 700;
  color: #6B4226;
  letter-spacing: 0.4px;
  border: 1px solid transparent;
  border-radius: 4px 4px 0 0;
}
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  background-color: #A0522D;
  color: #FAF3EA;
  border-color: #A0522D;
}
.nav-tabs > li > a:hover {
  background-color: #EDD9C0;
  border-color: #C4956A;
  color: #3D2B1F;
}
.tab-content {
  background-color: #FAF3EA;
  border: 1px solid #C4956A;
  border-top: none;
  padding: 24px;
  border-radius: 0 0 4px 4px;
}

/* ── Stat cards ── */
.stat-card {
  background: linear-gradient(135deg, #6B4226, #A0522D);
  border-radius: 6px;
  padding: 15px 18px;
  color: #FAF3EA;
  margin-bottom: 14px;
  box-shadow: 2px 3px 10px #3D2B1F25;
}
.stat-card-dark { background: linear-gradient(135deg, #3D2B1F, #6B4226); }
.stat-card-tan  { background: linear-gradient(135deg, #C4956A, #D4783E); }
.stat-card-tan .stat-label { color: #3D2B1F80; }
.stat-card-tan .stat-value { color: #3D2B1F; }
.stat-label {
  font-size: 9px;
  text-transform: uppercase;
  letter-spacing: 1.8px;
  color: #C4956A;
  margin-bottom: 5px;
}
.stat-value {
  font-family: 'Playfair Display', serif;
  font-size: 22px;
  font-weight: 700;
  line-height: 1.2;
}

/* ── Section headings ── */
.section-head {
  font-family: 'Playfair Display', serif;
  color: #3D2B1F;
  font-size: 16px;
  border-left: 4px solid #D4783E;
  padding-left: 10px;
  margin: 0 0 16px 0;
}

/* ── Badge ── */
.n-badge {
  display: inline-block;
  background-color: #2E5A3E;
  color: #FAF3EA;
  font-size: 11px;
  font-weight: 700;
  padding: 2px 11px;
  border-radius: 12px;
  margin-left: 8px;
  vertical-align: middle;
}
.weighted-note {
  margin-left: 12px;
  font-size: 11px;
  color: #6B4226;
  vertical-align: middle;
}

/* ── Data table ── */
.dataTables_wrapper { font-family: 'Lato', sans-serif; font-size: 12px; }
table.dataTable thead th {
  background-color: #6B4226 !important;
  color: #FAF3EA !important;
  border-bottom: 2px solid #A0522D !important;
}
table.dataTable tbody tr { background-color: #FAF3EA !important; }
table.dataTable tbody tr:nth-child(even) { background-color: #F5E6D3 !important; }
table.dataTable tbody tr:hover td { background-color: #EDD9C0 !important; }
.dataTables_filter input { border: 1px solid #C4956A; background: #FAF3EA; color: #3D2B1F; }
.dataTables_length select { border: 1px solid #C4956A; background: #FAF3EA; }

/* ── Summary table ── */
table.summary-tbl {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}
table.summary-tbl th {
  background-color: #6B4226;
  color: #FAF3EA;
  padding: 8px 12px;
  text-align: left;
  letter-spacing: 0.5px;
  font-size: 11px;
  text-transform: uppercase;
}
table.summary-tbl td {
  padding: 7px 12px;
  border-bottom: 1px solid #C4956A30;
  color: #3D2B1F;
}
table.summary-tbl tr:nth-child(even) td { background-color: #F5E6D3; }
table.summary-tbl tr:hover td { background-color: #EDD9C0; }

/* ── Main content ── */
.main-content { padding: 18px 22px; height: 100%; overflow-y: auto; }
.record-bar { margin-bottom: 16px; padding: 8px 0; border-bottom: 1px solid #C4956A40; }

/* ── Splash / start screen ── */
.splash-overlay {
  position: fixed;
  inset: 0;
  z-index: 9999;
  background-color: #F5E6D3;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  transition: opacity 0.6s ease;
}
.splash-overlay.splash-hidden {
  opacity: 0;
  pointer-events: none;
}
.splash-start-btn {
  width: auto;
  padding: 14px 52px;
  font-size: 16px;
  letter-spacing: 1.5px;
  border: 2px solid #2E5A3E;
}
.splash-start-btn:hover { background-color: #2E5A3E; border-color: #2E5A3E; }

/* ── Interactive calculation log (Model Fit & Gini modal) ── */
/* This app only ever opens one modal (the calc-log one), so it's safe to
   drop Bootstrap's default modal-body padding globally -- the step cards
   supply their own padding and should run edge-to-edge inside the modal. */
.modal-body { padding: 0; }
#calc-log-container {
  background-color: #FAF3EA;
  border: 1px solid #C4956A;
  border-radius: 6px;
  padding: 14px 16px;
}
.calc-step {
  border-left: 4px solid #8B7D6B;
  background-color: #F5E6D3;
  border-radius: 4px;
  padding: 10px 14px;
  margin-bottom: 12px;
}
.calc-step-best { border-left-color: #1F3D28; background-color: #A9CBAE; }
.calc-step-fail { border-left-color: #A0522D; background-color: #F7E9E2; opacity: 0.88; }
.calc-step-done { border-left-color: #3D2B1F; background-color: #EDD9C0; }
.calc-step-head {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  gap: 10px;
  margin-bottom: 6px;
}
.calc-step-num {
  font-size: 10px;
  text-transform: uppercase;
  letter-spacing: 1.2px;
  color: #8B7D6B;
  font-weight: 700;
  white-space: nowrap;
}
.calc-step-title {
  font-family: 'Playfair Display', serif;
  font-size: 14px;
  color: #3D2B1F;
  font-weight: 700;
}
.calc-step-eq {
  font-size: 13px;
  color: #3D2B1F;
  margin: 6px 0;
}
.calc-step-eq-line { margin: 6px 0; overflow-x: auto; }
.calc-step-detail { font-size: 12px; color: #6B4226; line-height: 1.5; }
.calc-line { margin: 2px 0; }
.calc-line.calc-best    { color: #1F3D28; font-weight: 700; }
.calc-line.calc-notbest { color: #8B7D6B; }
.calc-line.calc-fail    { color: #A0522D; font-weight: 700; }
.calc-ready-banner {
  background-color: #2E5A3E;
  color: #FAF3EA;
  padding: 10px 16px;
  border-radius: 6px;
  font-family: 'Playfair Display', serif;
  font-weight: 700;
  text-align: center;
  margin-bottom: 14px;
}

.gini-explainer {
  width: 100%;
  box-sizing: border-box;
  padding: 28px 32px;
  background-color: #EDD9C0;
  border-left: 4px solid #D4783E;
  border-radius: 4px;
  font-size: 16px;
  color: #3D2B1F;
  line-height: 1.7;
}

/* ── Computing placeholder shown while the real math runs (no artificial delay) ── */
.calc-computing {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  min-height: 45vh;
  gap: 16px;
  padding: 0 20px;
  font-family: 'Playfair Display', serif;
  font-size: 17px;
  color: #3D2B1F;
  text-align: center;
}
.calc-spinner { font-size: 36px; animation: calc-spin 1.4s linear infinite; }
@keyframes calc-spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }

/* ── Live one-step-at-a-time slideshow (replayed at a strict, uniform pace,
   decoupled from how long the underlying computation actually took) ── */
/* flex-direction:column (not the row default) so the single child stretches
   to fill the full width via the default cross-axis align-items:stretch --
   a row flexbox would let a narrower card (e.g. a short optim/gini card)
   shrink-wrap instead of reaching the edges like a wide one naturally would. */
#calc-log-live { min-height: 52vh; display: flex; flex-direction: column; justify-content: center; }
.calc-step-live {
  width: 100%;
  padding: 28px 34px;
  margin-bottom: 0;
}
.calc-step-live .calc-step-num   { font-size: 13px; }
.calc-step-live .calc-step-title { font-size: 21px; }
.calc-step-live .calc-step-eq    { font-size: 19px; margin: 12px 0; }
.calc-step-live .calc-step-eq-line { margin: 12px 0; }
.calc-step-live .calc-step-detail { font-size: 15px; line-height: 1.8; }
"

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  title = "IES 2022/23 Household Income Explorer",
  tags$head(tags$style(HTML(mcm_css))),
  withMathJax(),

  # ── Splash / start screen ──
  div(id = "splash-overlay", class = "splash-overlay",
    tags$button(id = "splash-start-btn", class = "btn-mcm splash-start-btn", "Start")
  ),
  tags$script(HTML("
    (function(){
      var startBtn = document.getElementById('splash-start-btn');
      var overlay  = document.getElementById('splash-overlay');

      startBtn.addEventListener('click', function(){
        overlay.classList.add('splash-hidden');
        setTimeout(function(){ overlay.style.display = 'none'; }, 650);
      });
    })();
  ")),

  # ── Live calculation log (Model Fit & Gini modal) ──────────────────────────
  # The modal's inner card is a renderUI()-bound uiOutput() that Shiny
  # patches in place on every "Next" click; this message tells MathJax to
  # typeset just that node (not the whole page -- cheap, and avoids racing
  # a removed/replaced node elsewhere). Deferred via setTimeout + a couple
  # of retries because this message and Shiny's own DOM patch for the
  # renderUI output are queued independently and can arrive in either order.
  tags$script(HTML(r"---(
    Shiny.addCustomMessageHandler('mathjax_retypeset', function(msg) {
      var id = msg && msg.id;
      var tries = 0;
      function attempt() {
        tries++;
        if (!(window.MathJax && window.MathJax.Hub)) {
          if (tries < 20) { setTimeout(attempt, 50); }
          return;
        }
        var target = id ? document.getElementById(id) : null;
        // The target element (a static div in the modal) can exist before
        // Shiny has actually patched in the new renderUI content -- wait
        // for the raw \(...\) / \[...\] markers to actually be present,
        // not just the container, or MathJax typesets an empty/stale node.
        var hasMath = target && /\\\(|\\\[/.test(target.innerHTML);
        if (id && !hasMath && tries < 20) { setTimeout(attempt, 50); return; }
        if (target) {
          MathJax.Hub.Queue(['Typeset', MathJax.Hub, target]);
        } else {
          MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
        }
      }
      setTimeout(attempt, 0);
    });
  )---")),

  # ── App header ──
  div(class = "app-header",
    tags$h1("\U0001F3E0  IES 2022/23 — Household Income Explorer"),
    tags$p("Income and Expenditure Survey · Statistics South Africa · Microdata Explorer")
  ),

  div(class = "row app-body-row",

    # ── Sidebar ──────────────────────────────────────────────────────────────
    column(3,
      div(class = "sidebar-wrap",
        uiOutput("filter_ui")
      )
    ),

    # ── Main panel ────────────────────────────────────────────────────────────
    column(9,
      div(class = "main-content",

        uiOutput("record_bar"),

        tabsetPanel(id = "tabs",

          # ── Tab 1: Summary stats ──
          tabPanel("📊  Summary Statistics",
            br(),
            p(class = "section-head", "Key Household Income & Expenditure Metrics"),

            fluidRow(
              column(4, uiOutput("card1")),
              column(4, uiOutput("card2")),
              column(4, uiOutput("card3"))
            ),
            fluidRow(
              column(4, uiOutput("card4")),
              column(4, uiOutput("card5")),
              column(4, uiOutput("card6"))
            ),
            fluidRow(
              column(6, uiOutput("card7")),
              column(6, uiOutput("card8"))
            ),

            tags$hr(style = "border-color:#C4956A; margin: 8px 0 22px;"),

            fluidRow(
              column(6,
                p(class = "section-head", "Household Income by Decile"),
                tableOutput("tbl_decile")
              ),
              column(6,
                p(class = "section-head", "Household Income by Population Group"),
                tableOutput("tbl_pop")
              )
            )
          ),

          # ── Tab 2: Data table ──
          tabPanel("📋  Data Table",
            br(),
            DTOutput("main_dt")
          ),

          # ── Tab 3: Visualisations ──
          tabPanel("📈  Visualisations",
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Household Income Distribution"),
                plotOutput("plt_hist", height = "270px")
              ),
              column(6,
                p(class = "section-head", "Mode Household Income by Decile"),
                plotOutput("plt_decile", height = "270px")
              )
            ),
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Household Income by Population Group"),
                plotOutput("plt_pop_box", height = "290px")
              ),
              column(6,
                p(class = "section-head", "Household Income per Capita by Household Size"),
                plotOutput("plt_hsize", height = "290px")
              )
            ),
            br(),
            fluidRow(
              column(12,
                p(class = "section-head", "Mode Household Income by Self-Reported Status"),
                plotOutput("plt_status", height = "250px")
              )
            )
          ),

          # ── Tab 4: Regional Profile ──
          tabPanel("🌍  Regional Profile",
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Household Income by Province"),
                plotOutput("plt_province_box", height = "290px")
              ),
              column(6,
                p(class = "section-head", "Mode Household Income by Settlement Type"),
                plotOutput("plt_settlement_bar", height = "290px")
              )
            ),
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Province Summary"),
                tableOutput("tbl_province")
              ),
              column(6,
                p(class = "section-head", "Settlement Type Summary"),
                tableOutput("tbl_settlement")
              )
            )
          ),

          # ── Tab 5: Financial & Lifestyle ──
          tabPanel("💰  Financial & Lifestyle",
            br(),
            p(class = "section-head", "Share of Households Answering “Yes”"),
            fluidRow(
              column(12,
                plotOutput("plt_indicators", height = "320px")
              )
            ),

            tags$hr(style = "border-color:#C4956A; margin: 8px 0 22px;"),

            p(class = "section-head", "At a Glance"),
            uiOutput("ind_cards")
          ),

          # ── Tab 6: Model Fit & Gini ──
          tabPanel("\U0001F9EE  Model Fit & Gini",
            br(),
            p(class = "section-head", "Mode-Parameterized Contaminated Log-Normal Fit"),
            p(style = "color:#6B4226;",
              "Fits a two-component mode-parameterized contaminated log-normal model to the currently filtered household income data using a multi-start log-likelihood search, then computes the Gini coefficient from the fitted parameters."),
            actionButton("fit_btn", "▶ Calculate Fit & Gini", class = "btn-mcm",
                        style = "width:260px;"),
            br(), br(),

            conditionalPanel(
              condition = "!input.fit_btn || input.fit_btn == 0",
              p(style = "color:#6B4226; font-style:italic;",
                "Set your filters in the sidebar, then click Calculate to run the model fit on the current selection.")
            ),

            conditionalPanel(
              condition = "input.fit_btn > 0",
              fluidRow(
                column(6,
                  p(class = "section-head", "Household Income Distribution (Filtered Sample)"),
                  plotOutput("plt_fit_hist", height = "280px")
                ),
                column(6,
                  p(class = "section-head", "Fitted Mode-Parameterized Contaminated Log-Normal Parameters"),
                  uiOutput("fit_params_ui")
                )
              ),
              br(),
              fluidRow(
                column(6,
                  p(class = "section-head", "Gini Coefficient"),
                  uiOutput("gini_ui")
                ),
                column(6,
                  div(class = "gini-explainer",
                    "The Gini coefficient measures income inequality on a scale from 0 to 1: ",
                    tags$b("0"), " means perfect equality (everyone earns the same), while ",
                    tags$b("1"), " means perfect inequality (one household earns everything)."
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  raw  <- reactiveVal(NULL)

  # ── Helper: load & join CSVs ─────────────────────────────────────────────────
  # Uqno is an 18-digit ID — too long to round-trip through a double, so it must
  # be read as character in both files or the join will silently corrupt IDs.
  load_data <- function(house_path, geo_path) {
    if (!file.exists(house_path)) {
      showNotification(paste0("File not found:\n", house_path),
                       type = "error", duration = 7)
      return(NULL)
    }
    withProgress(message = "Loading data…", value = 0.2, {
      house <- tryCatch(
        read.csv(house_path, stringsAsFactors = FALSE,
                 colClasses = c(UQNO = "character")),
        error = function(e) {
          showNotification(paste("Read error:", e$message), type="error")
          NULL
        }
      )
      if (is.null(house)) return(NULL)
      names(house) <- toupper(names(house))
      setProgress(0.6)

      if (nchar(geo_path) > 0) {
        if (file.exists(geo_path)) {
          geo <- tryCatch(
            read.csv(geo_path, stringsAsFactors = FALSE,
                     colClasses = c(Uqno = "character")),
            error = function(e) {
              showNotification(paste("Geography read error:", e$message), type="error")
              NULL
            }
          )
          if (!is.null(geo)) {
            names(geo) <- toupper(names(geo))
            geo <- geo[, c("UQNO","PROVINCE","SETTLEMENT_TYPE")]
            house <- dplyr::left_join(house, geo, by = "UQNO")
          }
        } else {
          showNotification(paste0("Geography file not found:\n", geo_path),
                           type = "warning", duration = 6)
        }
      }

      setProgress(1)
      showNotification(
        paste0("✓  Loaded ", format(nrow(house), big.mark=","), " households"),
        type = "message", duration = 4
      )
      house
    })
  }

  # Load on startup
  observe({
    if (is.null(isolate(raw()))) raw(load_data(HOUSEHOLDS_CSV, GEOGRAPHY_CSV))
  })

  # ── Dynamic filter UI ───────────────────────────────────────────────────────
  output$filter_ui <- renderUI({
    req(raw())
    df <- raw()

    age_max  <- if ("HEAD_AGE"  %in% names(df)) max(df$HEAD_AGE,  na.rm=TRUE) else 104
    size_max <- if ("HSIZE"     %in% names(df)) min(max(df$HSIZE, na.rm=TRUE), 20) else 20

    tagList(
      tags$h4("Household Income"),
      sliderInput("f_decile",   "Household Income Decile",   1, 10,       c(1,10),       step=1, ticks=FALSE),
      selectInput("f_quintile", "Household Income Quintile",
                  c("All"="all","1 — Lowest"="1","2"="2","3"="3",
                    "4"="4","5 — Highest"="5"), selected="all"),

      tags$h4("Head of Household"),
      selectInput("f_sex", "Sex",
                  c("All"="all","Male"="1","Female"="2"), selected="all"),
      selectInput("f_pop", "Population Group",
                  c("All"="all","Black African"="1","Coloured"="2",
                    "Indian/Asian"="3","White"="4"), selected="all"),
      sliderInput("f_age", "Age", 15, age_max, c(15, age_max),
                  step=1, ticks=FALSE),

      tags$h4("Household"),
      sliderInput("f_hsize", "Size", 1, size_max, c(1, size_max),
                  step=1, ticks=FALSE),
      selectInput("f_dwell", "Dwelling Type",
                  c("All"="all","Formal house/brick"="1","Traditional hut"="2",
                    "Flat/apartment"="3","Cluster house"="4",
                    "Town house"="5","Semi-detached"="6",
                    "House/flat in backyard"="7",
                    "Informal shack (backyard)"="8",
                    "Informal shack (other)"="9",
                    "Room/granny flat"="10","Caravan/tent"="11",
                    "Other"="12"), selected="all"),
      selectInput("f_status", "Present Status",
                  c("All"="all","Wealthy"="1","Very comfortable"="2",
                    "Reasonably comfortable"="3","Just getting along"="4",
                    "Poor"="5","Very poor"="6"), selected="all"),
      selectInput("f_elec", "Electricity Access",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),

      tags$h4("Geography"),
      selectInput("f_province", "Province",
                  c("All"="all",
                    "Western Cape"="1","Eastern Cape"="2","Northern Cape"="3",
                    "Free State"="4","KwaZulu-Natal"="5","North West"="6",
                    "Gauteng"="7","Mpumalanga"="8","Limpopo"="9"), selected="all"),
      selectInput("f_settlement", "Settlement Type",
                  c("All"="all","Urban"="1","Traditional"="2","Farms"="3"),
                  selected="all"),

      tags$h4("Financial & Lifestyle"),
      selectInput("f_recreation", "Recreation Equipment",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_recservices", "Recreation Services",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_pets", "Acquired Pets",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_away", "Overnight Trips Away",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_tshare", "Timeshare/Holiday Accom.",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_mortgage", "Mortgage Bond",
                  c("All"="all","Yes"="1","Not applicable"="8"), selected="all"),
      selectInput("f_credit", "Credit Card Debt",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_arrears", "Municipal Arrears",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),
      selectInput("f_food", "Worried About Food",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),

      actionButton("reset_btn", "\u21ba Reset All Filters",
                   class = "btn-mcm btn-reset")
    )
  })

  # Reset
  observeEvent(input$reset_btn, {
    df <- raw(); req(df)
    age_max  <- if ("HEAD_AGE" %in% names(df)) max(df$HEAD_AGE,  na.rm=TRUE) else 104
    size_max <- if ("HSIZE"    %in% names(df)) min(max(df$HSIZE, na.rm=TRUE), 20) else 20
    updateSliderInput(session, "f_decile",  value=c(1,10))
    updateSelectInput(session, "f_quintile",selected="all")
    updateSelectInput(session, "f_sex",     selected="all")
    updateSelectInput(session, "f_pop",     selected="all")
    updateSliderInput(session, "f_age",     value=c(15,age_max))
    updateSliderInput(session, "f_hsize",   value=c(1,size_max))
    updateSelectInput(session, "f_dwell",   selected="all")
    updateSelectInput(session, "f_status",  selected="all")
    updateSelectInput(session, "f_elec",    selected="all")
    updateSelectInput(session, "f_province",    selected="all")
    updateSelectInput(session, "f_settlement",  selected="all")
    updateSelectInput(session, "f_recreation",  selected="all")
    updateSelectInput(session, "f_recservices", selected="all")
    updateSelectInput(session, "f_pets",        selected="all")
    updateSelectInput(session, "f_away",        selected="all")
    updateSelectInput(session, "f_tshare",      selected="all")
    updateSelectInput(session, "f_mortgage",    selected="all")
    updateSelectInput(session, "f_credit",      selected="all")
    updateSelectInput(session, "f_arrears",     selected="all")
    updateSelectInput(session, "f_food",        selected="all")
  })

  # ── Filtered data ───────────────────────────────────────────────────────────
  filt <- reactive({
    req(raw())
    df <- raw()

    # Income decile range
    if (!is.null(input$f_decile) && "INCOME_DECILE" %in% names(df))
      df <- df[df$INCOME_DECILE >= input$f_decile[1] &
               df$INCOME_DECILE <= input$f_decile[2], ]

    # Income quintile
    if (!is.null(input$f_quintile) && input$f_quintile != "all" &&
        "INCOME_QUINTILE" %in% names(df))
      df <- df[df$INCOME_QUINTILE == as.integer(input$f_quintile), ]

    # Head sex
    if (!is.null(input$f_sex) && input$f_sex != "all" && "HEAD_SEX" %in% names(df))
      df <- df[df$HEAD_SEX == as.integer(input$f_sex), ]

    # Population group
    if (!is.null(input$f_pop) && input$f_pop != "all" &&
        "HEAD_POPULATION" %in% names(df))
      df <- df[df$HEAD_POPULATION == as.integer(input$f_pop), ]

    # Head age
    if (!is.null(input$f_age) && "HEAD_AGE" %in% names(df))
      df <- df[df$HEAD_AGE >= input$f_age[1] & df$HEAD_AGE <= input$f_age[2], ]

    # Household size
    if (!is.null(input$f_hsize) && "HSIZE" %in% names(df))
      df <- df[df$HSIZE >= input$f_hsize[1] & df$HSIZE <= input$f_hsize[2], ]

    # Dwelling type
    if (!is.null(input$f_dwell) && input$f_dwell != "all" &&
        "IRD_MAIND" %in% names(df))
      df <- df[df$IRD_MAIND == as.integer(input$f_dwell), ]

    # Present status
    if (!is.null(input$f_status) && input$f_status != "all" &&
        "PRESENT_STATUS" %in% names(df))
      df <- df[df$PRESENT_STATUS == as.integer(input$f_status), ]

    # Electricity
    if (!is.null(input$f_elec) && input$f_elec != "all" &&
        "ENG_ACCESS" %in% names(df))
      df <- df[df$ENG_ACCESS == as.integer(input$f_elec), ]

    # Province
    if (!is.null(input$f_province) && input$f_province != "all" &&
        "PROVINCE" %in% names(df))
      df <- df[df$PROVINCE == as.integer(input$f_province), ]

    # Settlement type
    if (!is.null(input$f_settlement) && input$f_settlement != "all" &&
        "SETTLEMENT_TYPE" %in% names(df))
      df <- df[df$SETTLEMENT_TYPE == as.integer(input$f_settlement), ]

    # Recreation equipment
    if (!is.null(input$f_recreation) && input$f_recreation != "all" &&
        "RES_RECREATION" %in% names(df))
      df <- df[df$RES_RECREATION == as.integer(input$f_recreation), ]

    # Recreation services
    if (!is.null(input$f_recservices) && input$f_recservices != "all" &&
        "RES_RECSERVICES" %in% names(df))
      df <- df[df$RES_RECSERVICES == as.integer(input$f_recservices), ]

    # Acquired pets
    if (!is.null(input$f_pets) && input$f_pets != "all" &&
        "RES_ACQPETS" %in% names(df))
      df <- df[df$RES_ACQPETS == as.integer(input$f_pets), ]

    # Overnight trips away
    if (!is.null(input$f_away) && input$f_away != "all" &&
        "AWA_AWAY" %in% names(df))
      df <- df[df$AWA_AWAY == as.integer(input$f_away), ]

    # Timeshare/holiday accommodation
    if (!is.null(input$f_tshare) && input$f_tshare != "all" &&
        "AWA_TSHARE" %in% names(df))
      df <- df[df$AWA_TSHARE == as.integer(input$f_tshare), ]

    # Mortgage bond
    if (!is.null(input$f_mortgage) && input$f_mortgage != "all" &&
        "FAB_MORT_BOND" %in% names(df))
      df <- df[df$FAB_MORT_BOND == as.integer(input$f_mortgage), ]

    # Credit card debt
    if (!is.null(input$f_credit) && input$f_credit != "all" &&
        "FAB_CRED_CARD" %in% names(df))
      df <- df[df$FAB_CRED_CARD == as.integer(input$f_credit), ]

    # Municipal arrears
    if (!is.null(input$f_arrears) && input$f_arrears != "all" &&
        "FAB_ARREAR_MUN" %in% names(df))
      df <- df[df$FAB_ARREAR_MUN == as.integer(input$f_arrears), ]

    # Worried about food
    if (!is.null(input$f_food) && input$f_food != "all" &&
        "LCF_ANOMONEY" %in% names(df))
      df <- df[df$LCF_ANOMONEY == as.integer(input$f_food), ]

    df
  })

  # ── Record bar ──────────────────────────────────────────────────────────────
  output$record_bar <- renderUI({
    req(raw())
    df  <- filt()
    n   <- nrow(df)
    tot <- nrow(raw())
    div(class = "record-bar",
      tags$span(
        style = "font-family:'Playfair Display',serif;font-size:14px;color:#3D2B1F;",
        "Filtered Households:"
      ),
      tags$span(class = "n-badge", format(n, big.mark=",")),
      tags$span(
        class = "weighted-note",
        paste0("of ", format(tot, big.mark=","), " total")
      )
    )
  })

  # ── Helpers ─────────────────────────────────────────────────────────────────
  # fmt_r() now lives at top-level (shared with the animated calc-log helpers).

  stat_card <- function(label, value, style="default") {
    cls <- switch(style,
      dark  = "stat-card stat-card-dark",
      tan   = "stat-card stat-card-tan",
      "stat-card"
    )
    div(class = cls,
      div(class = "stat-label", label),
      div(class = "stat-value", value)
    )
  }

  # ── Model Fit & Gini ─────────────────────────────────────────────────────────
  # Fits the mode-parameterized contaminated lognormal to the CURRENT filtered
  # selection only when the button is clicked (not on every filter tweak,
  # since the multi-start optimization is too slow for that).
  #
  # ALL the real work (15 optim() multi-starts + the Gini derivation) runs
  # up front, in the background, as fast as R can do it -- on_step() just
  # RECORDS each step, no UI, no delay. Once that's done, a "Next ->" button
  # in the modal pages through the already-computed steps one at a time, so
  # there is never any waiting on a click -- the click only reveals a card
  # that's already sitting there.
  modal_title <- "Fitting the Mode-Parameterized Contaminated Log-Normal Model"
  calc_steps  <- reactiveValues(list = NULL, idx = 0, x = NULL, fit = NULL, gini = NULL)
  fit_result  <- reactiveVal(NULL)

  # The modal's INNER content only -- bound to calc_steps$idx, so clicking
  # "Next" just re-renders this one output in place (a targeted DOM patch)
  # instead of the whole modal being torn down and rebuilt via showModal().
  output$calc_step_ui <- renderUI({
    req(calc_steps$list)
    i <- calc_steps$idx
    if (i < 1 || i > length(calc_steps$list)) return(NULL)
    s <- calc_steps$list[[i]]
    HTML(calc_step_card(i, s$title, s$eq_latex, s$detail_html, s$status, size = "live"))
  })

  # MathJax needs to be told explicitly to typeset calc_step_ui's new
  # content each time it changes (scoped to just that node, cheap).
  observeEvent(calc_steps$idx, {
    req(calc_steps$list, calc_steps$idx >= 1, calc_steps$idx <= length(calc_steps$list))
    session$sendCustomMessage("mathjax_retypeset", list(id = "calc-log-live"))
  }, ignoreInit = TRUE)

  # The modal's footer -- Back (hidden on step 1) + Next/"View Results" --
  # also bound to calc_steps$idx so both buttons' presence/label stay
  # correct as the user pages back and forth, without recreating the modal.
  output$calc_footer_ui <- renderUI({
    req(calc_steps$list)
    i <- calc_steps$idx
    n <- length(calc_steps$list)
    if (i < 1 || i > n) return(NULL)
    next_btn <- actionButton("calc_next_btn", if (i == n) "View Results ✓" else "Next →",
                             class = "btn-mcm", style = "flex:1; width:auto;")
    if (i <= 1) {
      div(style = "display:flex; gap:10px;", next_btn)
    } else {
      back_btn <- actionButton("calc_back_btn", "← Back", class = "btn-mcm btn-reset",
                               style = "flex:1; width:auto;")
      div(style = "display:flex; gap:10px;", back_btn, next_btn)
    }
  })

  observeEvent(input$fit_btn, {
    df <- isolate(filt())
    x  <- df$INCOME
    x  <- x[!is.na(x) & x > 0]
    validate(need(length(x) >= 30,
                  "Not enough filtered households (need at least 30) to fit a model. Widen your filters and try again."))

    showModal(modalDialog(
      title = modal_title, size = "l", easyClose = FALSE, footer = NULL,
      div(class = "calc-computing",
          div(class = "calc-spinner", "⚙"),
          div("Running the optimization…"))
    ))

    steps <- list()
    record <- function(title, eq_latex, detail_html = "", status = "info") {
      steps[[length(steps) + 1]] <<- list(title = title, eq_latex = eq_latex,
                                          detail_html = detail_html, status = status)
    }

    # A generic glossary of what each parameter means, shown before any of
    # the real computation results, so the equations that follow make sense.
    record(
      title = "What do the parameters mean?",
      eq_latex = r"(m,\quad \sigma^2,\quad \lambda,\quad \varepsilon)",
      detail_html = paste0(
        r"(<div class="calc-line"><b>\(m\)</b> — the <i>mode</i>: the single most common household income value. This is the "typical" income the whole model is anchored around.</div>)",
        r"(<div class="calc-line"><b>\(\sigma^2\)</b> — the <i>variance</i> of the typical-income group: how spread out "ordinary" incomes are around the mode.</div>)",
        r"(<div class="calc-line"><b>\(\varepsilon\)</b> — the <i>proportion typical</i>: the probability a household's income comes from the ordinary (reference) group rather than the outlying group. So \(1-\varepsilon\) is the share of outlying incomes.</div>)",
        r"(<div class="calc-line"><b>\(\lambda\)</b> — the <i>contamination inflation factor</i>, \(\lambda>1\): how much MORE variable the outlying incomes are than typical ones. A bigger \(\lambda\) means the outliers are more extreme.</div>)"
      ),
      status = "info"
    )

    fit <- fit_cmLN_mode_animated(x, on_step = record, n_steps = 15)
    g <- if (is.na(fit$logLik)) {
      NULL
    } else {
      tryCatch(gini_coeff_animated(fit$epsilon, fit$m, fit$sigma^2, fit$lambda, on_step = record),
               error = function(e) NULL)
    }

    calc_steps$list <- steps
    calc_steps$idx  <- 1
    calc_steps$x    <- x
    calc_steps$fit  <- fit
    calc_steps$gini <- g
    # fit_result is intentionally left at its previous value (if any) until
    # the user finishes paging through and clicks "View Results" -- so a
    # re-run doesn't blank out the still-valid results panel underneath
    # while the new one is being fitted/paged through.

    showModal(modalDialog(
      title = modal_title, size = "l", easyClose = FALSE,
      footer = uiOutput("calc_footer_ui"),
      div(id = "calc-log-live", uiOutput("calc_step_ui"))
    ))
  })

  observeEvent(input$calc_back_btn, {
    req(calc_steps$list)
    if (calc_steps$idx > 1) {
      calc_steps$idx <- calc_steps$idx - 1
    }
  })

  observeEvent(input$calc_next_btn, {
    # Guards against a stray/duplicate click (e.g. a double-click, or one
    # landing on the modal backdrop right as it swaps to the "results ready"
    # state) re-running the finalize branch or operating on stale state.
    req(calc_steps$list)
    if (calc_steps$idx > length(calc_steps$list)) return(invisible())

    if (calc_steps$idx < length(calc_steps$list)) {
      calc_steps$idx <- calc_steps$idx + 1
    } else {
      calc_steps$idx <- calc_steps$idx + 1   # marks this run as finalized

      log_acc <- vapply(seq_along(calc_steps$list), function(i) {
        s <- calc_steps$list[[i]]
        calc_step_card(i, s$title, s$eq_latex, s$detail_html, s$status, size = "normal")
      }, character(1))

      showModal(modalDialog(
        title = modal_title, size = "l", easyClose = TRUE, footer = modalButton("Close"),
        div(class = "calc-ready-banner",
            "✓ Results ready — scroll through the log below, or close this window to view the summary."),
        div(id = "calc-log-container", style = "max-height:65vh; overflow-y:auto; padding-right:8px;",
            HTML(paste(log_acc, collapse = "")))
      ))
      session$sendCustomMessage("mathjax_retypeset", list())

      fit_result(list(fit = calc_steps$fit, gini = calc_steps$gini, n = length(calc_steps$x), x = calc_steps$x))
    }
  })

  # ── Plot: income histogram for the data used in the fit ─────────────────────
  output$plt_fit_hist <- renderPlot({
    res <- fit_result()
    req(res, res$x)
    x   <- res$x
    cap <- quantile(x, 0.95, na.rm = TRUE)
    df2 <- data.frame(INCOME = x[x <= cap])

    ggplot(df2, aes(x = INCOME)) +
      geom_histogram(fill = "#A0522D", colour = "#3D2B1F", bins = 45, alpha = 0.9) +
      scale_x_continuous(labels = label_dollar(prefix = "R ", big.mark = ",",
                                                scale = 1e-3, suffix = "k")) +
      scale_y_continuous(labels = comma) +
      labs(x = "Annual Household Income", y = "Households",
           subtitle = paste0("n = ", format(res$n, big.mark = ","), " · capped at 95th percentile")) +
      theme_mcm()
  }, bg = "#F5E6D3")

  # ── Fitted cmLN parameter cards ──────────────────────────────────────────────
  output$fit_params_ui <- renderUI({
    res <- fit_result()
    req(res, res$fit)
    fit <- res$fit

    if (is.na(fit$logLik)) {
      return(stat_card("Fit Status", "Failed to converge — try widening your filters", "dark"))
    }

    tagList(
      fluidRow(
        column(6, stat_card("Mode (m)", fmt_r(fit$m))),
        column(6, stat_card("Sigma²", round(fit$sigma^2, 4), "dark"))
      ),
      fluidRow(
        column(6, stat_card("Lambda", round(fit$lambda, 4), "tan")),
        column(6, stat_card("Epsilon", round(fit$epsilon, 4)))
      ),
      fluidRow(
        column(6, stat_card("Log-Likelihood", format(round(fit$logLik, 2), big.mark = ","), "dark")),
        column(6, stat_card("Convergence Rate",
                            paste0(round(100 * fit$n_converged / fit$n_tries, 1), "%"), "tan"))
      )
    )
  })

  # ── Gini coefficient card ────────────────────────────────────────────────────
  output$gini_ui <- renderUI({
    res <- fit_result()
    req(res)
    g   <- res$gini

    if (is.null(g)) {
      return(stat_card("Gini Coefficient", "N/A — model fit failed", "dark"))
    }
    div(style = "max-width:320px;",
      stat_card("Gini Coefficient", round(g$G, 4), "dark")
    )
  })

  # ── Stat cards ───────────────────────────────────────────────────────────────
  output$card1 <- renderUI({
    req(filt())
    stat_card("Mode Annual Household Income",
              fmt_r(density_mode(filt()$INCOME, na.rm=TRUE)))
  })
  output$card2 <- renderUI({
    req(filt())
    stat_card("Mode Annual Expenditure",
              fmt_r(density_mode(filt()$EXPENDITURE, na.rm=TRUE)), "dark")
  })
  output$card3 <- renderUI({
    req(filt())
    stat_card("Households (sample)",
              format(nrow(filt()), big.mark=","), "tan")
  })
  output$card4 <- renderUI({
    req(filt())
    stat_card("Mean Annual Household Income",
              fmt_r(mean(filt()$INCOME, na.rm=TRUE)))
  })
  output$card5 <- renderUI({
    req(filt())
    stat_card("Mean Annual Expenditure",
              fmt_r(mean(filt()$EXPENDITURE, na.rm=TRUE)), "dark")
  })
  output$card6 <- renderUI({
    req(filt())
    stat_card("Mode Household Income Per Capita",
              fmt_r(density_mode(filt()$INCOME_PCP, na.rm=TRUE)), "tan")
  })
  output$card7 <- renderUI({
    req(filt())
    pct <- 100 * mean(filt()$FAB_MORT_BOND == 1, na.rm=TRUE)
    stat_card("Households with Mortgage Bond",
              paste0(round(pct, 1), "%"), "dark")
  })
  output$card8 <- renderUI({
    req(filt())
    pct <- 100 * mean(filt()$LCF_ANOMONEY == 1, na.rm=TRUE)
    stat_card("Worried About Food Security",
              paste0(round(pct, 1), "%"), "tan")
  })

  # ── Summary: Decile table ────────────────────────────────────────────────────
  output$tbl_decile <- renderTable({
    req(filt())
    filt() %>%
      group_by(Decile = INCOME_DECILE) %>%
      summarise(
        `Households` = n(),
        `Mode Household Income (R)` = round(density_mode(INCOME, na.rm=TRUE)),
        `Mean Household Income (R)`   = round(mean(INCOME,   na.rm=TRUE)),
        .groups = "drop"
      ) %>%
      mutate(across(c(`Mode Household Income (R)`,`Mean Household Income (R)`),
                    ~format(.x, big.mark=",")))
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="cccc",
  width="100%")

  # ── Summary: Population group table ─────────────────────────────────────────
  output$tbl_pop <- renderTable({
    req(filt())
    filt() %>%
      mutate(
        Group = dplyr::recode(as.character(HEAD_POPULATION),
                              !!!pop_lbl, .default="Other")
      ) %>%
      group_by(`Pop Group` = Group) %>%
      summarise(
        `N`                 = n(),
        `Mode Household Income (R)` = format(round(density_mode(INCOME,       na.rm=TRUE)), big.mark=","),
        `Mode Exp. (R)`   = format(round(density_mode(EXPENDITURE,  na.rm=TRUE)), big.mark=","),
        `Mode Inc/Cap (R)` = format(round(density_mode(INCOME_PCP,  na.rm=TRUE)), big.mark=","),
        .groups = "drop"
      )
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="lcccc",
  width="100%")

  # ── Data table ───────────────────────────────────────────────────────────────
  output$main_dt <- renderDT({
    req(filt())

    display_cols <- c("UQNO","HEAD_SEX","HEAD_AGE","HEAD_POPULATION","HSIZE",
                      "PRESENT_STATUS","ENG_ACCESS","INCOME","EXPENDITURE",
                      "INCOME_PCP","EXPENDITURE_PCP",
                      "INCOME_DECILE","INCOME_QUINTILE","HHOLD_WGT")

    avail <- display_cols[display_cols %in% names(filt())]
    df    <- filt()[, avail, drop=FALSE]

    # Decode factors
    if ("HEAD_SEX"        %in% names(df))
      df$HEAD_SEX        <- dplyr::recode(as.character(df$HEAD_SEX),        !!!sex_lbl,    .default="?")
    if ("HEAD_POPULATION" %in% names(df))
      df$HEAD_POPULATION <- dplyr::recode(as.character(df$HEAD_POPULATION), !!!pop_lbl,    .default="?")
    if ("PRESENT_STATUS"  %in% names(df))
      df$PRESENT_STATUS  <- dplyr::recode(as.character(df$PRESENT_STATUS),  !!!status_lbl, .default="?")
    if ("ENG_ACCESS"      %in% names(df))
      df$ENG_ACCESS      <- dplyr::recode(as.character(df$ENG_ACCESS),      !!!elec_lbl,   .default="?")

    names(df) <- gsub("_", " ", names(df))

    num_cols <- c("INCOME","EXPENDITURE","INCOME PCP","EXPENDITURE PCP")
    num_cols <- num_cols[num_cols %in% names(df)]

    dt <- datatable(df,
      rownames  = FALSE,
      filter    = "top",
      class     = "cell-border",
      options   = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "lftip",
        columnDefs = list(list(className="dt-center", targets="_all"))
      )
    )

    if (length(num_cols) > 0)
      dt <- formatCurrency(dt, num_cols,
                           currency="R ", interval=3, mark=",", digits=0)

    dt
  })

  # ── Plot: histogram ──────────────────────────────────────────────────────────
  output$plt_hist <- renderPlot({
    req(filt()); df <- filt()
    cap  <- quantile(df$INCOME, 0.95, na.rm=TRUE)
    med  <- density_mode(df$INCOME[df$INCOME <= cap], na.rm=TRUE)
    df2  <- df[df$INCOME <= cap & !is.na(df$INCOME), ]

    ggplot(df2, aes(x=INCOME)) +
      geom_histogram(fill="#A0522D", colour="#3D2B1F", bins=45, alpha=0.9) +
      geom_vline(xintercept=med, colour="#D4783E",
                 linewidth=1.3, linetype="dashed") +
      annotate("text", x=med*1.05, y=Inf, vjust=1.4, hjust=0,
               label=paste0("Mode\nR", format(round(med), big.mark=",")),
               colour="#D4783E", size=3, family="serif", fontface="bold") +
      scale_x_continuous(labels=label_dollar(prefix="R ", big.mark=",", scale=1e-3,
                                              suffix="k")) +
      scale_y_continuous(labels=comma) +
      labs(x="Annual Household Income", y="Households",
           subtitle="Capped at 95th percentile") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: decile bar ─────────────────────────────────────────────────────────
  output$plt_decile <- renderPlot({
    req(filt())
    filt() %>%
      group_by(Decile = factor(INCOME_DECILE)) %>%
      summarise(Mode=density_mode(INCOME, na.rm=TRUE), n=n(), .groups="drop") %>%
      ggplot(aes(x=Decile, y=Mode, fill=Decile)) +
      geom_col(colour="#3D2B1F", width=0.72, alpha=0.9) +
      geom_text(aes(label=paste0("R",format(round(Mode/1000), big.mark=","),"k")),
                vjust=-0.4, size=2.8, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=PAL, guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k"),
                         expand=expansion(mult=c(0,0.12))) +
      labs(x="Household Income Decile", y="Mode Annual Household Income") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: population group boxplot ───────────────────────────────────────────
  output$plt_pop_box <- renderPlot({
    req(filt()); df <- filt()
    cap <- quantile(df$INCOME, 0.95, na.rm=TRUE)
    df  %>%
      filter(INCOME <= cap) %>%
      mutate(Pop = dplyr::recode(as.character(HEAD_POPULATION),
                                 !!!pop_lbl, .default="Other")) %>%
      ggplot(aes(x=reorder(Pop, INCOME, FUN=density_mode), y=INCOME, fill=Pop)) +
      geom_boxplot(colour="#3D2B1F", outlier.colour="#D4783E",
                   outlier.size=0.6, alpha=0.85, width=0.6) +
      scale_fill_manual(
        values=c("Black African"="#A0522D","Coloured"="#D4783E",
                 "Indian/Asian"="#C4956A","White"="#8B7D6B","Other"="#6B4226"),
        guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k")) +
      coord_flip() +
      labs(x=NULL, y="Annual Household Income", subtitle="Capped at 95th percentile") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: income/capita by household size ────────────────────────────────────
  output$plt_hsize <- renderPlot({
    req(filt())
    filt() %>%
      filter(HSIZE <= 12, !is.na(INCOME_PCP)) %>%
      group_by(Size = HSIZE) %>%
      summarise(med=density_mode(INCOME_PCP, na.rm=TRUE), n=n(), .groups="drop") %>%
      ggplot(aes(x=Size, y=med)) +
      geom_area(fill="#A0522D30") +
      geom_line(colour="#A0522D", linewidth=1.4) +
      geom_point(aes(size=n), colour="#D4783E", fill="#FAF3EA",
                 shape=21, stroke=1.8) +
      scale_x_continuous(breaks=1:12) +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k")) +
      scale_size_continuous(range=c(3,11), name="N households") +
      labs(x="Household Size", y="Mode Household Income Per Capita",
           subtitle="Households with size 1–12 | point size = n households") +
      theme_mcm() +
      theme(legend.position="bottom")
  }, bg="#F5E6D3")

  # ── Plot: status bar ─────────────────────────────────────────────────────────
  output$plt_status <- renderPlot({
    req(filt())
    lvls <- c("Wealthy","Very comfortable","Reasonably comfortable",
              "Just getting along","Poor","Very poor")
    clrs <- c("#3D2B1F","#6B4226","#A0522D","#C4956A","#D4783E","#EDD9C0")

    filt() %>%
      mutate(Status = dplyr::recode(as.character(PRESENT_STATUS),
                                    !!!status_lbl, .default="Unknown")) %>%
      filter(Status %in% lvls) %>%
      group_by(Status) %>%
      summarise(med=density_mode(INCOME, na.rm=TRUE), n=n(), .groups="drop") %>%
      mutate(Status = factor(Status, levels=lvls)) %>%
      ggplot(aes(x=Status, y=med, fill=Status)) +
      geom_col(colour="#3D2B1F", width=0.65, alpha=0.92) +
      geom_text(aes(label=paste0("n=", format(n, big.mark=","))),
                vjust=-0.5, size=3, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=setNames(clrs, lvls), guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k"),
                         expand=expansion(mult=c(0,0.14))) +
      labs(x=NULL, y="Mode Annual Household Income") +
      theme_mcm() +
      theme(axis.text.x=element_text(angle=22, hjust=1, size=10))
  }, bg="#F5E6D3")

  # ── Plot: income by province (boxplot) ───────────────────────────────────────
  output$plt_province_box <- renderPlot({
    req(filt()); df <- filt()
    cap <- quantile(df$INCOME, 0.95, na.rm=TRUE)
    df %>%
      filter(INCOME <= cap, !is.na(PROVINCE)) %>%
      mutate(Province = dplyr::recode(as.character(PROVINCE),
                                      !!!province_lbl, .default="Other")) %>%
      ggplot(aes(x=reorder(Province, INCOME, FUN=density_mode), y=INCOME, fill=Province)) +
      geom_boxplot(colour="#3D2B1F", outlier.colour="#D4783E",
                   outlier.size=0.6, alpha=0.85, width=0.6) +
      scale_fill_manual(values=PAL, guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k")) +
      coord_flip() +
      labs(x=NULL, y="Annual Household Income", subtitle="Capped at 95th percentile") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: mode income by settlement type (bar) ───────────────────────────────
  output$plt_settlement_bar <- renderPlot({
    req(filt())
    filt() %>%
      filter(!is.na(SETTLEMENT_TYPE)) %>%
      mutate(Settlement = dplyr::recode(as.character(SETTLEMENT_TYPE),
                                        !!!settlement_lbl, .default="Other")) %>%
      filter(Settlement %in% names(SETTLE_PAL)) %>%
      mutate(Settlement = factor(Settlement, levels=names(SETTLE_PAL))) %>%
      group_by(Settlement) %>%
      summarise(Mode=density_mode(INCOME, na.rm=TRUE), n=n(), .groups="drop") %>%
      ggplot(aes(x=Settlement, y=Mode, fill=Settlement)) +
      geom_col(colour="#3D2B1F", width=0.6, alpha=0.9) +
      geom_text(aes(label=paste0("R",format(round(Mode/1000), big.mark=","),"k")),
                vjust=-0.4, size=3, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=SETTLE_PAL, guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k"),
                         expand=expansion(mult=c(0,0.14))) +
      labs(x=NULL, y="Mode Annual Household Income") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Table: province summary ──────────────────────────────────────────────────
  output$tbl_province <- renderTable({
    req(filt())
    filt() %>%
      filter(!is.na(PROVINCE)) %>%
      mutate(Province = dplyr::recode(as.character(PROVINCE),
                                      !!!province_lbl, .default="Other")) %>%
      group_by(Province) %>%
      summarise(
        `N`               = n(),
        `Mode Household Income (R)` = format(round(density_mode(INCOME, na.rm=TRUE)), big.mark=","),
        `Mean Household Income (R)` = format(round(mean(INCOME, na.rm=TRUE)), big.mark=","),
        .groups = "drop"
      ) %>%
      arrange(desc(N))
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="lccc",
  width="100%")

  # ── Table: settlement type summary ───────────────────────────────────────────
  output$tbl_settlement <- renderTable({
    req(filt())
    filt() %>%
      filter(!is.na(SETTLEMENT_TYPE)) %>%
      mutate(Settlement = dplyr::recode(as.character(SETTLEMENT_TYPE),
                                        !!!settlement_lbl, .default="Other")) %>%
      filter(Settlement %in% names(SETTLE_PAL)) %>%
      mutate(Settlement = factor(Settlement, levels=names(SETTLE_PAL))) %>%
      group_by(Settlement) %>%
      summarise(
        `N`                    = n(),
        `Mode Household Income (R)`      = format(round(density_mode(INCOME, na.rm=TRUE)), big.mark=","),
        `Mean Household Income (R)`      = format(round(mean(INCOME, na.rm=TRUE)), big.mark=","),
        `Mean Expenditure (R)` = format(round(mean(EXPENDITURE, na.rm=TRUE)), big.mark=","),
        .groups = "drop"
      )
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="lcccc",
  width="100%")

  # ── Plot: financial & lifestyle indicators overview ──────────────────────────
  output$plt_indicators <- renderPlot({
    req(filt()); df <- filt()
    pct_df <- data.frame(
      Indicator = names(ind_vars),
      PctYes    = sapply(ind_vars, function(col) 100 * mean(df[[col]] == 1, na.rm=TRUE))
    )
    pct_df$Indicator <- factor(pct_df$Indicator, levels=pct_df$Indicator[order(pct_df$PctYes)])

    ggplot(pct_df, aes(x=Indicator, y=PctYes, fill=Indicator)) +
      geom_col(colour="#3D2B1F", width=0.68, alpha=0.9) +
      geom_text(aes(label=paste0(round(PctYes,1), "%")),
                hjust=-0.15, size=3.2, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=PAL, guide="none") +
      scale_y_continuous(labels=label_percent(scale=1),
                         limits=c(0, max(pct_df$PctYes)*1.18),
                         expand=expansion(mult=c(0,0))) +
      coord_flip() +
      labs(x=NULL, y="Share of Households Answering “Yes”") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Cards: financial & lifestyle indicators at a glance ─────────────────────
  output$ind_cards <- renderUI({
    req(filt()); df <- filt()
    styles <- rep(c("default","dark","tan"), length.out = length(ind_vars))

    cards <- lapply(seq_along(ind_vars), function(i) {
      lbl <- names(ind_vars)[i]
      col <- ind_vars[[i]]
      pct <- 100 * mean(df[[col]] == 1, na.rm=TRUE)
      column(4,
        stat_card(paste(ind_icons[[lbl]], lbl),
                  paste0(round(pct, 1), "%"), styles[i])
      )
    })
    fluidRow(cards)
  })
}

shinyApp(ui, server)
