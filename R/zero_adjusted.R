# zero_adjusted_discrete <- function(distrib, link_za = logit_link()) {
#   # --- 1. Validation ---
#   if (!inherits(distrib, "distrib")) stop("Input must be a 'distrib' object.")
#   if (distrib$type != "discrete") stop("This function requires a discrete distribution.")
#   if (distrib$bounds[1] > 0) stop("Original distribution must include 0 in support.")
#
#   o <- list()
#   class(o) <- "distrib"
#
#   # --- 2. Metadata ---
#   o$distrib_name <- paste0("zero-adjusted_", distrib$distrib_name)
#   o$type <- "discrete"
#   o$dimension <- distrib$dimension
#   o$bounds <- c(0, distrib$bounds[2])
#
#   o$params <- c(distrib$params, "za")
#   o$n_params <- distrib$n_params + 1
#   o$params_interpretation <- c(distrib$params_interpretation, za = "prob. zero")
#
#   o$params_bounds <- distrib$params_bounds
#   o$params_bounds$za <- c(0, 1)
#
#   o$link_params <- distrib$link_params
#   o$link_params$za <- link_za
#
#   split_theta <- function(theta) {
#     n <- length(theta)
#     list(orig = theta[1:(n - 1)], za = theta[[n]])
#   }
#
#   # --- 3. PMF ---
#   # P(Y=0) = za
#   # P(Y=y) = (1-za) * f(y) / (1 - f(0))  per y > 0
#   o$pdf <- function(y, theta, log = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     # Pre-calcolo f(0) e log(1-f(0)) per normalizzare
#     pdf0 <- distrib$pdf(0, pars$orig, log = FALSE)
#     log_den <- log1p(-pdf0) # log(1 - f(0)) più stabile
#
#     # Log-PMF per parte positiva
#     log_res_pos <- log(1 - za) + distrib$pdf(y, pars$orig, log = TRUE) - log_den
#
#     # Gestione vettoriale
#     n <- length(y)
#     log_res <- numeric(n)
#
#     is_zero <- (y == 0)
#     if (any(is_zero)) log_res[is_zero] <- log(za)
#     if (any(!is_zero)) log_res[!is_zero] <- log_res_pos[!is_zero]
#
#     if (log) log_res else exp(log_res)
#   }
#
#   # --- 4. CDF ---
#   o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     f0 <- distrib$pdf(0, pars$orig, log = FALSE)
#     F_orig <- distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)
#
#     # F_za(q) = za + (1-za) * (F(q) - f0) / (1 - f0)  se q >= 0
#     # Nota: (F(q) - f0) è la massa accumulata da 1 a q
#
#     # Numeratore normalizzato
#     cdf_trunc <- (F_orig - f0) / (1 - f0)
#     cdf_trunc <- pmax(0, cdf_trunc) # Safety
#
#     res <- za + (1 - za) * cdf_trunc
#
#     # Se q < 0, CDF è 0 (gestito implicitamente se F_orig da 0, ma forziamo)
#     res[q < 0] <- 0
#     res <- pmin(res, 1)
#
#     if (!lower.tail) res <- 1 - res
#     if (log.p) log(res) else res
#   }
#
#   # --- 5. Quantile Function ---
#   o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     if (log.p) p <- exp(p)
#     if (!lower.tail) p <- 1 - p
#     p <- pmin(pmax(p, 0), 1)
#
#     # Hurdle Quantile Logic:
#     # Se p <= za, restituisci 0.
#     # Se p > za, inverti la truncated: Q_trunc( (p - za)/(1 - za) )
#
#     # Per invertire la truncated count Q_trunc(u):
#     # Cerca x tale che F_trunc(x) = u
#     # (F(x) - f0)/(1-f0) = u  => F(x) = u(1-f0) + f0
#
#     f0 <- distrib$pdf(0, pars$orig, log = FALSE)
#
#     q_vals <- numeric(length(p))
#     is_pos <- (p > za)
#
#     if (any(is_pos)) {
#       # Probabilità trasformata per la CDF originale
#       # u = (p - za) / (1 - za)
#       # target_prob = u*(1-f0) + f0
#
#       p_curr <- p[is_pos]
#       za_curr <- if (length(za) > 1) za[is_pos] else za
#       f0_curr <- if (length(f0) > 1) f0[is_pos] else f0
#       pars_orig_curr <- lapply(pars$orig, function(x) {
#         if (length(x) > 1) x[is_pos] else x
#       })
#
#       u <- (p_curr - za_curr) / (1 - za_curr)
#       target_prob <- u * (1 - f0_curr) + f0_curr
#
#       # numerical safety clamp
#       target_prob <- pmin(target_prob, 1 - 1e-10)
#
#       q_vals[is_pos] <- distrib$quantile(target_prob, pars_orig_curr)
#     }
#     q_vals
#   }
#
#   # --- 6. Random Number Generator ---
#   o$rng <- function(n, theta) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     # 1. Decidi chi è zero (Bernoulli)
#     is_zero <- stats::runif(n) < za
#
#     # Inizializza vettore risultato
#     y <- numeric(n)
#
#     # Se ci sono valori positivi da generare...
#     if (any(!is_zero)) {
#       n_pos <- sum(!is_zero)
#
#       # --- SUBSETTING CRUCIALE ---
#       # Dobbiamo usare i parametri corretti per le osservazioni che NON sono zero
#
#       # 1. Calcola f0 per TUTTI, poi prendi solo i !is_zero
#       f0_all <- distrib$pdf(0, pars$orig, log = FALSE)
#       f0_curr <- if (length(f0_all) > 1) f0_all[!is_zero] else f0_all
#
#       # 2. Subset parametri originali
#       pars_orig_curr <- lapply(pars$orig, function(x) {
#         if (length(x) > 1) x[!is_zero] else x
#       })
#
#       # --- GENERAZIONE TRONCATA ---
#       # Inverse Transform Sampling per la troncata:
#       # Genera u in [0, 1] e scalalo nell'intervallo [f0, 1] della CDF originale
#       u <- stats::runif(n_pos)
#       u_scaled <- f0_curr + u * (1 - f0_curr)
#
#       # Inverti la CDF originale usando i parametri corretti
#       y[!is_zero] <- distrib$quantile(u_scaled, pars_orig_curr)
#     }
#
#     y
#   }
#
#   o$loglik <- function(y, theta) {
#     o$pdf(y, theta, log = TRUE)
#   }
#
#   # --- 7. Gradient ---
#   o$gradient <- function(y, theta, par = NULL) {
#     if (is.null(par)) par <- o$params
#     pars <- split_theta(theta)
#     za <- pars$za
#     par_or <- par[par %in% distrib$params]
#
#     # Calcoli ausiliari su f(0)
#     f0 <- distrib$pdf(0, pars$orig, log = FALSE)
#     score_0 <- distrib$gradient(0, pars$orig, par_or) # Score S(0) = f'(0)/f(0)
#
#     # Termine di correzione per il troncamento:
#     # d/dTheta [ -log(1 - f(0)) ] = f'(0) / (1 - f(0))
#     # = (f(0) * S(0)) / (1 - f(0))
#     correction_factor <- f0 / (1 - f0)
#
#     grad_orig <- distrib$gradient(y, pars$orig, par_or)
#     res_grad <- list()
#
#     # 1. Gradiente rispetto a Theta (Parametri Originali)
#     # Se y=0: Gradiente è 0 (Theta non impatta la probabilità za)
#     # Se y>0: Score(y) + Correction
#     for (nm in names(grad_orig)) {
#       if (nm %in% par) {
#         # grad_orig[[nm]] è lo Score S(y)
#         # score_0[[nm]] è lo Score S(0)
#
#         term_pos <- grad_orig[[nm]] + correction_factor * score_0[[nm]]
#         res_grad[[nm]] <- ifelse(y == 0, 0, term_pos)
#       }
#     }
#
#     # 2. Gradiente rispetto a ZA
#     if ("za" %in% par) {
#       res_grad$za <- ifelse(y == 0, 1 / za, -1 / (1 - za))
#     }
#
#     res_grad
#   }
#
#   # --- 8. Hessian ---
#   o$hessian <- function(y, theta, expected = FALSE) {
#     pars <- split_theta(theta)
#     th_orig <- pars$orig
#     za <- pars$za
#
#     # --- 1. PRE-CALCOLI SU THETA (Indipendenti da y) ---
#     # Questi termini servono sia per expected che per observed
#
#     # PDF e Derivate in 0
#     f0 <- distrib$pdf(0, th_orig, log = FALSE)
#     grad_0 <- distrib$gradient(0, th_orig) # f'(0)/f(0)
#     hess_0_obs <- distrib$hessian(0, th_orig, expected = FALSE) # H(0)
#
#     # Costanti utili
#     denom <- 1 - f0
#     denom2 <- denom^2
#
#     # Prepariamo la struttura per i risultati
#     res_hess <- list()
#
#     # --- 2. BLOCCHI ZA (Ortogonali) ---
#     if (expected) {
#       # Expected: Fisher Info (negativa)
#       res_hess[["za_za"]] <- -1 / (za * (1 - za))
#     } else {
#       # Observed
#       val_0 <- -1 / (za^2)
#       val_pos <- -1 / ((1 - za)^2)
#       res_hess[["za_za"]] <- ifelse(y == 0, val_0, val_pos)
#     }
#
#     # I termini misti sono sempre 0
#     for (nm in names(th_orig)) {
#       res_hess[[paste0(nm, "_za")]] <- rep(0, length(y))
#     }
#
#     # --- 3. BLOCCO THETA-THETA (Il cuore del calcolo) ---
#
#     # Iteriamo sui parametri originali (es. mu_mu, theta_theta, mu_theta)
#     for (nm in names(hess_0_obs)) {
#       parts <- strsplit(nm, "_")[[1]]
#       p1 <- parts[1]
#       p2 <- parts[length(parts)]
#
#       # A. Calcolo del termine di correzione H_corr
#       #    H_corr = [ (1-f)*f'' + (f')^2 ] / (1-f)^2
#
#       # Recuperiamo f' e f'' dai gradienti e hessiani in 0
#       s1 <- grad_0[[p1]]
#       s2 <- grad_0[[p2]]
#       h_log_0 <- hess_0_obs[[nm]]
#
#       f_prime_1 <- f0 * s1
#       f_prime_2 <- f0 * s2
#       f_prime_prime <- f0 * (h_log_0 + s1 * s2)
#
#       # Formula corretta con segno positivo per entrambi i termini al numeratore
#       hess_correction <- ((1 - f0) * f_prime_prime + f_prime_1 * f_prime_2) / denom2
#
#
#       # B. Branching netto tra Expected e Observed
#       if (expected) {
#         # --- CASO EXPECTED ---
#         # E_za[H] = (1-za) * E_trunc[H_obs_pos]
#         # H_obs_pos = H_orig(y) + H_corr
#         # E_trunc[H_orig] = (E_orig[H] - f0*H_orig(0)) / (1-f0)
#
#         # Hessiana Attesa Originale (Fisher Info completa)
#         I_orig_full <- distrib$hessian(y, th_orig, expected = TRUE)[[nm]]
#
#         # Aspettazione troncata di H_orig
#         E_trunc_H_orig <- (I_orig_full - f0 * h_log_0) / denom
#
#         # Totale pesato per la probabilità di non-zero (1-za)
#         # Nota: hess_correction è costante, quindi E[H_corr] = H_corr
#         val_expected <- (1 - za) * (E_trunc_H_orig + hess_correction)
#
#         res_hess[[nm]] <- val_expected
#       } else {
#         # --- CASO OBSERVED ---
#         # Se y=0: 0
#         # Se y>0: H_orig(y) + H_corr
#
#         h_orig_y <- distrib$hessian(y, th_orig, expected = FALSE)[[nm]]
#
#         # SOMMIAMO la correzione (non sottraiamo!)
#         term_pos <- h_orig_y + hess_correction
#
#         res_hess[[nm]] <- ifelse(y == 0, 0, term_pos)
#       }
#     }
#
#     expand_params(res_hess[hess_names(o$params)], length(y))
#   }
#
#   # --- 9. Moments ---
#   # Rinormalizzazione corretta per troncamento
#   o$raw_moment <- function(n, theta) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     # E_orig[Y^n]
#     raw_orig <- moment(distrib, pars$orig, p = n, central = FALSE)
#     f0 <- distrib$pdf(0, pars$orig, log = FALSE)
#
#     # E_trunc[Y^n] = E_orig[Y^n] / (1 - f0)  (poiché 0^n = 0 non contribuisce alla somma originale)
#     trunc_moment <- raw_orig / (1 - f0)
#
#     # E_za[Y^n] = za*0 + (1-za)*E_trunc
#     (1 - za) * trunc_moment
#   }
#
#   o$mean <- function(theta) o$raw_moment(1, theta)
#
#   o$variance <- function(theta) {
#     m1 <- o$raw_moment(1, theta)
#     m2 <- o$raw_moment(2, theta)
#     m2 - m1^2
#   }
#
#   o$skewness <- function(theta) {
#     m1 <- o$raw_moment(1, theta)
#     m2 <- o$raw_moment(2, theta)
#     m3 <- o$raw_moment(3, theta)
#     sigma <- sqrt(m2 - m1^2)
#     (m3 - 3 * m1 * m2 + 2 * m1^3) / (sigma^3)
#   }
#
#   o$kurtosis <- function(theta) {
#     m1 <- o$raw_moment(1, theta)
#     m2 <- o$raw_moment(2, theta)
#     m3 <- o$raw_moment(3, theta)
#     m4 <- o$raw_moment(4, theta)
#     var_val <- m2 - m1^2
#     num <- m4 - 4 * m1 * m3 + 6 * (m1^2) * m2 - 3 * (m1^4)
#     (num / (var_val^2)) - 3
#   }
#
#   o$initialize <- function(y) {
#     prop_zeros <- mean(y == 0)
#     za_init <- min(max(prop_zeros, 0.01), 0.99)
#     y_pos <- y[y > 0]
#     if (length(y_pos) == 0) y_pos <- y
#     res <- distrib$initialize(y_pos)
#     res$za <- za_init
#     res
#   }
#
#   o
# }
#
#
#
#
#
# zero_adjusted_continuous <- function(distrib, link_za = logit_link()) {
#   # --- 1. Validation ---
#   if (!inherits(distrib, "distrib")) stop("Input must be a 'distrib' object.")
#   if (distrib$type != "continuous") stop("This function requires a continuous distribution.")
#
#   o <- list()
#   class(o) <- "distrib"
#
#   # --- 2. Metadata ---
#   o$distrib_name <- paste0("zero-adjusted_", distrib$distrib_name)
#   o$type <- "continuous" # Remains continuous (with a spike)
#   o$dimension <- distrib$dimension
#
#   # Bounds: Expand to include 0 if not present (e.g. for LogNormal)
#   o$bounds <- c(min(0, distrib$bounds[1]), distrib$bounds[2])
#
#   o$params <- c(distrib$params, "za")
#   o$n_params <- distrib$n_params + 1
#   o$params_interpretation <- c(distrib$params_interpretation, za = "prob. zero")
#
#   o$params_bounds <- distrib$params_bounds
#   o$params_bounds$za <- c(0, 1)
#
#   o$link_params <- distrib$link_params
#   o$link_params$za <- link_za
#
#   split_theta <- function(theta) {
#     n <- length(theta)
#     list(orig = theta[1:(n - 1)], za = theta[[n]])
#   }
#
#   # --- 3. PDF ---
#   # P(Y=0) = za (Dirac Delta mass)
#   # f(y)   = (1-za) * f_orig(y)   for y != 0
#   o$pdf <- function(y, theta, log = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     # Calculate original log-pdf
#     log_f_orig <- distrib$pdf(y, pars$orig, log = TRUE)
#
#     n <- length(y)
#     log_res <- numeric(n)
#
#     is_zero <- (y == 0)
#
#     # 1. Mass at zero
#     if (any(is_zero)) {
#       log_res[is_zero] <- log(za)
#     }
#
#     # 2. Continuous density scaled by (1-za)
#     if (any(!is_zero)) {
#       # log( (1-za) * f(y) ) = log(1-za) + log(f(y))
#       log_res[!is_zero] <- log(1 - za) + log_f_orig[!is_zero]
#     }
#
#     if (log) log_res else exp(log_res)
#   }
#
#   # --- 4. CDF ---
#   # F_za(q) = (1-za)F_orig(q) + za * I(q >= 0)
#   o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     F_orig <- distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)
#
#     # Base component scaled
#     res <- (1 - za) * F_orig
#
#     # Add jump at 0
#     # Note: For Gaussian, F_orig(0)=0.5. So jump happens from (1-za)*0.5 to (1-za)*0.5 + za.
#     res[q >= 0] <- res[q >= 0] + za
#
#     res <- pmin(pmax(res, 0), 1)
#
#     if (!lower.tail) res <- 1 - res
#     if (log.p) log(res) else res
#   }
#
#   # --- 5. Quantile Function ---
#   o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     if (log.p) p <- exp(p)
#     if (!lower.tail) p <- 1 - p
#     p <- pmin(pmax(p, 0), 1)
#
#     # Calculate critical points for the jump at 0
#     F0_orig <- distrib$cdf(0, pars$orig) # E.g., 0.5 for Gaussian, 0 for Gamma
#
#     p_lower <- (1 - za) * F0_orig
#     p_upper <- p_lower + za
#
#     q_vals <- numeric(length(p))
#
#     # Case A: Left Tail (e.g. Gaussian negatives) -> p < p_lower
#     idx_left <- (p < p_lower)
#     if (any(idx_left)) {
#       # (1-za)F(q) = p  =>  F(q) = p / (1-za)
#       p_trans <- p[idx_left] / (1 - za)
#
#       # Parametri subsetting (per evitare warning di lunghezza)
#       pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[idx_left] else x)
#       q_vals[idx_left] <- distrib$quantile(p_trans, pars_sub)
#     }
#
#     # Case B: The Jump (Zero Mass) -> p_lower <= p <= p_upper
#     # q_vals is already 0 initialized, so we do nothing.
#
#     # Case C: Right Tail (Positives) -> p > p_upper
#     idx_right <- (p > p_upper)
#     if (any(idx_right)) {
#       # (1-za)F(q) + za = p  =>  F(q) = (p - za) / (1 - za)
#       p_curr <- p[idx_right]
#       za_curr <- if (length(za) > 1) za[idx_right] else za
#
#       p_trans <- (p_curr - za_curr) / (1 - za_curr)
#       p_trans <- pmin(p_trans, 1) # Safety clamp
#
#       pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[idx_right] else x)
#       q_vals[idx_right] <- distrib$quantile(p_trans, pars_sub)
#     }
#
#     q_vals
#   }
#
#   # --- 6. RNG ---
#   o$rng <- function(n, theta) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     is_zero <- stats::runif(n) < za
#     y <- numeric(n)
#
#     if (any(!is_zero)) {
#       # Per la parte continua non serve Inverse Sampling complesso o troncamento.
#       # Basta generare dalla distribuzione originale.
#       # MA: Se la distribuzione originale può generare negativi (Gaussiana), è tutto OK.
#       # Se può generare solo positivi (Gamma), è tutto OK.
#       # L'unica cosa: il modello teorico dice che se Y != 0, Y ~ f_orig.
#
#       pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[!is_zero] else x)
#       y[!is_zero] <- distrib$rng(sum(!is_zero), pars_sub)
#     }
#     y
#   }
#
#   o$loglik <- function(y, theta) {
#     o$pdf(y, theta, log = TRUE)
#   }
#
#   # --- 7. Gradient ---
#   # Molto più semplice del discreto!
#   o$gradient <- function(y, theta, par = NULL) {
#     if (is.null(par)) par <- o$params
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     res_grad <- list()
#
#     # A. Gradiente Theta (Parametri Originali)
#     # Se y=0: 0
#     # Se y!=0: Gradiente originale standard
#     grad_orig <- distrib$gradient(y, pars$orig)
#
#     for (nm in names(grad_orig)) {
#       if (nm %in% par) {
#         val <- grad_orig[[nm]]
#         res_grad[[nm]] <- ifelse(y == 0, 0, val)
#       }
#     }
#
#     # B. Gradiente ZA
#     if ("za" %in% par) {
#       res_grad$za <- ifelse(y == 0, 1 / za, -1 / (1 - za))
#     }
#
#     res_grad
#   }
#
#   # --- 8. Hessian ---
#   # Estremamente pulita: niente termini correttivi, niente misti.
#   o$hessian <- function(y, theta, expected = FALSE) {
#     pars <- split_theta(theta)
#     za <- pars$za
#
#     res_hess <- list()
#
#     # A. Blocco ZA-ZA
#     if (expected) {
#       res_hess[["za_za"]] <- -1 / (za * (1 - za))
#     } else {
#       val_0 <- -1 / (za^2)
#       val_pos <- -1 / ((1 - za)^2)
#       res_hess[["za_za"]] <- ifelse(y == 0, val_0, val_pos)
#     }
#
#     # B. Blocco Misto (Sempre 0)
#     for (nm in names(pars$orig)) {
#       res_hess[[paste0(nm, "_za")]] <- rep(0, length(y))
#     }
#
#     # C. Blocco Theta-Theta
#     if (expected) {
#       # E_za[H] = (1-za) * E_orig[H]
#       # Non c'è troncamento, quindi E[H | y!=0] = E_orig[H]
#       h_exp_orig <- distrib$hessian(y, pars$orig, expected = TRUE)
#
#       for (nm in names(h_exp_orig)) {
#         res_hess[[nm]] <- (1 - za) * h_exp_orig[[nm]]
#       }
#     } else {
#       # Observed:
#       # Se y=0: 0
#       # Se y!=0: H_orig(y)
#       h_obs_orig <- distrib$hessian(y, pars$orig, expected = FALSE)
#
#       for (nm in names(h_obs_orig)) {
#         res_hess[[nm]] <- ifelse(y == 0, 0, h_obs_orig[[nm]])
#       }
#     }
#
#     expand_params(res_hess[hess_names(o$params)], length(y))
#   }
#
#   # --- 9. Moments ---
#   o$raw_moment <- function(n, theta) {
#     pars <- split_theta(theta)
#     # E[Y^n] = za*0 + (1-za)*E_orig[Y^n]
#     # Semplice scalatura
#     (1 - pars$za) * moment(distrib, pars$orig, p = n, central = FALSE)
#   }
#
#   o$mean <- function(theta) o$raw_moment(1, theta)
#
#   o$variance <- function(theta) {
#     m1 <- o$raw_moment(1, theta)
#     m2 <- o$raw_moment(2, theta)
#     m2 - m1^2
#   }
#
#   o$initialize <- function(y) {
#     prop_zeros <- mean(y == 0)
#     za_init <- min(max(prop_zeros, 0.01), 0.99)
#     # Inizializza sui non-zero
#     y_pos <- y[y != 0]
#     if (length(y_pos) == 0) y_pos <- y
#     res <- distrib$initialize(y_pos)
#     res$za <- za_init
#     res
#   }
#
#   o
# }
#
#
# zero_adjusted <- function(distrib, link_za = logit_link()) {
#   if (!inherits(distrib, "distrib")) {
#     stop("Input must be a 'distrib' object.")
#   }
#
#   if (distrib$type == "discrete") {
#     zero_adjusted_discrete(distrib, link_za)
#   } else if (distrib$type == "continuous") {
#     zero_adjusted_continuous(distrib, link_za)
#   }
# }
