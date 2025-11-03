############################################################
# TP2
############################################################

# Chargement des packages

library(actuar)
library(ggplot2)
library(dplyr)

## === EXERCICE 2.1 ========================================

# --- Paramètres ---
lambda <- 5
q <- 0.3
n_sim <- 1e4

# --- Simulation empirique ---
set.seed(123)
N <- rpois(n_sim, lambda)
S_sim <- sapply(N, function(n) sum(rbinom(n, 1, q)))

# --- Panjer manuel ---

# Paramètre
k_max <- 2  # calcul de P(S=0), P(S=1), P(S=2)

# Initialisation
p <- numeric(k_max + 1)
f <- c(1-q, q)  # Bernoulli : f[0]=P(X=0), f[1]=P(X=1)

# p0
p[1] <- exp(-lambda * q)  # P(S=0)

# Récursion de Panjer
for (k in 1:k_max) {
  p[k+1] <- (lambda * q / k) * p[k]  # formule de Panjer pour Poisson - Bernoulli
}

names(p) <- 0:k_max
p

# --- Panjer avec actuar ---
sev_probs <- c(1 - q, q)

Fs <- aggregateDist("recursive",
                    model.freq = "poisson", lambda = lambda,
                    model.sev = sev_probs,
                    x.scale = 1)

# --- Comparaisons graphiques ---
par(mfrow = c(1, 2))

# Densité (loi empirique vs Panjer)
hist(S_sim, freq = FALSE, breaks = seq(-0.5, max(S_sim)+0.5, 1),
     main = "Loi empirique vs théorique", xlab = "S")
points(0:8, diff(Fs(-1:8)), col = "red", pch = 19)
legend("topright", legend = c("Empirique", "Panjer"),
       col = c("grey", "red"), pch = 19)

# Fonction de répartition (empirique vs Panjer)
plot(ecdf(S_sim), do.points = FALSE, main = "Fonction de répartition",
     ylab = "F(x)", xlab = "x")
lines(0:8, Fs(0:8), col = "red", lwd = 2)

# --- Probabilités demandées ---
P0 <- Fs(0)
P1 <- Fs(1) - Fs(0)
P2 <- Fs(2) - Fs(1)
cat("P[S=0] =", P0, "\nP[S=1] =", P1, "\nP[S=2] =", P2, "\n")




## === EXERCICE 2.2 ========================================

# --- Paramètres ---
n <- 10
p <- 0.1
alpha <- 1/2
tau <- 1/3
set.seed(123)

# --- Simulation empirique ---
n_sim <- 1e5
N <- rbinom(n_sim, n, p)
S_sim <- sapply(N, function(m) sum(rgamma(m, shape = alpha, rate = tau)))

# --- Fonctions paramétrées ---
Fx_gamma  <- function(x) pgamma(x, shape = alpha, rate = tau)
lev_gamma <- function(x) levgamma(x, shape = alpha, rate = tau)

# --- Discrétisation sans biais ---
q_hi <- qgamma(0.9995, shape = alpha, rate = tau)
step <- 0.05
fx_discr <- discretize(Fx_gamma, from = 0, to = q_hi, step = step,
                       method = "unbiased", lev = lev_gamma)
m0 <- fx_discr[1]
fx_discr[1] <- 0
fx_discr[2] <- fx_discr[2] + m0
fx_discr <- fx_discr / sum(fx_discr)

# --- Distribution agrégée ---
Fs <- aggregateDist("recursive",
                    model.freq = "binomial", size = n, prob = p,
                    model.sev = fx_discr,
                    x.scale = step)

# --- Comparaison graphique ---
plot(ecdf(S_sim), do.points = FALSE,
     main = "F(S) empirique vs Panjer", xlab = "S", ylab = "F(S)")
curve(Fs(x), add = TRUE, col = "red", lwd = 2)
legend("bottomright", legend = c("Empirique", "Panjer"),
       col = c("black", "red"), lwd = 2)

# --- Probabilités ---
cuts <- seq(0, ceiling(max(S_sim)), by = 1)  # borne supérieure dynamique
intervals <- paste0("[", cuts[-length(cuts)], ",", cuts[-1], ")")

probs_Panjer <- diff(Fs(cuts))
probs_emp  <- as.numeric(table(cut(S_sim, breaks = cuts))) / length(S_sim)

x_grid <- seq(0, q_hi, by = step)
dens_values <- diff(Fs(x_grid))
idx_1 <- which.min(abs(x_grid[-1] - 1))
idx_2 <- which.min(abs(x_grid[-1] - 2))
P_S_eq1_Panjer <- dens_values[idx_1]
P_S_eq2_Panjer <- dens_values[idx_2]
P_S_eq1_emp <- mean(abs(S_sim - 1) < step/2)
P_S_eq2_emp <- mean(abs(S_sim - 2) < step/2)

# --- Résumé ---

# Tableau 1 : Probabilités par intervalle
res <- data.frame(
  Intervalle_I = intervals[1:3],
  `P[S ∈ I] (Panjer)` = round(probs_Panjer[1:3], 6),
  `P[S ∈ I] (Empirique)` = round(probs_emp[1:3], 6),
  check.names = FALSE
)

cat("\nProbabilités par intervalles [k, k+1):\n")
print(res, row.names = FALSE, right = FALSE)

# Tableau 2 : Probabilités ponctuelles
res_point <- data.frame(
  Valeur = c(1, 2),
  `P[S = x] (Panjer discrétisé)` = round(c(P_S_eq1_Panjer, P_S_eq2_Panjer), 6),
  `P[S ≈ x] (Emp., ±step/2)` = round(c(P_S_eq1_emp, P_S_eq2_emp), 6),
  check.names = FALSE
)

cat("\nProbabilités ponctuelles approchées :\n")
print(res_point, row.names = FALSE, right = FALSE)




## === EXERCICE 2.3 ========================================

set.seed(123)

library(actuar)

set.seed(123)
Nsim <- 1e5
S_coll <- numeric(Nsim)

p_i <- c(0.80, 0.90, 0.95)

# --- Paramètres Poisson (approximation collective conservatrice) ---
lambda <- -log(p_i)
# lambda <- 1-p_i  # option plus précise pour la prime stop-loss (voir remarque slide 61 chapitre 5)
vals <- list(S1=c(1,2), S2=2, S3=1)
probs <- list(S1=c(0.15,0.05)/0.2, S2=1, S3=1)

# --- Simulation collective ---
for(i in 1:3){
  Ni <- rpois(Nsim, lambda[i])
  Xi <- numeric(Nsim)
  
  idx <- which(Ni >= 1)
  if(length(idx) > 0){
    if(length(vals[[i]]) == 1){
      Xi[idx] <- Ni[idx] * vals[[i]]
    } else {
      Xi[idx] <- sapply(Ni[idx], function(n)
        sum(sample(vals[[i]], n, replace=TRUE, prob=probs[[i]]))
      )
    }
  }
  S_coll <- S_coll + Xi
}

# --- Mesures de risque (approx. collective) ---
VaR_coll <- quantile(S_coll, 0.95)
E_stoploss_coll <- mean(pmax(S_coll - 2, 0))

cat("=== Approximation collective ===\n")
cat(sprintf("VaR[S; 0.95]    = %.3f\n", VaR_coll))
cat(sprintf("E[(S - 2)+]     = %.5f\n", E_stoploss_coll))

# --- Loi exacte (par convolution discrète) ---
pS1 <- c(0.80, 0.15, 0.05)
pS2 <- c(0.90, 0.00, 0.10)
pS3 <- c(0.95, 0.05, 0.00)
conv <- function(p, q) stats::convolve(p, rev(q), type = "open")
pS <- conv(conv(pS1, pS2), pS3)
pS <- pS / sum(pS)
valeurs_S <- 0:(length(pS) - 1)
cdf_S <- cumsum(pS)

VaR_exact <- valeurs_S[min(which(cdf_S >= 0.95))]
E_stoploss_exact <- sum(pmax(valeurs_S - 2, 0) * pS)

cat("\n=== Valeurs exactes ===\n")
cat(sprintf("VaR[S; 0.95]    = %d\n", VaR_exact))
cat(sprintf("E[(S - 2)+]     = %.5f\n", E_stoploss_exact))

# --- Graphique comparatif ---
plot(ecdf(S_coll),
     main = "CDF : approximation collective (simulation) vs exacte",
     xlab = "S", ylab = "F(S)", col = "red", lwd = 2,
     do.points = FALSE)
lines(valeurs_S, cdf_S, type = "s", col = "black", lwd = 2)
legend("bottomright",
       legend = c("CDF exacte (convolution)", "CDF collective (simulation)"),
       col = c("black", "red", "black", "red"),
       lty = c(1,1,2,2), lwd = 2)


## === EXERCICE 2.4 ========================================

# Paramètres
alpha <- 0.99       # niveau de risque
theta <- 1         # paramètre de l’exponentielle
n <- 1e6            # taille de l'échantillon simulé

# ---------------------------
# 1. Valeurs théoriques
# ---------------------------
VaR_th <- -log(1 - alpha) / theta
TVaR_th <- VaR_th + 1/theta   # pour l’exponentielle, TVaR = VaR + 1/theta

# ---------------------------
# 2. Simulation Monte-Carlo
# ---------------------------
set.seed(123)
X <- rexp(n, rate = theta)

# VaR empirique (quantile empirique)
VaR_emp <- quantile(X, probs = alpha)

# TVaR empirique = moyenne conditionnelle au-delà de la VaR
TVaR_emp <- mean(X[X > VaR_emp])

# ---------------------------
# 3. Comparaison
# ---------------------------
results <- data.frame(
  Theorique = c(VaR_th, TVaR_th),
  Empirique = c(VaR_emp, TVaR_emp),
  row.names = c("VaR", "TVaR")
)

print("Comparaison entre valeurs théoriques et empiriques :")
print(results)

# ---------------------------
# 4. Visualisation
# ---------------------------
hist(X, breaks = 100, probability = TRUE, col = "lightblue",
     main = "Loi exponentielle et mesures de risque",
     xlab = "Valeur")
curve(dexp(x, rate = theta), col = "red", lwd = 2, add = TRUE)

# Ajout de lignes pour la VaR et la TVaR
abline(v = VaR_th, col = "blue", lwd = 2, lty = 2)
abline(v = TVaR_th, col = "darkgreen", lwd = 2, lty = 2)

legend("topright",
       legend = c("Densité théorique", "VaR théorique", "TVaR théorique"),
       col = c("red", "blue", "darkgreen"), lty = 1:2, lwd = 2, bty = "n")


###################################

### Convergence VaR et TVaR pour une exponentielle ###

library(ggplot2)
library(dplyr)

# Paramètres
alpha <- 0.99
theta <- 1
N <- 1e5  # taille max de l'échantillon

# ---------------------------
# 1. Valeurs théoriques
# ---------------------------
VaR_th <- -log(1 - alpha) / theta
TVaR_th <- VaR_th + 1/theta

# ---------------------------
# 2. Simulation
# ---------------------------
set.seed(123)
X <- rexp(N, rate = theta)

# Tailles d’échantillons intermédiaires
n_seq <- unique(floor(exp(seq(log(100), log(N), length.out = 50))))

# Fonction pour calculer VaR et TVaR empiriques
compute_risks <- function(n) {
  x <- X[1:n]
  VaR_emp <- quantile(x, probs = alpha)
  TVaR_emp <- mean(x[x > VaR_emp])
  data.frame(n = n, VaR = as.numeric(VaR_emp), TVaR = TVaR_emp)
}

# Calcul pour toutes les tailles
res <- do.call(rbind, lapply(n_seq, compute_risks))

# Passage au format long pour ggplot
res_long <- res %>%
  tidyr::pivot_longer(cols = c("VaR", "TVaR"), names_to = "Mesure", values_to = "Valeur")

# ---------------------------
# 3. Visualisation ggplot
# ---------------------------
ggplot(res_long, aes(x = n, y = Valeur, color = Mesure)) +
  geom_line(size = 1) +
  geom_hline(yintercept = VaR_th, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = TVaR_th, color = "darkgreen", linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Convergence de la VaR et de la TVaR (Exponential)",
       x = "Taille de l'échantillon (log)", y = "Valeur",
       caption = paste("Paramètre theta =", theta, " - niveau alpha =", alpha)) +
  theme_minimal(base_size = 14)


## === EXERCICE 2.5 ========================================
library(ggplot2)
set.seed(123)

# ---- Simulation ----
n <- 1e6
X1 <- rnorm(n)
Y1 <- rnorm(n)   # indépendant
V  <- sample(c(-1,1), n, replace = TRUE)
X2 <- X1
Y2 <- V * X1

# ---- Sommes ----
S1 <- X1 + Y1
S2 <- X2 + Y2

# ---- Corrélations ----
rho1 <- cor(X1, Y1)
rho2 <- cor(X2, Y2)

cat("Corrélations :\n")
cat(sprintf("rho(X1, Y1) = %.3f\n", rho1))
cat(sprintf("rho(X2, Y2) = %.3f\n", rho2))

# ---- Niveaux de probabilité ----
p <- seq(0.8, 0.999, by = 0.001)

# ---- VaR théoriques ----
VaR1_theo <- sqrt(2) * qnorm(p)
VaR2_theo <- ifelse(p <= 0.5, 0, 2 * qnorm(2*p - 1))

# ---- Données pour ggplot ----
df_var <- data.frame(
  p = rep(p, 2),
  VaR = c(VaR1_theo, VaR2_theo),
  Portefeuille = factor(rep(c("X1+Y1", "X2+Y2"), each = length(p)))
)
# ---- Graphique ggplot ----
ggplot(df_var, aes(x = p, y = VaR, color = Portefeuille, linetype = Portefeuille)) +
  geom_line(linewidth = 1) +
  labs(title = "VaR des portefeuilles X1+Y1 et X2+Y2",
       x = "Niveau de probabilité p",
       y = "VaR(p)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(legend.position = "top")