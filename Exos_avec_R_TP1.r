# TP1

# Exercice 1.1

# v1 avec segments verticaux et horizontaux pour la CDF

# Valeurs et probabilités
N_values <- 1:4
N_probs  <- c(0.2, 0.15, 0.25, 0.4)
F_cum    <- cumsum(N_probs)

# Fonction CDF pour tout x
F_N_func <- function(x) {
  sapply(x, function(z) {
    if (z < N_values[1]) return(0)
    else if (z >= N_values[length(N_values)]) return(1)
    else return(F_cum[max(which(N_values <= z))])
  })
}

# Exemple d'utilisation
x_test <- seq(0, 5, 0.1)
F_vals <- F_N_func(x_test)

# Visualisation avec ggplot2
library(ggplot2)
df <- data.frame(x = x_test, F = F_vals)

ggplot(df, aes(x = x, y = F)) +
  geom_step(direction = "hv", color = "blue", size = 1.2) +
  #geom_point(data = data.frame(N=N_values, F=F_cum), aes(x=N, y=F), color="red", size=3) +
  labs(title = expression("Fonction de répartition"~F[N](x)~"pour tout réel x"),
       x = expression(x), y = expression(F[N](x))) +
  theme_minimal(base_size = 14)

# v2 avec segments horizontaux uniquement

# Valeurs et probabilités
N_values <- 1:4
N_probs  <- c(0.2, 0.15, 0.25, 0.4)
F_cum    <- cumsum(N_probs)

# Points exacts pour N
points_df <- data.frame(N = N_values, F = F_cum)

# Segments horizontaux incluant avant 0 et après 4
segments_df <- data.frame(
  x_start = c(-1, N_values[-length(N_values)], N_values[length(N_values)]),
  x_end   = c(N_values[1], N_values[-1], 5),
  y       = c(0, F_cum[-length(F_cum)], F_cum[length(F_cum)])
)

library(ggplot2)

ggplot() +
  # segments horizontaux
  geom_segment(data = segments_df,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "blue", size = 1.2) +
  # points exacts
  geom_point(data = points_df, aes(x = N, y = F), color="red", size=3) +
  labs(title = expression("Fonction de répartition"~F[N](x)),
       x = expression(x), y = expression(F[N](x))) +
  theme_minimal(base_size = 14)


# Fonction quantile de N
qN <- function(p) {
  sapply(p, function(z) {
    N_values[min(which(F_cum >= z))]
  })
}

# Grille de probabilités
p_vals <- seq(0, 1, 0.001)
q_vals <- qN(p_vals)

# Dataframe pour ggplot
df <- data.frame(p = p_vals, Q = q_vals)

# Points exacts pour les discontinuités
points_df <- data.frame(
  p = F_cum,
  Q = N_values
)

library(ggplot2)

ggplot(df, aes(x = p, y = Q)) +
  geom_step(direction = "hv", color = "blue", size = 1.2) +  # escalier horizontal
  geom_point(data = points_df, aes(x = p, y = Q), color = "red", size = 3) +
  labs(title = "Fonction quantile de N",
       x = "Probabilité p", y = "Quantile Q(p)") +
  theme_minimal(base_size = 14)

# Simulation de la loi de N

# Fonction pour simuler N
simuler_N <- function(n) {
  sample(N_values, size = n, replace = TRUE, prob = N_probs)
}

# Simulation de 100 tirages
set.seed(123)
n <- 100
tirages <- simuler_N(n)

# Préparer le dataframe pour ggplot2
df <- data.frame(N = tirages)

library(ggplot2)

# Histogramme / distribution empirique
ggplot(df, aes(x = factor(N), group = 1)) +  # <- ajout de group = 1
  geom_bar(aes(y = after_stat(prop)), fill = "skyblue", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = paste("Distribution empirique de N sur", n, "tirages"),
       x = "Valeurs de N", y = "Fréquence (%)") +
  theme_minimal(base_size = 14)

################################
# Exercice 1.2

# Packages
library(ggplot2)

# Paramètres
alpha <- 20
lambda <- 100

# Définition de la densité
f <- function(x, alpha, lambda) {
  (alpha * lambda^alpha) / ( (x + lambda)^(alpha + 1) )
}

# Définition de la fonction de répartition (intégration numérique)
F_num <- function(x, alpha, lambda) {
  sapply(x, function(t) {
    integrate(f, lower = 0, upper = t, alpha = alpha, lambda = lambda)$value
  })
}

# Grille de points
x_vals <- seq(0, 40, length.out = 200)

# Calcul de la CDF numérique de X
cdf_num  <- F_num(x_vals, alpha, lambda)

# Assemblage en data frame
df <- data.frame(x = x_vals, F = cdf_num)

# Visualisation avec ggplot
ggplot(df, aes(x = x, y = F)) +
  geom_line(color = "blue", linewidth = 1.2) +
  labs(title = "Fonction de répartition (CDF)",
       subtitle = paste("alpha =", alpha, ", lambda =", lambda),
       x = "x", y = "F(x)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



# Fonction quantile de la loi X
q_X <- function(p, alpha, lambda) {
  lambda * ( (1 / (1 - p)^(1/alpha)) - 1 )
}

# Exemple : quelques quantiles
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
quantiles <- q_X(probs, alpha, lambda)

# Résultats
data.frame(Probabilite = probs, Quantile = quantiles)

# Simulation de la loi X

n <- 1000

# Fonction de simulation (inversion de la CDF)
r_X <- function(n, alpha, lambda) {
  u <- runif(n)
  lambda * ((1 / (1 - u)^(1/alpha)) - 1)
}

# Fonction CDF théorique
p_X <- function(x, alpha, lambda) {
  1 - (lambda / (x + lambda))^alpha
}

# Simulation
set.seed(123)
x <- r_X(n, alpha, lambda)
df <- data.frame(x = x)

# Tracé de la CDF empirique + théorique
ggplot(df, aes(x = x)) +
  stat_ecdf(geom = "step", color = "blue", linewidth = 1.2) +  # CDF empirique
  stat_function(fun = function(t) p_X(t, alpha, lambda),
                color = "red", size = 1.2, linetype = "dashed") +  # CDF théorique
  labs(title = "CDF empirique vs CDF théorique de la loi X",
       subtitle = paste("alpha =", alpha, ", lambda =", lambda, ", n =", n),
       x = "x", y = "F(x)",
       caption = "Bleu : CDF empirique | Rouge pointillé : CDF théorique") +
  theme_minimal(base_size = 14)


# Moyenne de l'échantillon
moyenne <- mean(x)

# Écart-type de l'échantillon
ecart_type <- sd(x)

# Affichage
cat("Moyenne de l'échantillon :", moyenne, "\n")
cat("Écart-type de l'échantillon :", ecart_type, "\n")

# package actuar
library(actuar)
help(rpareto)

################################

# Exercice 1.3

library(ggplot2)

set.seed(123)
lambda <- 5     # paramètre Poisson
q_ber  <- 0.3   # paramètre Bernoulli
n_sims <- 1000

# Simulation de S
S <- numeric(n_sims)
for (i in 1:n_sims) {
  N <- rpois(1, lambda)          # tirage de N
  if (N == 0) {
    S[i] <- 0
  } else {
    X <- rbinom(N, size = 1, prob = q_ber)  # N Bernoulli
    S[i] <- sum(X)
  }
}

# Dataframe pour ggplot2
df <- data.frame(S = S)

# Distribution empirique et théorique
ggplot(df, aes(x = S)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = function(x) dpois(round(x), lambda * q_ber), color = "red", size = 1.2) +
  scale_x_continuous(breaks = 0:max(S)) +
  labs(title = paste("Distribution de S ~ Poisson composée (N~Poisson, X~Bernoulli)"),
       x = "S", y = "Densité",
       caption = "Histogramme : simulation | Rouge : loi théorique Poisson(lambda*q_ber)") +
  theme_minimal(base_size = 14)

# CDF empirique vs CDF théorique

# Valeurs entières pour la CDF théorique
S_max <- max(S)
x_vals <- 0:S_max
cdf_theo <- ppois(x_vals, lambda * q_ber)
df_theo <- data.frame(S = x_vals, F = cdf_theo)

# Tracé
ggplot() +
  stat_ecdf(data = df, aes(x = S), geom = "step", color = "blue", size = 1.2) +
  geom_step(data = df_theo, aes(x = S, y = F), color = "red", size = 1.2, linetype = "dashed") +
  scale_x_continuous(breaks = 0:S_max) +
  labs(title = "CDF empirique vs CDF théorique de S",
       x = "S", y = "F(S)",
       caption = "Bleu : empirique | Rouge pointillé : théorique Poisson(lambda*theta)") +
  theme_minimal(base_size = 14)

# Calcul de certaines probabilités
mu <- lambda * q_ber  # paramètre Poisson

# Probabilités exactes
dpois(0:2, mu)
# Résultat attendu :
# [1] 0.2231302 0.3346953 0.2510215

###################################

# Exercice 1.4

library(ggplot2)

# 1. Paramètres
n <- 10
p <- 0.1
alpha <- 1/2
tau <- 1/3
n_sims <- 10000

# 2. Simulation directe
set.seed(123)
S_sim <- numeric(n_sims)

for (i in 1:n_sims) {
  N <- rbinom(1, size = n, prob = p)
  if (N == 0) {
    S_sim[i] <- 0
  } else {
    C <- rgamma(N, shape = alpha, rate = tau)
    S_sim[i] <- sum(C)
  }
}

df_sim <- data.frame(S = S_sim)

# CDF empirique
df_sim_cdf <- data.frame(
  S = sort(S_sim),
  CDF = ecdf(S_sim)(sort(S_sim))
)

ggplot() +
  geom_step(data = df_sim_cdf, aes(x = S, y = CDF),
            color = "blue", size = 1, direction = "hv") +
  labs(title = "CDF empirique de S (simulation)",
       x = "S", y = "CDF") +
  theme_minimal(base_size = 14)


# Calcul de certaines probabilités
F_S <- function(s) {
  sapply(s, function(x) {
    if (x < 0) return(0)
    k_vals <- 1:n
    pk <- dbinom(k_vals, size = n, prob = p)  # P(N=k)
    terms <- sum(pk * pgamma(x, shape = k_vals * alpha, rate = tau))
    
    # Ajouter la masse en N=0
    dbinom(0, size = n, prob = p) + terms
  })
}

# exemples
F_S(c(0, 1, 2))