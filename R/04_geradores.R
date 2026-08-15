# =============================================================================
#  04_geradores.R -- cenarios de simulacao e metricas de avaliacao
# =============================================================================

## Preditores gaussianos com correlacao AR(1): cor(x_i, x_j) = rho^|i-j|.
## E o cenario do Exemplo 1 de Tibshirani (1996), o que permite comparar com
## os resultados originais.
gera_ar1 <- function(N, p, rho, beta, sigma = 1, n_teste = 0) {
  Sig <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  R   <- chol(Sig)
  X   <- matrix(stats::rnorm(N * p), N, p) %*% R
  y   <- as.vector(X %*% beta) + stats::rnorm(N, sd = sigma)
  out <- list(X = X, y = y, beta = beta, Sigma = Sig, sigma = sigma)
  if (n_teste > 0) {
    Xt <- matrix(stats::rnorm(n_teste * p), n_teste, p) %*% R
    out$X_teste  <- Xt
    out$mu_teste <- as.vector(Xt %*% beta)            # media verdadeira
    out$y_teste  <- out$mu_teste + stats::rnorm(n_teste, sd = sigma)
  }
  out
}


## Preditores equicorrelacionados, delineamento de tempos da Secao 5.1 de
## Friedman, Hastie & Tibshirani (2010): x_j = sqrt(rho) z + sqrt(1-rho) e_j.
gera_equicor <- function(N, p, rho, snr = 3) {
  z <- stats::rnorm(N)
  X <- sqrt(rho) * matrix(z, N, p) +
       sqrt(1 - rho) * matrix(stats::rnorm(N * p), N, p)
  beta <- (-1)^seq_len(p) * exp(-2 * (seq_len(p) - 1) / 20)
  eta  <- as.vector(X %*% beta)
  y    <- eta + stats::rnorm(N, sd = stats::sd(eta) / snr)
  list(X = X, y = y, beta = beta)
}


## g grupos de `tam` preditores quase identicos (correlacao rho_g dentro do
## grupo, zero entre grupos) mais `n_ruido` preditores independentes. Somente o
## PRIMEIRO grupo tem efeito.
##
## Este e o cenario critico para a contribuicao original: dentro de um grupo
## quase identico a solucao do lasso pode nao ser unica, e o desempate cabe ao
## algoritmo -- nao aos dados.
gera_grupos <- function(N, g = 3L, tam = 5L, n_ruido = 15L, rho_g = 0.99,
                        beta_grupo = 1, sigma = 1) {
  p <- g * tam + n_ruido
  Z <- matrix(stats::rnorm(N * g), N, g)
  X <- matrix(0, N, p)
  for (j in seq_len(g)) {
    X[, ((j - 1) * tam + 1):(j * tam)] <-
      sqrt(rho_g) * Z[, j] +
      sqrt(1 - rho_g) * matrix(stats::rnorm(N * tam), N, tam)
  }
  if (n_ruido > 0)
    X[, (g * tam + 1):p] <- matrix(stats::rnorm(N * n_ruido), N, n_ruido)
  beta <- numeric(p); beta[seq_len(tam)] <- beta_grupo
  list(X = X, y = as.vector(X %*% beta) + stats::rnorm(N, sd = sigma),
       beta = beta, idx_efeito = seq_len(tam), p = p, tam = tam, g = g)
}


# ---------------------------------------------------------------------------
# Metricas
# ---------------------------------------------------------------------------

metricas_selecao <- function(beta_est, beta_verdadeiro) {
  sel  <- beta_est != 0
  verd <- beta_verdadeiro != 0
  VP <- sum(sel & verd); FP <- sum(sel & !verd)
  FN <- sum(!sel & verd); VN <- sum(!sel & !verd)
  tpr  <- if ((VP + FN) > 0) VP / (VP + FN) else NA_real_
  fpr  <- if ((FP + VN) > 0) FP / (FP + VN) else NA_real_
  prec <- if ((VP + FP) > 0) VP / (VP + FP) else NA_real_
  f1   <- if (is.finite(prec) && is.finite(tpr) && (prec + tpr) > 0)
            2 * prec * tpr / (prec + tpr) else 0
  c(TPR = tpr, FPR = fpr, F1 = f1, tamanho = sum(sel))
}

## Indice de Jaccard entre dois conjuntos de indices selecionados.
jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

## Jaccard medio sobre todos os pares de uma lista de conjuntos.
## Vale 1 quando o mesmo conjunto e sempre escolhido.
estabilidade <- function(lista) {
  B <- length(lista)
  if (B < 2) return(NA_real_)
  pares <- utils::combn(B, 2)
  mean(apply(pares, 2, function(ij) jaccard(lista[[ij[1]]], lista[[ij[2]]])))
}

## Media e erro-padrao de Monte Carlo por celula do delineamento.
## Toda media reportada no trabalho vem acompanhada do seu EP: sem ele nao se
## distingue efeito de ruido de simulacao.
resume_mc <- function(dados, resposta, por) {
  f <- stats::as.formula(paste(resposta, "~", paste(por, collapse = " + ")))
  m <- stats::aggregate(f, data = dados, FUN = mean)
  s <- stats::aggregate(f, data = dados,
                        FUN = function(z) stats::sd(z) / sqrt(length(z)))
  n <- stats::aggregate(f, data = dados, FUN = length)
  names(m)[ncol(m)] <- "media"
  m$ep <- s[[ncol(s)]]
  m$R  <- n[[ncol(n)]]
  m
}
