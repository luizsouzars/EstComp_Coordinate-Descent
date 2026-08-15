# =============================================================================
#  03_referencia_1996.R
#
#  A rota de Tibshirani (1996) para o mesmo problema, implementada para que a
#  comparacao com 2010 seja MEDIDA e nao afirmada.
#
#  Na Secao 4 daquele artigo, o lasso e resolvido aproximando a penalidade
#  |b_j| por uma forma quadratica. Usando a desigualdade
#
#        |b| <= b^2 / (2|b0|) + |b0| / 2,          valida para todo b0 != 0,
#
#  com igualdade em b = b0, o problema penalizado vira, a cada iteracao, um
#  problema de ridge com pesos: (X'X + N lambda W) b = X'y, W = diag(1/|b0|).
#  E, em linguagem moderna, um algoritmo de majorizacao-minimizacao -- o mesmo
#  principio do IRLS. A ideia e boa; o problema esta em outro lugar.
#
#  TRES DEFEITOS ESTRUTURAIS, que a implementacao torna visiveis:
#
#   1. NENHUM ZERO EXATO nasce do algoritmo. Como 1/|b_j| explode quando b_j
#      se aproxima de zero, e preciso um limiar `eps` escolhido pelo usuario
#      para declarar um coeficiente nulo. Os zeros do modelo sao consequencia
#      de uma decisao arbitraria, nao da otimizacao.
#   2. O DESCARTE E IRREVERSIVEL. Zerada, uma variavel sai de W e nao pode
#      voltar. Num caminho de regularizacao isso e fatal: variaveis precisam
#      poder reentrar quando lambda diminui.
#   3. CUSTO O(|A|^3) por iteracao, por causa do sistema linear, contra O(pN)
#      de um ciclo completo de descida coordenada.
#
#  Compare com o que a reformulacao penalizada permite: a penalidade separavel
#  da forma fechada por coordenada, a forma fechada produz zeros exatos, e os
#  zeros exatos permitem tanto o conjunto ativo quanto a reentrada de
#  variaveis. O contraste nao e de eficiencia; e de estrutura.
# =============================================================================

## Resolve o lasso pela aproximacao de ridge iterado de Tibshirani (1996).
## Devolve tambem a trajetoria, para que se possa ver ONDE ela estabiliza.
ridge_iterado <- function(X, y, lambda, maxit = 50L, eps = 1e-4) {
  X <- as.matrix(X); y <- as.numeric(y)
  N <- nrow(X); p <- ncol(X)

  ## partida: ridge comum
  b <- as.vector(solve(crossprod(X) + N * lambda * diag(p), crossprod(X, y)))
  traj <- matrix(0, maxit, p)

  for (i in seq_len(maxit)) {
    b[abs(b) <= eps] <- 0                     # <- o limiar arbitrario
    A <- which(b != 0)
    if (!length(A)) { traj[i, ] <- b; next }

    Xa <- X[, A, drop = FALSE]
    W  <- diag(1 / abs(b[A]), nrow = length(A))
    ba <- solve(crossprod(Xa) + N * lambda * W, crossprod(Xa, y))

    b <- numeric(p)
    b[A] <- as.vector(ba)                     # variaveis fora de A nao voltam
    traj[i, ] <- b
  }
  list(beta = b, traj = traj, eps = eps, lambda = lambda)
}


## Objetivo do lasso na escala padronizada, para comparar rotas.
obj_lasso <- function(X, y, b, lambda) {
  sum((y - X %*% b)^2) / (2 * nrow(X)) + lambda * sum(abs(b))
}


## Trajetoria do objetivo sob descida coordenada, para o mesmo lambda e o mesmo
## numero de passos -- assim as duas rotas sao comparadas em pe de igualdade.
traco_cd <- function(X, y, lambda, maxit = 50L) {
  p <- ncol(X); N <- nrow(X)
  beta <- numeric(p); r <- y
  ob <- numeric(maxit)
  for (i in seq_len(maxit)) {
    s <- ciclo_gauss(X, r, beta, lambda, 1, rep(1, p), seq_len(p), N)
    r <- s$r; beta <- s$beta
    ob[i] <- obj_lasso(X, y, beta, lambda)
  }
  list(beta = beta, obj = ob)
}
