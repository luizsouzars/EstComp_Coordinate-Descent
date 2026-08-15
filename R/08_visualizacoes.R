# =============================================================================
#  08_visualizacoes.R -- a geometria do problema, em duas e tres dimensoes
#
#  Com p = 2 o problema inteiro cabe num grafico: o objetivo e uma superficie
#  sobre o plano (beta_1, beta_2), e o algoritmo e um caminho sobre ela. Tudo o
#  que este trabalho discute -- por que o soft-thresholding produz zeros, por que
#  a correlacao atrapalha, o que o warm start economiza, o que significa a
#  solucao nao ser unica -- aparece nessa figura.
#
#  Nada aqui e um caso especial inventado para a ilustracao: sao as MESMAS
#  funcoes de 01_caminho.R, rodando em p = 2.
#
#  Estas funcoes so CALCULAM. Quem desenha e o relatorio.
# =============================================================================


## Um problema de duas variaveis com correlacao controlada.
## rho alto deixa as curvas de nivel muito alongadas -- e e disso que a descida
## coordenada, que so anda paralelamente aos eixos, sofre.
dados_2d <- function(N = 100L, rho = 0.85, beta = c(1.4, 0.5), sigma = 1,
                     semente = 42) {
  set.seed(semente)
  z1 <- stats::rnorm(N)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(N)
  X  <- cbind(z1, z2)
  X  <- sweep(X, 2, colMeans(X), "-")
  X  <- sweep(X, 2, sqrt(colMeans(X^2)), "/")     # mesma padronizacao do en_path
  y  <- as.vector(X %*% beta) + stats::rnorm(N, sd = sigma)
  list(X = X, y = as.vector(y - mean(y)), beta = beta, rho = rho)
}


## Superficie do objetivo penalizado sobre uma grade de (beta_1, beta_2).
##
## O termo quadratico faz uma tigela elíptica; o termo l1 acrescenta um cone com
## QUINAS ao longo dos eixos. E a quina -- nao a inclinacao -- que produz zeros
## exatos: num ponto onde a funcao nao e diferenciavel, o subgradiente e um
## intervalo, e o minimo pode ficar preso ali para uma faixa inteira de valores
## de lambda.
superficie <- function(d, lambda, alpha = 1, lim = NULL, n = 90L) {
  N <- nrow(d$X)
  if (is.null(lim)) {
    bmqo <- as.vector(stats::coef(stats::lm(d$y ~ d$X + 0)))
    lim  <- c(min(-0.35, min(bmqo) - 0.55), max(0.35, max(bmqo) + 0.55))
  }
  b1 <- seq(lim[1], lim[2], length.out = n)
  b2 <- seq(lim[1], lim[2], length.out = n)
  pen <- function(b) alpha * abs(b) + (1 - alpha) * b^2 / 2
  Z <- outer(b1, b2, Vectorize(function(u, v) {
    r <- d$y - d$X %*% c(u, v)
    sum(r^2) / (2 * N) + lambda * (pen(u) + pen(v))
  }))
  list(b1 = b1, b2 = b2, Z = Z, lambda = lambda, alpha = alpha, lim = lim)
}


## DECOMPOSICAO DO OBJETIVO em suas duas partes, na mesma grade.
##
## Motivacao: plotar so "sem penalidade" ao lado de "com penalidade" nao funciona.
## Com lambda em faixa util, o termo l1 responde por uma fracao pequena da
## variacao total da superficie, e os dois paineis saem visualmente iguais. E
## preciso mostrar a PENALIDADE SOZINHA -- que e um cone com quinas sobre os
## eixos -- para que a origem dos vincos fique visivel.
decomposicao <- function(d, lambda, alpha = 1, lim = c(-1.2, 2.2), n = 45L) {
  N  <- nrow(d$X)
  b1 <- seq(lim[1], lim[2], length.out = n)
  b2 <- b1
  pen1 <- function(b) alpha * abs(b) + (1 - alpha) * b^2 / 2

  Zq <- outer(b1, b2, Vectorize(function(u, v)
    sum((d$y - d$X %*% c(u, v))^2) / (2 * N)))
  Zp <- outer(b1, b2, Vectorize(function(u, v)
    lambda * (pen1(u) + pen1(v))))

  mk <- function(Z) list(b1 = b1, b2 = b2, Z = Z, lambda = lambda,
                         alpha = alpha, lim = lim)
  list(quadratico = mk(Zq), penalidade = mk(Zp), soma = mk(Zq + Zp))
}


## CORTE UNIDIMENSIONAL: o objetivo como funcao de UM coeficiente, com o outro
## fixado no seu valor otimo. E a figura que prova o argumento do soft-threshold:
## a curva tem um BICO em zero, e conforme lambda cresce o minimo gruda no bico
## e fica la -- nao passa perto, fica exatamente em zero.
corte_1d <- function(d, lambdas, alpha = 1, j = 2L, lim = c(-1.0, 1.4),
                     n = 400L) {
  N <- nrow(d$X); outro <- if (j == 1L) 2L else 1L
  pen1 <- function(b) alpha * abs(b) + (1 - alpha) * b^2 / 2
  do.call(rbind, lapply(lambdas, function(lam) {
    bo <- en_path(d$X, d$y, alpha = alpha, lambda = lam,
                  padronizar = FALSE)$beta[, 1]
    grade <- seq(lim[1], lim[2], length.out = n)
    val <- vapply(grade, function(u) {
      b <- numeric(2); b[j] <- u; b[outro] <- bo[outro]
      sum((d$y - d$X %*% b)^2) / (2 * N) + lam * (pen1(b[1]) + pen1(b[2]))
    }, numeric(1))
    data.frame(beta = grade, objetivo = val, lambda = lam,
               otimo = bo[j], rotulo = sprintf("lambda = %.2f", lam))
  }))
}


## PERFIL AO LONGO DO VALE, para o caso de duas colunas identicas.
##
## Fixada a soma s = b1 + b2 no seu valor otimo, percorremos as divisoes
## possiveis: b1 = t*s, b2 = (1-t)*s, com t indo de 0 a 1. No lasso a curva e
## EXATAMENTE plana -- toda divisao e igualmente otima. Com alpha < 1 vira uma
## parabola com minimo em t = 1/2. Esta figura de uma linha so diz o que dois
## paineis 3D nao conseguem dizer.
perfil_vale <- function(N = 120L, lambda = 0.25, alphas = c(1, 0.9, 0.5),
                        semente = 3, n = 201L) {
  set.seed(semente)
  x <- stats::rnorm(N); x <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  X <- cbind(x, x); y <- as.vector(2 * x + stats::rnorm(N)); y <- y - mean(y)
  tt <- seq(0, 1, length.out = n)
  do.call(rbind, lapply(alphas, function(a) {
    b <- en_path(X, y, alpha = a, lambda = lambda, padronizar = FALSE)$beta[, 1]
    s_tot <- sum(b)
    val <- vapply(tt, function(t) {
      bb <- c(t * s_tot, (1 - t) * s_tot)
      sum((y - X %*% bb)^2) / (2 * N) +
        lambda * sum(a * abs(bb) + (1 - a) * bb^2 / 2)
    }, numeric(1))
    data.frame(t = tt, objetivo = val, alpha = a,
               rotulo = ifelse(a == 1, "lasso (alpha = 1)",
                               sprintf("elastic net (alpha = %.2f)", a)),
               relativo = val - min(val))
  }))
}


## Trajetoria da descida coordenada, um COEFICIENTE por passo.
##
## O caminho e sempre paralelo a um dos eixos: essa e a assinatura visual do
## metodo. Cada passo resolve exatamente o problema unidimensional na coordenada
## corrente -- vai ate o fundo do vale NAQUELA direcao -- e depois troca de
## direcao. Por isso a trajetoria "escorrega" pela parede do vale quando as
## curvas de nivel sao alongadas.
traj_2d <- function(d, lambda, alpha = 1, inicio = c(0, 0), n_passos = 24L) {
  N <- nrow(d$X); b <- inicio; r <- d$y - d$X %*% b
  pts <- matrix(NA_real_, n_passos + 1L, 2); pts[1, ] <- b
  for (t in seq_len(n_passos)) {
    j   <- 1L + (t - 1L) %% 2L                       # varredura ciclica: 1, 2, 1, 2...
    rho <- sum(d$X[, j] * r) / N + b[j]
    bn  <- S(rho, lambda * alpha) / (1 + lambda * (1 - alpha))
    r   <- r - d$X[, j] * (bn - b[j])
    b[j] <- bn
    pts[t + 1L, ] <- b
  }
  obj <- apply(pts, 1, function(bb) {
    rr <- d$y - d$X %*% bb
    sum(rr^2) / (2 * N) +
      lambda * sum(alpha * abs(bb) + (1 - alpha) * bb^2 / 2)
  })
  data.frame(passo = 0:n_passos, b1 = pts[, 1], b2 = pts[, 2], objetivo = obj)
}


## Quantos passos de coordenada ate convergir, partindo de um ponto dado.
## E a medida usada para comparar warm start com partida fria em p = 2.
passos_ate_convergir <- function(d, lambda, alpha = 1, inicio = c(0, 0),
                                 tol = 1e-10, max_passos = 2000L) {
  N <- nrow(d$X); b <- inicio; r <- d$y - d$X %*% b
  for (t in seq_len(max_passos)) {
    j   <- 1L + (t - 1L) %% 2L
    rho <- sum(d$X[, j] * r) / N + b[j]
    bn  <- S(rho, lambda * alpha) / (1 + lambda * (1 - alpha))
    dif <- bn - b[j]
    r   <- r - d$X[, j] * dif; b[j] <- bn
    if (t %% 2L == 0L && abs(dif) < sqrt(tol)) return(list(passos = t, beta = b))
  }
  list(passos = max_passos, beta = b)
}


## WARM START contra PARTIDA FRIA, sobre a grade de lambda.
##
## Para cada lambda da grade, conta quantos passos de coordenada sao necessarios
## partindo (a) da solucao do lambda anterior e (b) da origem. A diferenca e,
## literalmente, o que o caminho de regularizacao economiza.
comparacao_warm <- function(d, nlambda = 18L, alpha = 1, lambda_min_ratio = 0.02) {
  N <- nrow(d$X)
  lmax <- max(abs(crossprod(d$X, d$y))) / (N * max(alpha, 1e-3))
  lam  <- exp(seq(log(lmax * 0.98), log(lmax * lambda_min_ratio),
                  length.out = nlambda))
  sol  <- matrix(NA_real_, nlambda, 2)
  quente <- fria <- integer(nlambda)
  b_ant <- c(0, 0)
  for (k in seq_len(nlambda)) {
    q <- passos_ate_convergir(d, lam[k], alpha, inicio = b_ant)
    f <- passos_ate_convergir(d, lam[k], alpha, inicio = c(0, 0))
    quente[k] <- q$passos; fria[k] <- f$passos
    sol[k, ]  <- q$beta;   b_ant   <- q$beta
  }
  list(lambda = lam,
       solucoes = data.frame(lambda = lam, b1 = sol[, 1], b2 = sol[, 2]),
       passos = rbind(
         data.frame(lambda = lam, passos = quente, partida = "warm start (solucao anterior)"),
         data.frame(lambda = lam, passos = fria,   partida = "partida fria (do zero)")))
}


## O VALE DE NAO UNICIDADE.
##
## Com duas colunas IDENTICAS, qualquer par (b1, b2) com b1 + b2 = c e mesmo
## sinal produz o mesmo ajuste e a mesma norma l1 -- logo o mesmo objetivo. O
## conjunto de otimos e um SEGMENTO DE RETA, nao um ponto, e a superficie tem um
## vale de fundo perfeitamente plano. Com alpha < 1 o termo quadratico inclina o
## fundo do vale e cria um minimo unico, no meio, por simetria.
vale_nao_unicidade <- function(N = 120L, lambda = 0.25, semente = 3, n = 90L) {
  set.seed(semente)
  x <- stats::rnorm(N); x <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  X <- cbind(x, x)                                  # colunas exatamente iguais
  y <- as.vector(2 * x + stats::rnorm(N)); y <- y - mean(y)
  d <- list(X = X, y = y)
  list(lasso = superficie(d, lambda, alpha = 1,   lim = c(-0.6, 2.6), n = n),
       enet  = superficie(d, lambda, alpha = 0.5, lim = c(-0.6, 2.6), n = n),
       dados = d)
}


## Caminho de coeficientes em formato longo, para o grafico que da nome ao metodo.
caminho_para_df <- function(fit, verdadeiro = NULL) {
  p <- nrow(fit$beta); nl <- length(fit$lambda)
  out <- data.frame(
    j      = rep(seq_len(p), times = nl),
    lambda = rep(fit$lambda, each = p),
    coef   = as.vector(fit$beta))
  out$grupo <- if (is.null(verdadeiro)) "todas"
               else rep(ifelse(verdadeiro != 0, "efeito verdadeiro", "ruido"), times = nl)
  out
}
