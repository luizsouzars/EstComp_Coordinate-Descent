# =============================================================================
#  01_caminho.R -- A IMPLEMENTACAO PRINCIPAL DO TRABALHO
#
#  Friedman, Hastie & Tibshirani (2010), "Regularization Paths for Generalized
#  Linear Models via Coordinate Descent", JSS 33(1), 1-22.
#
#  ESCOPO. Implementamos UM objeto: o caminho de regularizacao para o caso
#  gaussiano, percorrido por descida coordenada ciclica com INICIALIZACAO
#  QUENTE (warm start). Nao implementamos GLMs, atualizacoes por covariancia,
#  multinomial nem resposta agrupada -- essas sao extensoes construidas SOBRE
#  este objeto, e o recorte esta justificado no relatorio.
#
#  Problema resolvido, para cada lambda de uma grade decrescente:
#
#     min_{b0, b}  (1/2N) sum_i (y_i - b0 - x_i'b)^2 + lambda * P_alpha(b)
#     P_alpha(b) = sum_j [ (1-alpha) b_j^2 / 2 + alpha |b_j| ] * pf_j
#
#  DEPENDENCIAS: nenhuma. So 'base' e 'stats'.
#
#  TRES INTERRUPTORES DE ABLACAO, usados no estudo de Monte Carlo. Eles NAO
#  mudam a solucao, so o caminho percorrido ate ela:
#     warm  = TRUE/FALSE   inicializacao quente ligada/desligada
#     ativo = TRUE/FALSE   estrategia de conjunto ativo ligada/desligada
#     ordem = "ciclica" | "aleatoria" | vetor de indices
#
#  CONTABILIDADE DE CUSTO. Alem do tempo de relogio, contamos VISITAS: o numero
#  de vezes que uma coordenada foi avaliada. Cada visita custa O(N) e essa
#  contagem nao depende da maquina, do sistema operacional nem da carga -- e a
#  unidade certa para comparar algoritmos.
# =============================================================================


## Operador de soft-thresholding, S(z, g) = sinal(z) (|z| - g)_+.
## E a solucao em forma fechada do lasso univariado; aparece na equacao (3) de
## Tibshirani (1996) como um resultado analitico isolado para o caso
## ortonormal, e em 2010 vira o passo elementar executado milhoes de vezes.
S <- function(z, gamma) {
  sign(z) * pmax(abs(z) - gamma, 0)
}


## lambda_max: o menor lambda que zera TODOS os coeficientes penalizados.
## Acima dele o soft-thresholding zera toda coordenada, e beta = 0 e a solucao
## exata -- e e por isso que o caminho pode partir de um ponto conhecido.
## A versao com sufixo _pad recebe os dados ja centrados e padronizados; existe
## para que en_path nao recalcule a padronizacao e para que nenhuma outra funcao
## precise duplicar esta conta.
.lambda_max_pad <- function(Xs, yc, alpha, pf) {
  a_ef <- max(alpha, 1e-3)                     # evita divisao por zero no ridge
  g    <- abs(as.vector(crossprod(Xs, yc))) / nrow(Xs)
  pen  <- pf > 0
  max(if (any(pen)) max(g[pen] / (a_ef * pf[pen])) else 1, .Machine$double.eps)
}

lambda_max <- function(X, y, alpha = 1, pf = rep(1, ncol(X)), padronizar = TRUE) {
  X   <- as.matrix(X)
  mu  <- colMeans(X); Xc <- sweep(X, 2, mu, "-")
  sdn <- if (padronizar) sqrt(colMeans(Xc^2)) else rep(1, ncol(X))
  sdn[sdn == 0 | !is.finite(sdn)] <- 1
  .lambda_max_pad(sweep(Xc, 2, sdn, "/"), y - mean(y), alpha, pf)
}


## Um ciclo de descida coordenada sobre os indices em `idx`.
## Para cada j: coeficiente de MQ simples do residuo parcial, depois
## soft-thresholding, depois encolhimento proporcional (a parte l2).
## O residuo r e atualizado em O(N) apenas quando o coeficiente muda -- e daqui
## que vem o ganho com solucoes esparsas.
ciclo_gauss <- function(X, r, beta, lambda, alpha, pf, idx, N) {
  dmax <- 0
  for (j in idx) {
    bj  <- beta[j]
    rho <- sum(X[, j] * r) / N + bj          # <x_j, r>/N + beta_j
    den <- 1 + lambda * (1 - alpha) * pf[j]
    bn  <- S(rho, lambda * alpha * pf[j]) / den
    d   <- bn - bj
    if (d != 0) {
      r       <- r - X[, j] * d
      beta[j] <- bn
      if (d * d > dmax) dmax <- d * d
    }
  }
  list(r = r, beta = beta, dmax = dmax, visitas = length(idx))
}


## Itera ciclos ate a convergencia, para UM valor de lambda.
##
## Com `ativo = TRUE` alterna-se entre dois regimes:
##   (a) um ciclo completo sobre as p coordenadas, que pode admitir variaveis
##       novas no modelo;
##   (b) ciclos restritos ao conjunto ativo A = {j : beta_j != 0}, muito mais
##       baratos (custo |A| em vez de p por ciclo).
## A convergencia so e declarada apos um ciclo COMPLETO que nao altera nada --
## e essa exigencia que garante que nenhuma variavel de fora deveria ter
## entrado. Com `ativo = FALSE`, o regime (b) e desligado: e a descida
## coordenada ciclica pura.
## `tol2` e o limiar aplicado a max_j (delta beta_j)^2, ja elevado ao quadrado
## por quem chama (en_path), que decide se ele escala com lambda ou nao.
## `paciencia`: numero de ciclos sem melhora antes de desistir. Sem esse freio,
## um lambda patologico -- tipicamente no fim do caminho, com colunas
## exatamente dependentes, onde coeficientes minusculos oscilam em torno de zero
## sem nunca satisfazer o criterio -- consome `maxit` ciclos inteiros. Numa
## simulacao com centenas de ajustes isso e a diferenca entre minutos e horas.
## Desistir cedo e registrar `convergiu = FALSE` e mais honesto e muito mais
## barato do que insistir.
cd_ate_convergir <- function(X, r, beta, lambda, alpha, pf, tol2, maxit,
                             ativo, ordem, N, paciencia = 50L) {
  p <- ncol(X)
  ciclos <- 0L; visitas <- 0; convergiu <- FALSE
  melhor <- Inf; sem_melhora <- 0L

  ordem_de <- function(ids) {
    if (identical(ordem, "ciclica")) ids
    ## Atencao: sample(ids) seria um bug. Quando ids tem UM elemento, digamos 7,
    ## sample(7) devolve uma permutacao de 1:7 em vez do proprio 7 -- e o ciclo
    ## visitaria coordenadas erradas. Isso acontece sempre que o conjunto ativo
    ## tem tamanho 1, o que e comum no inicio do caminho.
    else if (identical(ordem, "aleatoria")) ids[sample.int(length(ids))]
    else {                                   # ordem fixa fornecida pelo usuario
      o <- as.integer(ordem); o[o %in% ids]
    }
  }

  repeat {
    s <- ciclo_gauss(X, r, beta, lambda, alpha, pf, ordem_de(seq_len(p)), N)
    r <- s$r; beta <- s$beta
    ciclos <- ciclos + 1L; visitas <- visitas + s$visitas

    if (s$dmax < tol2) { convergiu <- TRUE; break }  # ciclo completo sem mudanca
    if (s$dmax < melhor * 0.99) { melhor <- s$dmax; sem_melhora <- 0L }
    else sem_melhora <- sem_melhora + 1L
    if (sem_melhora >= paciencia || ciclos >= maxit) break

    if (ativo) {
      A <- which(beta != 0)
      while (length(A) && ciclos < maxit) {
        sa <- ciclo_gauss(X, r, beta, lambda, alpha, pf, ordem_de(A), N)
        r <- sa$r; beta <- sa$beta
        ciclos <- ciclos + 1L; visitas <- visitas + sa$visitas
        if (sa$dmax < tol2) break
      }
    }
  }
  list(r = r, beta = beta, ciclos = ciclos, visitas = visitas, convergiu = convergiu)
}


## O CAMINHO DE REGULARIZACAO -- a contribuicao central do artigo de 2010.
##
## A grade comeca em lambda_max, o menor lambda que zera todos os coeficientes.
## Como em lambda_max a solucao e conhecida em forma fechada (beta = 0), o
## caminho parte de um ponto exato e cada solucao serve de ponto de partida
## para o proximo lambda. Nada e resolvido "do zero" a nao ser a primeira.
## NOTA SOBRE O CRITERIO DE PARADA -- e a decisao mais delicada da implementacao.
##
## Paramos quando um ciclo completo muda todo coeficiente por menos que `tol`.
## Duas armadilhas, ambas encontradas na pratica ao rodar este trabalho:
##
## (1) O criterio controla o PASSO, nao a DISTANCIA AO OTIMO. Sob correlacao
##     alta a descida coordenada converge linearmente com razao proxima de 1, e
##     a distancia ao otimo e da ordem do passo dividido por (1 - razao). Passo
##     pequeno nao e sinonimo de estar perto.
##
## (2) Pior: um limiar ABSOLUTO nao e equivariante em escala. As condicoes KKT
##     comparam o gradiente com lambda*alpha, e lambda varia quatro ordens de
##     grandeza ao longo do caminho. Um erro de 1e-5 e desprezivel quando
##     lambda*alpha = 1 e catastrofico quando lambda*alpha = 1e-5. Na primeira
##     execucao deste trabalho, com limiar absoluto, a violacao KKT RELATIVA
##     chegou a 4,6 no fim do caminho -- ou seja, o gradiente estava mais longe
##     da condicao de otimalidade do que o proprio lambda.
##
## Por isso o padrao e `escala_tol = "lambda"`: exige-se |delta beta_j| < tol *
## lambda_k. O criterio fica proporcional a escala local do problema e a
## violacao KKT relativa passa a ser aproximadamente constante ao longo do
## caminho. `escala_tol = "absoluta"` reproduz o comportamento usual, e serve
## para exibir o contraste.
##
## CUSTO. Escalar por lambda encarece o FIM do caminho, onde o limiar efetivo
## vira tol * lambda_min_ratio * lambda_max. Com tol = 1e-5 e razao 1e-4, exige-se
## |delta beta| < 1e-9 nos ultimos valores de lambda -- e ali o problema esta
## proximo de minimos quadrados sem penalizacao, regime em que a descida
## coordenada converge devagar se as colunas forem quase colineares. Na pratica
## isso multiplica o tempo do caminho por algo entre 1,5 e 3. Se um ajuste
## demorar demais, as saidas sao, nesta ordem: aumentar lambda_min_ratio (o fim
## do caminho raramente interessa), afrouxar tol, ou aceitar `convergiu = FALSE`
## em alguns lambdas -- que o objeto devolvido registra honestamente.
en_path <- function(X, y, alpha = 1, lambda = NULL, nlambda = 100L,
                    lambda_min_ratio = if (nrow(X) < ncol(X)) 0.01 else 1e-4,
                    pf = rep(1, ncol(X)), padronizar = TRUE,
                    warm = TRUE, ativo = TRUE, ordem = "ciclica",
                    tol = 1e-5, escala_tol = c("lambda", "absoluta"),
                    maxit = 2000L) {

  escala_tol <- match.arg(escala_tol)

  X <- as.matrix(X); y <- as.numeric(y)
  N <- nrow(X); p <- ncol(X)
  stopifnot(length(y) == N, alpha >= 0, alpha <= 1, length(pf) == p, all(pf >= 0))

  ## --- padronizacao: sem ela, lambda significaria coisas diferentes para
  ## --- variaveis em escalas diferentes
  mu  <- colMeans(X)
  Xc  <- sweep(X, 2, mu, "-")
  sdn <- if (padronizar) sqrt(colMeans(Xc^2)) else rep(1, p)
  sdn[sdn == 0 | !is.finite(sdn)] <- 1
  Xs  <- sweep(Xc, 2, sdn, "/")
  ybar <- mean(y); yc <- y - ybar
  ## o intercepto nao e penalizado e, com tudo centrado, sai de graca no fim

  ## --- grade de lambda ------------------------------------------------------
  if (is.null(lambda)) {
    lmax   <- .lambda_max_pad(Xs, yc, alpha, pf)
    lambda <- exp(seq(log(lmax), log(lmax * lambda_min_ratio), length.out = nlambda))
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  nl <- length(lambda)

  ## --- percurso -------------------------------------------------------------
  B <- matrix(0, p, nl)
  ciclos <- integer(nl); visitas <- numeric(nl); conv <- logical(nl)
  beta <- numeric(p); r <- yc

  for (k in seq_len(nl)) {
    if (!warm) { beta <- numeric(p); r <- yc }   # ABLACAO: partida fria
    ## o criterio compara (delta beta)^2, entao o limiar entra ao quadrado
    tol2 <- if (escala_tol == "lambda") (tol * lambda[k])^2 else tol^2
    f <- cd_ate_convergir(Xs, r, beta, lambda[k], alpha, pf,
                          tol2, maxit, ativo, ordem, N)
    beta <- f$beta; r <- f$r
    B[, k] <- beta
    ciclos[k] <- f$ciclos; visitas[k] <- f$visitas; conv[k] <- f$convergiu
  }

  ## --- volta a escala original ---------------------------------------------
  Bo <- B / sdn
  a0 <- as.vector(ybar - crossprod(Bo, mu))

  structure(list(
    lambda = lambda, beta = Bo, a0 = a0, df = colSums(B != 0),
    alpha = alpha, pf = pf, beta_pad = B, escala = sdn, centro = mu,
    ciclos = ciclos, visitas = visitas, convergiu = conv,
    warm = warm, ativo = ativo, ordem = ordem,
    X = X, y = y
  ), class = "enpath")
}


predict.enpath <- function(object, newx, ...) {
  sweep(as.matrix(newx) %*% object$beta, 2, object$a0, "+")
}


## AUDITORIA DE OTIMALIDADE (condicoes KKT).
##
## Com g_j = (1/N) <x_j, y - X beta> na escala padronizada, o otimo satisfaz
##    beta_j != 0 :  g_j = lambda alpha pf_j sinal(beta_j) + lambda(1-alpha) pf_j beta_j
##    beta_j == 0 :  |g_j| <= lambda alpha pf_j
##
## Esta funcao NAO faz parte do algoritmo: e o instrumento de verificacao.
## Ela permite afirmar que a implementacao resolve o problema sem comparar com
## outra implementacao -- comparar com o glmnet seria circular, porque nao
## diria qual das duas esta certa.
## `relativa = TRUE` divide a violacao por lambda*alpha, que e a escala natural
## do problema: a restricao para coeficientes nulos e |g_j| <= lambda*alpha*pf_j.
## Uma violacao absoluta de 1e-3 e desprezivel se lambda*alpha = 1 e enorme se
## lambda*alpha = 1e-4, no fim do caminho. Sem normalizar, a leitura engana.
viol_kkt <- function(fit, relativa = FALSE) {
  Xs <- sweep(sweep(fit$X, 2, fit$centro, "-"), 2, fit$escala, "/")
  yc <- fit$y - mean(fit$y)
  N  <- nrow(Xs); al <- fit$alpha
  vapply(seq_along(fit$lambda), function(k) {
    b   <- fit$beta_pad[, k]; lam <- fit$lambda[k]
    g   <- as.vector(crossprod(Xs, yc - Xs %*% b)) / N
    z   <- b == 0
    v   <- numeric(length(b))
    v[z]  <- pmax(abs(g[z]) - lam * al * fit$pf[z], 0)
    v[!z] <- abs(g[!z] - lam * al * fit$pf[!z] * sign(b[!z]) -
                   lam * (1 - al) * fit$pf[!z] * b[!z])
    if (relativa) max(v) / max(lam * al, .Machine$double.eps) else max(v)
  }, numeric(1))
}


## Valor do objetivo penalizado (escala padronizada) em todo o caminho.
##
## Esta e a quantidade certa para comparar duas execucoes do algoritmo. Comparar
## os vetores beta pode enganar: quando a solucao nao e unica -- ou quase nao e --
## dois pontos distintos podem ter o MESMO valor do objetivo. Diferenca em beta
## grande com diferenca em objetivo nula significa degenerescencia, nao erro.
objetivo_caminho <- function(fit) {
  Xs <- sweep(sweep(fit$X, 2, fit$centro, "-"), 2, fit$escala, "/")
  yc <- fit$y - mean(fit$y)
  al <- fit$alpha; N <- nrow(Xs)
  vapply(seq_along(fit$lambda), function(k) {
    b <- fit$beta_pad[, k]; lam <- fit$lambda[k]
    sum((yc - Xs %*% b)^2) / (2 * N) +
      lam * sum(fit$pf * (al * abs(b) + (1 - al) * b^2 / 2))
  }, numeric(1))
}

objetivo <- function(fit, k) {
  objetivo_caminho(fit)[k]
}
