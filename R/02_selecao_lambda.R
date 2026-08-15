# =============================================================================
#  02_selecao_lambda.R
#
#  Escolher lambda e o problema que o caminho barato REDEFINIU. Em 1996, cada
#  valor do parametro exigia uma execucao independente do otimizador: uma
#  validacao cruzada de 5 dobras sobre 15 valores custava 75 execucoes. Por
#  isso Tibshirani (1996) discute, em pe de igualdade, criterios ANALITICOS
#  (validacao cruzada generalizada e uma estimativa nao viesada de risco tipo
#  Stein) que dispensam reamostragem.
#
#  Depois de 2010, o caminho INTEIRO custa da ordem de um unico ajuste, e a
#  validacao cruzada K-fold passa a custar K+1 caminhos. Este arquivo
#  implementa as duas rotas para que a troca possa ser TESTADA, e nao apenas
#  narrada (experimento EC2).
# =============================================================================


## Validacao cruzada K-fold.
##
## Detalhe que importa: a grade de lambda e fixada UMA VEZ, no ajuste com a
## amostra completa, e reusada em todas as dobras. Se cada dobra construisse a
## propria grade (cada uma com seu lambda_max), as curvas nao seriam
## comparaveis ponto a ponto e a media entre dobras nao faria sentido.
cv_en <- function(X, y, alpha = 1, nfolds = 10L, nlambda = 60L, lambda = NULL,
                  folds = NULL, ...) {
  X <- as.matrix(X); y <- as.numeric(y); N <- nrow(X)

  fit <- en_path(X, y, alpha = alpha, lambda = lambda, nlambda = nlambda, ...)
  lam <- fit$lambda

  if (is.null(folds)) folds <- sample(rep(seq_len(nfolds), length.out = N))
  ids <- sort(unique(folds)); nfolds <- length(ids)

  E <- matrix(NA_real_, nfolds, length(lam))
  for (i in seq_along(ids)) {
    tr <- folds != ids[i]
    fk <- en_path(X[tr, , drop = FALSE], y[tr], alpha = alpha, lambda = lam, ...)
    eta <- predict(fk, X[!tr, , drop = FALSE])
    E[i, ] <- colMeans((y[!tr] - eta)^2)
  }

  m <- colMeans(E)
  s <- apply(E, 2, stats::sd) / sqrt(nfolds)
  i_min <- which.min(m)
  ## regra de um erro-padrao: o modelo mais simples cujo erro esta a menos de
  ## um erro-padrao do minimo
  i_1se <- min(which(m <= m[i_min] + s[i_min]))

  list(lambda = lam, cvm = m, cvse = s, erros = E, folds = folds, ajuste = fit,
       df = fit$df, i_min = i_min, i_1se = i_1se,
       lambda_min = lam[i_min], lambda_1se = lam[i_1se])
}


## Coeficientes no lambda escolhido (intercepto na primeira posicao).
coef_em <- function(cv, regra = c("min", "1se")) {
  k <- if (match.arg(regra) == "min") cv$i_min else cv$i_1se
  c(intercepto = cv$ajuste$a0[k], cv$ajuste$beta[, k])
}


## Criterios analiticos ao longo do caminho.
##
## Todos usam graus de liberdade estimados pelo NUMERO DE COEFICIENTES NAO
## NULOS. Isso nao e uma aproximacao grosseira: Zou, Hastie & Tibshirani (2007)
## mostraram que |A| e um estimador nao viesado dos graus de liberdade do
## lasso. Curiosamente, essa justificativa e POSTERIOR ao uso do criterio.
##
##   AIC = N log(RSS/N) + 2 gl
##   BIC = N log(RSS/N) + log(N) gl
##   GCV = (RSS/N) / (1 - gl/N)^2
##
## Quando gl >= N o criterio deixa de ser definido; marcamos Inf, o que na
## pratica proibe esses lambdas -- uma limitacao real dos criterios analiticos
## no regime p > N, que a validacao cruzada nao tem.
criterios_analiticos <- function(fit) {
  N   <- nrow(fit$X)
  eta <- predict(fit, fit$X)
  rss <- colSums((fit$y - eta)^2)
  gl  <- fit$df + 1                                  # + intercepto
  logml <- N * log(pmax(rss / N, .Machine$double.eps))

  aic <- logml + 2 * gl
  bic <- logml + log(N) * gl
  gcv <- (rss / N) / pmax(1 - gl / N, 1e-8)^2
  invalido <- gl >= N
  aic[invalido] <- Inf; bic[invalido] <- Inf; gcv[invalido] <- Inf

  data.frame(lambda = fit$lambda, df = fit$df, rss = rss,
             AIC = aic, BIC = bic, GCV = gcv)
}


## Indice do lambda escolhido por um criterio analitico.
escolhe_por_criterio <- function(fit, criterio = c("AIC", "BIC", "GCV")) {
  criterio <- match.arg(criterio)
  tab <- criterios_analiticos(fit)
  v <- tab[[criterio]]
  if (all(!is.finite(v))) return(length(v))          # degenerado: devolve o menor lambda
  which.min(v)
}
