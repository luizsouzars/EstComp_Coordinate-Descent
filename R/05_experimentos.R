# =============================================================================
#  05_experimentos.R -- ESTUDO DE MONTE CARLO
#
#  Dois experimentos, ambos sobre o MESMO objeto (a funcao en_path).
#
#  EC1  ABLACAO. Liga e desliga os dois ingredientes que fazem do caminho algo
#       barato -- inicializacao quente e conjunto ativo -- DENTRO da mesma
#       implementacao. Isso e mais informativo que comparar com outro software:
#       nao ha diferenca de linguagem, de compilador ou de estilo de codigo
#       para confundir. Duas perguntas:
#         (a) as quatro configuracoes dao a MESMA solucao? (elas devem: sao
#             aceleracoes, nao estimadores diferentes)
#         (b) quanto cada ingrediente economiza, e como isso varia com N, p e
#             a correlacao?
#       O artigo afirma que o warm start e essencial. Nunca mede o quanto.
#
#  EC2  O QUE O CAMINHO BARATO PERMITIU. Em 1996, criterios analiticos (GCV,
#       risco de Stein) eram atrativos porque a validacao cruzada custava uma
#       execucao completa do otimizador por valor do parametro. Depois de 2010
#       a CV custa K+1 caminhos. A troca compensou? Comparamos CV-min, regra de
#       1 EP, AIC, BIC e GCV pelo erro fora da amostra e pela selecao.
#
#  Ambos devolvem data.frames; quem plota e o relatorio.
# =============================================================================


# ---------------------------------------------------------------------------
# EC1 -- ablacao dos ingredientes do caminho
# ---------------------------------------------------------------------------
## <<ec1
## lambda_min_ratio = 1e-3: o criterio de parada agora escala com lambda, e o
## fim do caminho (lambda quatro ordens abaixo do maximo) ficaria caro sem
## acrescentar nada ao que o experimento quer medir.
## PARES (N, p) explicitos, em vez de cruzamento. Cruzar N em {120, 320} com p
## em {30, 90} produzia quatro celulas TODAS com N > p, e o cenario p > N
## desaparecia da tabela sem aviso. Com pares, os dois regimes ficam
## representados duas vezes cada.
##
## lambda_min_ratio = 0,05: o criterio de parada escala com lambda, e a ultima
## decada do caminho e onde ele fica mais caro sem acrescentar nada ao que este
## experimento mede -- a comparacao entre configuracoes do MESMO algoritmo.
experimento_ablacao <- function(R = 3L,
                                pares = list(c(300, 60), c(300, 150),
                                             c(100, 200), c(60, 300)),
                                rhos = c(0, 0.9), nlambda = 20L,
                                lambda_min_ratio = 0.05, tol = 1e-5,
                                semente = 2026) {
  grade <- expand.grid(rep = seq_len(R), rho = rhos, par = seq_along(pares))
  grade$N <- vapply(pares[grade$par], `[`, numeric(1), 1)
  grade$p <- vapply(pares[grade$par], `[`, numeric(1), 2)
  cfg <- data.frame(
    warm  = c(TRUE,  TRUE,  FALSE, FALSE),
    ativo = c(TRUE,  FALSE, TRUE,  FALSE),
    rotulo = c("warm start + active set", "so warm start",
               "so active set", "nenhum dos dois")
  )

  linhas <- vector("list", nrow(grade) * nrow(cfg))
  z <- 0L
  for (i in seq_len(nrow(grade))) {
    progresso(i, nrow(grade), "EC1 celula")
    g <- grade[i, ]
    set.seed(semente + i)                       # mesma amostra para as 4 configuracoes
    d <- gera_equicor(g$N, g$p, g$rho)

    ## a grade de lambda e definida uma unica vez e reusada: sem isso as
    ## configuracoes resolveriam problemas diferentes
    ref <- en_path(d$X, d$y, alpha = 1, nlambda = nlambda,
                   lambda_min_ratio = lambda_min_ratio, tol = tol,
                   warm = TRUE, ativo = TRUE)
    lam <- ref$lambda

    obj_ref <- objetivo_caminho(ref)

    for (j in seq_len(nrow(cfg))) {
      t0 <- Sys.time()
      f <- en_path(d$X, d$y, alpha = 1, lambda = lam, tol = tol,
                   warm = cfg$warm[j], ativo = cfg$ativo[j])
      tempo <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      z <- z + 1L
      linhas[[z]] <- data.frame(
        g, warm = cfg$warm[j], ativo = cfg$ativo[j], rotulo = cfg$rotulo[j],
        visitas = sum(f$visitas),              # custo independente de maquina
        ciclos  = sum(f$ciclos),
        tempo   = tempo,
        kkt     = max(viol_kkt(f)),
        kkt_rel = max(viol_kkt(f, relativa = TRUE)),
        ## DUAS medidas de concordancia com a referencia, e a distincao importa:
        ##  difer   -- distancia entre os VETORES de coeficientes;
        ##  dif_obj -- distancia entre os VALORES DO OBJETIVO (relativa).
        ## Se difer for grande e dif_obj for nula, os dois pontos sao igualmente
        ## otimos: a solucao nao e unica. Isso e degenerescencia do problema, nao
        ## erro do algoritmo -- e so o objetivo distingue os dois casos.
        difer   = max(abs(f$beta - ref$beta)),
        dif_obj = max(abs(objetivo_caminho(f) - obj_ref) / pmax(abs(obj_ref), 1e-12)),
        convergiu_tudo = all(f$convergiu)
      )
    }
  }
  out <- do.call(rbind, linhas)
  out$cenario <- ifelse(out$N > out$p, "N > p", "p > N")
  out$rotulo  <- factor(out$rotulo, levels = cfg$rotulo)
  out
}
## ec1>>


# ---------------------------------------------------------------------------
# EC2 -- validacao cruzada contra os criterios analiticos de 1996
# ---------------------------------------------------------------------------
## <<ec2
experimento_criterios <- function(R = 20L, N = 150L, p = 60L,
                                  rhos = c(0, 0.5, 0.9), nlambda = 50L,
                                  lambda_min_ratio = 1e-3,
                                  nfolds = 5L, snr = 2, n_teste = 2000L,
                                  semente = 4090) {
  ## Efeitos em ESCADA, do forte ao quase indistinguivel do ruido. Com seis
  ## efeitos todos fortes, qualquer criterio acerta todos e a taxa de verdadeiros
  ## positivos fica presa em 1, sem poder de discriminacao. O gradiente de
  ## magnitudes e o que torna a comparacao entre criterios informativa.
  beta <- numeric(p)
  beta[seq(1, by = 10, length.out = 6)] <- c(2, -1.5, 1, -0.6, 0.35, -0.25)

  linhas <- list(); z <- 0L
  n_total <- length(rhos) * R; feitos <- 0L
  for (rho in rhos) {
    Sig <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
    ## sigma fixado para manter a razao sinal-ruido constante entre cenarios:
    ## sem isso, mudar rho mudaria duas coisas ao mesmo tempo
    sigma <- sqrt(as.numeric(t(beta) %*% Sig %*% beta)) / snr

    for (r in seq_len(R)) {
      feitos <- feitos + 1L; progresso(feitos, n_total, "EC2 replica")
      set.seed(semente + 977 * r + round(100 * rho))
      d <- gera_ar1(N, p, rho, beta, sigma = sigma, n_teste = n_teste)

      cv  <- cv_en(d$X, d$y, alpha = 1, nfolds = nfolds, nlambda = nlambda,
                   lambda_min_ratio = lambda_min_ratio)
      fit <- cv$ajuste

      escolhas <- c(
        "CV (minimo)"  = cv$i_min,
        "CV (1 EP)"    = cv$i_1se,
        "AIC"          = escolhe_por_criterio(fit, "AIC"),
        "BIC"          = escolhe_por_criterio(fit, "BIC"),
        "GCV"          = escolhe_por_criterio(fit, "GCV")
      )

      for (nm in names(escolhas)) {
        k   <- escolhas[[nm]]
        b   <- fit$beta[, k]
        eta <- as.vector(fit$a0[k] + d$X_teste %*% b)
        z <- z + 1L
        linhas[[z]] <- data.frame(
          rho = rho, rep = r, criterio = nm,
          as.list(metricas_selecao(b, beta)),
          ## erro contra a media verdadeira E(y|x): mede o ajuste, sem o ruido
          ## irredutivel do teste, o que reduz a variancia da comparacao
          eqm_mu  = mean((d$mu_teste - eta)^2),
          erro_l2 = sqrt(sum((b - beta)^2)),
          lambda  = fit$lambda[k],
          ## um criterio que escolhe a ULTIMA casa da grade esta dizendo que
          ## queria ir mais longe: a escolha e da grade, nao do criterio
          no_limite = as.integer(k == length(fit$lambda))
        )
      }
    }
  }
  out <- do.call(rbind, linhas)
  out$criterio <- factor(out$criterio,
    levels = c("CV (minimo)", "CV (1 EP)", "AIC", "BIC", "GCV"))
  out
}
## ec2>>


# ---------------------------------------------------------------------------
# Ilustracao: as tres rotas para o mesmo lasso (2010 x 1996)
# ---------------------------------------------------------------------------
experimento_rotas <- function(N = 150L, p = 40L, rho = 0.5, lambda = 0.12,
                              passos = 40L, semente = 7) {
  set.seed(semente)
  d  <- gera_ar1(N, p, rho, c(3, -2, 1.5, 0, 0, 1, rep(0, p - 6)), sigma = 2)
  Xc <- scale(d$X, TRUE, FALSE)
  Xs <- sweep(Xc, 2, sqrt(colMeans(Xc^2)), "/")
  yc <- as.vector(d$y - mean(d$y))

  ## REFERENCIA INDEPENDENTE. Usar min(cd$obj) como otimo faria a rota de 2010
  ## servir de referencia para si mesma: a curva dela cairia a zero por
  ## construcao, nao por merito. O otimo vem de um ajuste separado com tolerancia
  ## de 1e-12, e assim as duas rotas sao medidas contra o mesmo alvo externo.
  b_ref <- en_path(Xs, yc, alpha = 1, lambda = lambda, padronizar = FALSE,
                   tol = 1e-12, maxit = 100000L)$beta[, 1]
  otimo <- obj_lasso(Xs, yc, b_ref, lambda)

  cd <- traco_cd(Xs, yc, lambda, maxit = passos)
  ri <- ridge_iterado(Xs, yc, lambda = lambda, maxit = passos, eps = 1e-4)
  obj_ri <- apply(ri$traj, 1, function(b) obj_lasso(Xs, yc, b, lambda))

  traj <- rbind(
    data.frame(passo = seq_len(passos), gap = cd$obj - otimo,
               rota = "descida coordenada (2010)"),
    data.frame(passo = seq_len(passos), gap = obj_ri - otimo,
               rota = "ridge iterado (1996, Secao 4)")
  )
  traj$gap <- pmax(traj$gap, 1e-16)

  ## Sensibilidade ao limiar arbitrario da rota de 1996. Todas as distancias sao
  ## medidas contra b_ref, a solucao de referencia -- inclusive a do proprio
  ## coordinate descent, que assim tambem e avaliada e nao apenas usada como
  ## regua.
  eps_grade <- c(1e-2, 1e-3, 1e-4, 1e-6)
  sens <- do.call(rbind, lapply(eps_grade, function(e) {
    b <- ridge_iterado(Xs, yc, lambda = lambda, maxit = passos, eps = e)$beta
    data.frame(eps = e, nao_nulos = sum(b != 0),
               dist_ao_otimo = max(abs(b - b_ref)),
               objetivo = obj_lasso(Xs, yc, b, lambda),
               excesso = obj_lasso(Xs, yc, b, lambda) - otimo)
  }))
  sens <- rbind(sens, data.frame(
    eps = NA, nao_nulos = sum(cd$beta != 0),
    dist_ao_otimo = max(abs(cd$beta - b_ref)),
    objetivo = min(cd$obj), excesso = min(cd$obj) - otimo))

  ## invariante: nenhuma rota pode ficar ABAIXO do otimo de referencia
  stopifnot(all(sens$excesso > -1e-9))

  list(traj = traj, sensibilidade = sens, otimo = otimo,
       nao_nulos_ref = sum(b_ref != 0))
}
