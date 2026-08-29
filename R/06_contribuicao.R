# =============================================================================
#  06_contribuicao.R -- CONTRIBUICAO ORIGINAL
#
#  A LACUNA. O artigo de 2010 varre as coordenadas na ordem j = 1, 2, ..., p e
#  nunca justifica essa escolha; ela aparece como detalhe de implementacao. Mas
#  a ordem de visita tem consequencias que o artigo nao avalia.
#
#  O QUE A TEORIA DIZ. Com colunas exatamente linearmente dependentes -- copias,
#  variaveis derivadas, indicadoras redundantes, ou p > N -- o conjunto de
#  solucoes otimas do lasso e uma FACE do poliedro, nao um ponto: qualquer
#  divisao do coeficiente entre as colunas dependentes tem o mesmo ajuste e a
#  mesma norma l1, logo o mesmo objetivo. A descida coordenada devolve UM ponto
#  dessa face, e qual ponto depende de quem foi visitado primeiro.
#
#  Quando as colunas sao apenas MUITO correlacionadas, mas distintas, a solucao
#  e unica com probabilidade 1 (R. J. Tibshirani, 2013). Nossos proprios
#  resultados confirmam isso: com correlacao intragrupo 0,99 a reordenacao quase
#  nao muda nada. A degenerescencia nao e um fenomeno de "correlacao alta"; e um
#  fenomeno de DEPENDENCIA EXATA.
#
#  MAS existe um segundo mecanismo, que se confunde com o primeiro na pratica:
#  com tolerancia finita, coeficientes numericamente minusculos oscilam entre
#  zero e nao zero conforme a ordem de visita. Para quem le a saida, isso e
#  indistinguivel de degenerescencia verdadeira -- e o remedio e completamente
#  diferente (apertar a tolerancia, e nao trocar de estimador).
#
#  A PROPOSTA. Um diagnostico que SEPARA os dois. Reajustamos B vezes com ordens
#  de visita sorteadas e olhamos duas coisas:
#
#    - o CONJUNTO selecionado varia?
#    - o VALOR DO OBJETIVO varia?
#
#  objetivo identico + conjunto diferente  ->  nao unicidade verdadeira
#  objetivo diferente                      ->  convergencia incompleta
#
#  Essa distincao e o ponto original do trabalho. Reamostragem bootstrap, o
#  instrumento usual para avaliar estabilidade de selecao, nao consegue faze-la:
#  ela muda os dados, e portanto mistura degenerescencia, erro numerico e
#  incerteza amostral numa unica medida.
#
#  Testamos ainda o remedio conhecido para o primeiro mecanismo -- reduzir alpha,
#  o que torna o objetivo estritamente convexo e a solucao unica -- medindo o
#  preco pago em erro de predicao.
# =============================================================================


## DIAGNOSTICO DE NAO UNICIDADE POR REORDENACAO.
##
## Ajusta B vezes o mesmo caminho, cada vez com uma permutacao sorteada da ordem
## de visita, e resume a variacao no conjunto selecionado E no valor do objetivo.
##
## `limiar_rel` define o que conta como "selecionada": |beta_j| acima de
## limiar_rel vezes o maior coeficiente em modulo. Sem esse limiar, um
## coeficiente de 1e-12 conta como variavel do modelo e o diagnostico mede ruido
## de arredondamento em vez de estrutura.
##
## Leitura do resultado:
##   jaccard_limiar = 1                      -> selecao estavel
##   jaccard_limiar < 1 e amplitude_obj ~ 0  -> NAO UNICIDADE: pontos distintos,
##                                              igualmente otimos
##   jaccard_limiar < 1 e amplitude_obj > 0  -> CONVERGENCIA INCOMPLETA: aperte a
##                                              tolerancia antes de concluir
##                                              qualquer coisa
## A grade TERMINA no lambda alvo: nao ha razao para percorrer o caminho alem
## dele. Isso importa muito mais do que parece. Com colunas exatamente
## dependentes, os lambdas do fim do caminho sao justamente os que nao
## convergem, e o diagnostico refaz o caminho B vezes -- prolongar a grade
## multiplica um custo patologico por B. Terminar no alvo tornou este
## experimento cerca de dez vezes mais rapido.
diagnostico_ordem <- function(X, y, alpha = 1, lambda_alvo = NULL, B = 20L,
                              nlambda = 20L, limiar_rel = 1e-6, semente = 1L,
                              nomes = NULL, ...) {
  X <- as.matrix(X)
  p <- ncol(X)
  if (is.null(nomes)) nomes <- colnames(X)
  if (is.null(nomes)) nomes <- paste0("V", seq_len(p))

  lmax <- lambda_max(X, y, alpha = alpha)
  alvo <- if (is.null(lambda_alvo)) {
    0.05 * lmax
  } else {
    min(lambda_alvo, 0.98 * lmax)
  }
  lam <- exp(seq(log(lmax), log(alvo), length.out = nlambda))
  k <- nlambda # o alvo e, por construcao, o ultimo ponto

  ## caminho de referencia (ordem ciclica), na mesma grade
  ref <- en_path(X, y, alpha = alpha, lambda = lam, ...)

  selecionadas <- function(b) {
    corte <- limiar_rel * max(abs(b))
    which(abs(b) > corte)
  }

  set.seed(semente)
  conj_exato <- conj_limiar <- vector("list", B)
  coefs <- matrix(0, p, B)
  objs <- numeric(B)
  for (b in seq_len(B)) {
    ord <- sample.int(p) # uma ordem fixa por reajuste
    f <- en_path(X, y, alpha = alpha, lambda = lam, ordem = ord, ...)
    bb <- f$beta[, k]
    coefs[, b] <- bb
    conj_exato[[b]] <- which(bb != 0)
    conj_limiar[[b]] <- selecionadas(bb)
    objs[b] <- objetivo_caminho(f)[k]
  }

  freq <- rowMeans(apply(coefs, 2, function(b) {
    z <- numeric(p)
    z[selecionadas(b)] <- 1
    z
  }))
  amp_obj <- (max(objs) - min(objs)) / max(abs(mean(objs)), 1e-12)

  list(
    lambda = lam[k],
    alpha = alpha,
    freq = stats::setNames(freq, nomes),
    jaccard_exato = estabilidade(conj_exato),
    jaccard_limiar = estabilidade(conj_limiar),
    amplitude_obj = amp_obj, # o discriminante entre os dois mecanismos
    tamanho_medio = mean(vapply(conj_limiar, length, numeric(1))),
    sempre = sum(freq == 1),
    nunca = sum(freq == 0),
    intercambiaveis = sum(freq > 0 & freq < 1),
    dp_coef = stats::setNames(apply(coefs, 1, stats::sd), nomes),
    coefs = coefs, objs = objs, conjuntos = conj_limiar, ref = ref, k = k
  )
}


# ---------------------------------------------------------------------------
# (i) A ordem afeta a VELOCIDADE?
# ---------------------------------------------------------------------------
## <<ordem_velocidade
## A pergunta e se ALEATORIZAR a ordem de varredura compensa, e a resposta
## depende de rho -- entao a grade de rho precisa ser fina o bastante para
## mostrar ONDE o efeito aparece e onde ele se apaga. Com tres pontos so se ve
## que o sinal muda; com oito se ve a forma da curva.
##
## PAREAMENTO EM DOIS NIVEIS, e e dele que vem toda a precisao deste experimento:
##   (i) dentro de uma replica, as duas ordens rodam sobre os MESMOS dados e a
##       MESMA grade de lambda. A razao aleatoria/ciclica existe POR REPLICA --
##       calcule-a assim, nunca dividindo as medias, que e o que descarta o
##       pareamento e deixa a razao sem erro-padrao.
##  (ii) a semente depende so da replica, NAO de rho: para um mesmo r, todos os
##       valores de rho partem do mesmo sorteio de base e diferem apenas na
##       estrutura de correlacao. Isso torna a CURVA em rho pareada tambem, e
##       elimina a chance de colisao de sementes que a formula antiga tinha ao
##       adensar a grade.
##
## `obj` e `convergiu` sao instrumentacao de auditoria: sem eles, comparar
## visitas entre um ajuste convergido e outro que estourou maxit nao mede nada.
## O mesmo cuidado que o EC1 ja exige das quatro configuracoes.
experimento_ordem_velocidade <- function(R = 5L,
                                         rhos = c(
                                           0, 0.3, 0.5, 0.7,
                                           0.8, 0.9, 0.95, 0.99
                                         ),
                                         N = 100L, p = 200L, nlambda = 30L,
                                         semente = 8100) {
  linhas <- list()
  z <- 0L
  n_total <- length(rhos) * R
  feitos <- 0L
  for (rho in rhos) {
    for (r in seq_len(R)) {
      feitos <- feitos + 1L
      progresso(feitos, n_total, "velocidade")
      set.seed(semente + r) # ver (ii) acima
      d <- gera_equicor(N, p, rho)

      ## grade fixa: as duas ordens resolvem exatamente os mesmos problemas
      lam <- en_path(d$X, d$y, alpha = 1, nlambda = nlambda)$lambda

      for (ord in c("ciclica", "aleatoria")) {
        t0 <- Sys.time()
        f <- en_path(d$X, d$y, alpha = 1, lambda = lam, ordem = ord)
        z <- z + 1L
        linhas[[z]] <- data.frame(
          rho = rho, rep = r, ordem = ord,
          visitas = sum(f$visitas), ciclos = sum(f$ciclos),
          tempo = as.numeric(difftime(Sys.time(), t0, units = "secs")),
          kkt = max(viol_kkt(f, relativa = TRUE)),
          obj = sum(objetivo_caminho(f)),
          convergiu = all(f$convergiu)
        )
      }
    }
  }
  do.call(rbind, linhas)
}
## ordem_velocidade>>


# ---------------------------------------------------------------------------
# (ii) A ordem afeta QUAL VARIAVEL e selecionada? E o remedio, quanto custa?
# ---------------------------------------------------------------------------
## <<ordem_selecao
## rho_g = 1 produz colunas EXATAMENTE identicas dentro de cada grupo (o caso em
## que a teoria garante nao unicidade); rho_g < 1 produz colunas distintas, em
## que a solucao e unica com probabilidade 1. Cruzar os dois com alpha separa o
## efeito da dependencia exata do efeito da correlacao alta.
experimento_ordem_selecao <- function(R = 10L, alphas = c(1, 0.95, 0.8, 0.5),
                                      rho_gs = c(0.9, 0.99, 1), N = 100L,
                                      B = 15L, n_teste = 2000L, semente = 5150) {
  linhas <- list()
  z <- 0L
  n_total <- R * length(rho_gs)
  feitos <- 0L
  for (r in seq_len(R)) {
    for (rg in rho_gs) {
      feitos <- feitos + 1L
      progresso(feitos, n_total, "selecao")
      set.seed(semente + 61 * r + round(100 * rg))
      d <- gera_grupos(N, g = 3L, tam = 5L, n_ruido = 15L, rho_g = rg)
      dt <- gera_grupos(n_teste, g = 3L, tam = 5L, n_ruido = 15L, rho_g = rg)
      mu_teste <- as.vector(dt$X %*% d$beta)

      for (a in alphas) {
        dg <- diagnostico_ordem(d$X, d$y,
          alpha = a, B = B, nlambda = 30L,
          semente = semente + r
        )
        b <- dg$ref$beta[, dg$k]
        eta <- as.vector(dg$ref$a0[dg$k] + dt$X %*% b)

        z <- z + 1L
        linhas[[z]] <- data.frame(
          rep = r, rho_g = rg, alpha = a,
          jaccard = dg$jaccard_limiar, # 1 = mesmo conjunto sempre
          jaccard_exato = dg$jaccard_exato,
          amplitude_obj = dg$amplitude_obj, # ~0 => os pontos sao igualmente otimos
          intercambiaveis = dg$intercambiaveis,
          tamanho = dg$tamanho_medio,
          recuperacao = mean(d$idx_efeito %in% which(b != 0)),
          eqm_mu = mean((mu_teste - eta)^2) # o preco do remedio
        )
      }
    }
  }
  do.call(rbind, linhas)
}
## ordem_selecao>>


## Ilustracao minima e completamente transparente: tres COPIAS EXATAS de um
## mesmo preditor. Aqui a nao unicidade e demonstravel no papel -- qualquer
## divisao (b1, b2, b3) com b1+b2+b3 = c e mesmo sinal tem o mesmo l1 e o mesmo
## ajuste, logo o mesmo objetivo. O algoritmo devolve uma delas.
## tol apertado de proposito: este exemplo tem N = 150 e p = 8, custa
## milissegundos, e a mensagem dele e uma IGUALDADE exata entre os coeficientes
## das copias quando alpha < 1. Com a tolerancia padrao a igualdade sai com erro
## da ordem de 1e-5 -- correto, mas confuso numa tabela que quer mostrar
## simetria perfeita.
exemplo_copias <- function(N = 150L, semente = 11, tol = 1e-11) {
  set.seed(semente)
  x <- stats::rnorm(N)
  X <- cbind(x, x, x, matrix(stats::rnorm(N * 5), N, 5))
  colnames(X) <- c("copia 1", "copia 2", "copia 3", paste("ruido", 1:5))
  y <- 3 * x + stats::rnorm(N)

  ordens <- list(`1,2,3` = 1:8, `3,2,1` = c(3, 2, 1, 4:8), `2,1,3` = c(2, 1, 3, 4:8))
  tab <- sapply(ordens, function(o) {
    en_path(X, y, alpha = 1, lambda = 0.5, ordem = o, tol = tol)$beta[1:3, 1]
  })
  rownames(tab) <- colnames(X)[1:3]

  ## com alpha < 1 o problema fica estritamente convexo: solucao unica
  tab_en <- sapply(ordens, function(o) {
    en_path(X, y, alpha = 0.5, lambda = 0.5, ordem = o, tol = tol)$beta[1:3, 1]
  })
  rownames(tab_en) <- colnames(X)[1:3]

  list(lasso = tab, elastic_net = tab_en, X = X, y = y)
}
