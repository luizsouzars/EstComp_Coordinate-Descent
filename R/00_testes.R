# =============================================================================
#  00_testes.R -- VERIFICACAO RAPIDA (~20 segundos)
#
#  Rode ANTES de qualquer coisa:   source("R/00_testes.R")
#
#  Cada teste checa uma propriedade que so vale se a implementacao estiver
#  correta. Os testes 2 e 4 sao as validacoes citadas no relatorio: nenhum
#  compara com o glmnet, porque validar uma implementacao contra outra nao diz
#  qual das duas esta certa. A validacao e contra teoria.
# =============================================================================

source("R/00_setup.R")

## Confere que as funcoes carregadas sao mesmo as desta versao. Se a pasta
## estiver misturada com a versao anterior do projeto, os nomes coincidem mas o
## comportamento nao, e os testes passariam a medir outra coisa.
if (!"escala_tol" %in% names(formals(en_path)))
  stop("A funcao en_path carregada nao e a desta versao do projeto. ",
       "Provavelmente ha arquivos de outra versao na pasta R/.", call. = FALSE)
if (exists("huber_path") || exists("en_path_cov_resync"))
  stop("Funcoes da versao anterior do projeto estao carregadas ",
       "(huber_path / en_path_cov_resync). Extraia numa pasta vazia.", call. = FALSE)

.falhas <- 0L
verifica <- function(nome, cond, valor = NULL, esperado = "") {
  ok <- isTRUE(cond)
  if (!ok) .falhas <<- .falhas + 1L
  cat(sprintf("  [%s] %-56s %s\n", if (ok) "OK " else "FALHA", nome,
              if (is.null(valor)) "" else
                sprintf("(%s%s)", paste(format(valor, digits = 3), collapse = ", "),
                        if (nzchar(esperado)) paste0(" ", esperado) else "")))
  invisible(ok)
}

cat("\n=== VERIFICACAO ===\n\n")
set.seed(1234)

cat("1. Soft-thresholding\n")
verifica("S(3, 1) = 2",     isTRUE(all.equal(S(3, 1), 2)))
verifica("S(-3, 1) = -2",   isTRUE(all.equal(S(-3, 1), -2)))
verifica("S(0.5, 1) = 0",   isTRUE(all.equal(S(0.5, 1), 0)))
verifica("S e vetorizado",  length(S(c(-2, 0, 2), 1)) == 3)

cat("\n2. VALIDACAO A -- caso ortonormal reproduz a eq. (3) de 1996\n")
N <- 60L; p <- 8L
Z <- scale(matrix(rnorm(N * p), N, p), TRUE, FALSE)
Q <- qr.Q(qr(Z)) * sqrt(N)             # colunas centradas, ortogonais, norma^2 = N
yq <- as.vector(Q %*% c(3, -2, 1, rep(0, p - 3))) + rnorm(N, sd = 0.5)
lam <- 0.4
bq  <- en_path(Q, yq, alpha = 1, lambda = lam)$beta[, 1]
mqo <- as.vector(crossprod(Q, yq - mean(yq))) / N
verifica("beta = S(MQO, lambda)", max(abs(bq - S(mqo, lam))) < 1e-10,
         max(abs(bq - S(mqo, lam))), "< 1e-10")

cat("\n3. Grade de lambda\n")
d <- gera_ar1(80, 20, rho = 0.5, beta = c(2, -1.5, 1, rep(0, 17)), sigma = 1)
f <- en_path(d$X, d$y, alpha = 1, nlambda = 40)
verifica("df = 0 no maior lambda", f$df[1] == 0, f$df[1], "= 0")
verifica("lambda decrescente", all(diff(f$lambda) < 0))
verifica("df cresce ao longo do caminho", f$df[40] > f$df[1], f$df[40])
verifica("tudo convergiu", all(f$convergiu))

cat("\n4. VALIDACAO B -- condicoes KKT ao longo do caminho\n")
## A violacao KKT fica da ORDEM da tolerancia (tol = 1e-6 sobre |delta beta|),
## nao da ordem de tol^2: o criterio limita o passo, e o residuo KKT herda essa
## escala. Exigir menos que isso seria um teste que nunca passa.
for (a in c(1, 0.5, 0.1)) {
  f_a <- en_path(d$X, d$y, alpha = a, nlambda = 30)
  vr <- max(viol_kkt(f_a, relativa = TRUE))
  verifica(sprintf("KKT relativa a lambda*alpha, alpha = %.1f", a),
           vr < 1e-3, vr, "< 1e-3")
}
## E o ponto do criterio proporcional: a violacao RELATIVA nao deve piorar no
## fim do caminho. Com limiar absoluto ela piora, e muito.
f_prop <- en_path(d$X, d$y, alpha = 1, nlambda = 40)
f_abs  <- en_path(d$X, d$y, alpha = 1, nlambda = 40,
                  tol = sqrt(1e-7), escala_tol = "absoluta")
r_prop <- viol_kkt(f_prop, relativa = TRUE)
r_abs  <- viol_kkt(f_abs,  relativa = TRUE)
verifica("criterio proporcional nao degrada ao longo do caminho",
         max(r_prop) / max(stats::median(r_prop), 1e-300) < 1e3,
         c(mediana = stats::median(r_prop), max = max(r_prop)))
verifica("criterio proporcional e melhor que o absoluto no fim do caminho",
         r_prop[40] < r_abs[40],
         c(proporcional = r_prop[40], absoluto = r_abs[40]))

cat("\n5. ABLACAO -- os interruptores nao mudam a SOLUCAO, so o custo\n")
ref <- en_path(d$X, d$y, alpha = 1, nlambda = 40)
lam40 <- ref$lambda
obj_ref <- objetivo_caminho(ref)
for (cfg in list(c(TRUE, FALSE), c(FALSE, TRUE), c(FALSE, FALSE))) {
  g <- en_path(d$X, d$y, alpha = 1, lambda = lam40,
               warm = cfg[1], ativo = cfg[2])
  ## O criterio de correcao e o VALOR DO OBJETIVO, nao o vetor de coeficientes:
  ## dois pontos distintos podem ser igualmente otimos se a solucao nao for unica.
  do <- max(abs(objetivo_caminho(g) - obj_ref) / pmax(abs(obj_ref), 1e-12))
  verifica(sprintf("warm=%s, ativo=%s atinge o mesmo objetivo", cfg[1], cfg[2]),
           do < 1e-6, do, "< 1e-6")
  verifica(sprintf("warm=%s, ativo=%s da o mesmo beta", cfg[1], cfg[2]),
           max(abs(g$beta - ref$beta)) < 1e-3, max(abs(g$beta - ref$beta)),
           "< 1e-3")
}
sem_ativo <- en_path(d$X, d$y, alpha = 1, lambda = lam40, ativo = FALSE)
verifica("conjunto ativo reduz o numero de visitas",
         sum(ref$visitas) < sum(sem_ativo$visitas),
         c(com = sum(ref$visitas), sem = sum(sem_ativo$visitas)))
frio <- en_path(d$X, d$y, alpha = 1, lambda = lam40, warm = FALSE)
verifica("warm start reduz o numero de ciclos",
         sum(ref$ciclos) < sum(frio$ciclos),
         c(quente = sum(ref$ciclos), frio = sum(frio$ciclos)))

cat("\n6. ORDEM de visita\n")
oa <- en_path(d$X, d$y, alpha = 1, lambda = lam40, ordem = "aleatoria")
verifica("ordem aleatoria resolve o mesmo problema",
         max(abs(oa$beta - ref$beta)) < 1e-5, max(abs(oa$beta - ref$beta)))
op <- en_path(d$X, d$y, alpha = 1, lambda = lam40, ordem = rev(seq_len(20)))
verifica("ordem fixa invertida tambem", max(viol_kkt(op, relativa = TRUE)) < 1e-3,
         max(viol_kkt(op, relativa = TRUE)))
## Armadilha classica: sample(x) com length(x) == 1 amostra de 1:x. Se o
## conjunto ativo tiver uma unica variavel, uma implementacao ingenua visitaria
## coordenadas erradas. Um lambda quase maximo forca |A| = 1.
lam_quase_max <- ref$lambda[1] * 0.97
a1 <- en_path(d$X, d$y, alpha = 1, lambda = lam_quase_max, ordem = "aleatoria")
c1 <- en_path(d$X, d$y, alpha = 1, lambda = lam_quase_max, ordem = "ciclica")
verifica("conjunto ativo unitario com ordem aleatoria",
         max(abs(a1$beta - c1$beta)) < 1e-8 && sum(c1$beta != 0) <= 2,
         c(dif = max(abs(a1$beta - c1$beta)), ativas = sum(c1$beta != 0)))

cat("\n7. Fatores de penalidade e intercepto\n")
pf0 <- rep(1, 20); pf0[5] <- 0
fp <- en_path(d$X, d$y, alpha = 1, nlambda = 20, pf = pf0)
verifica("variavel com pf = 0 entra em todo lambda", all(fp$beta[5, ] != 0))
verifica("KKT respeitam os pesos", max(viol_kkt(fp)) < 1e-5, max(viol_kkt(fp)))
pr <- predict(f, d$X)
verifica("predict tem dimensao N x nlambda", all(dim(pr) == c(80, 40)))
verifica("residuo medio ~ 0 (intercepto correto)",
         abs(mean(d$y - pr[, 30])) < 1e-8, abs(mean(d$y - pr[, 30])))

cat("\n8. Selecao de lambda\n")
cv <- cv_en(d$X, d$y, alpha = 1, nfolds = 5, nlambda = 25)
verifica("lambda_1se >= lambda_min", cv$lambda_1se >= cv$lambda_min)
verifica("modelo do 1 EP nao e maior que o do minimo",
         sum(coef_em(cv, "1se")[-1] != 0) <= sum(coef_em(cv, "min")[-1] != 0))
ff <- sample(rep(1:5, length.out = 80))
verifica("folds fixos dao resultado reprodutivel",
         isTRUE(all.equal(cv_en(d$X, d$y, folds = ff, nlambda = 15)$cvm,
                          cv_en(d$X, d$y, folds = ff, nlambda = 15)$cvm)))
ca <- criterios_analiticos(f)
verifica("criterios analiticos sem NA", !anyNA(ca$AIC) && !anyNA(ca$BIC))
verifica("BIC penaliza mais que AIC quando N > e^2",
         which.min(ca$BIC) <= which.min(ca$AIC),
         c(BIC = which.min(ca$BIC), AIC = which.min(ca$AIC)))

cat("\n9. Rota de 1996\n")
Xs <- sweep(scale(d$X, TRUE, FALSE), 2,
            sqrt(colMeans(scale(d$X, TRUE, FALSE)^2)), "/")
yc <- as.vector(d$y - mean(d$y))
ri <- ridge_iterado(Xs, yc, lambda = 0.15, maxit = 40)
cd <- traco_cd(Xs, yc, lambda = 0.15, maxit = 40)
verifica("descida coordenada decresce o objetivo", all(diff(cd$obj) <= 1e-12))
verifica("ridge iterado nao atinge objetivo menor que o otimo",
         obj_lasso(Xs, yc, ri$beta, 0.15) >= min(cd$obj) - 1e-10,
         c(ridge = obj_lasso(Xs, yc, ri$beta, 0.15), cd = min(cd$obj)))
b1 <- ridge_iterado(Xs, yc, lambda = 0.15, maxit = 40, eps = 1e-2)$beta
b2 <- ridge_iterado(Xs, yc, lambda = 0.15, maxit = 40, eps = 1e-6)$beta
verifica("o limiar eps muda o modelo do ridge iterado",
         sum(b1 != 0) != sum(b2 != 0) || max(abs(b1 - b2)) > 1e-8,
         c(eps_1e2 = sum(b1 != 0), eps_1e6 = sum(b2 != 0)))

cat("\n10. Nao unicidade: copias exatas\n")
ex <- exemplo_copias()
somas <- colSums(ex$lasso)
verifica("a SOMA dos coeficientes das copias e estavel",
         diff(range(somas)) < 1e-6, diff(range(somas)), "< 1e-6")
verifica("a DIVISAO entre as copias depende da ordem",
         diff(range(ex$lasso[1, ])) > 1e-6, diff(range(ex$lasso[1, ])), "> 1e-6")
## exemplo_copias usa tol = 1e-11, entao a igualdade entre as copias deve sair
## na precisao de maquina. Com a tolerancia PADRAO (1e-5) o resto seria da ordem
## de 1e-5 -- correto, mas o teste passaria a medir a tolerancia, nao a unicidade.
verifica("com alpha < 1 a divisao passa a ser unica",
         diff(range(ex$elastic_net[1, ])) < 1e-8,
         diff(range(ex$elastic_net[1, ])), "< 1e-8")

cat("\n11. Visualizacoes em 2 dimensoes\n")
d2t <- dados_2d(N = 80, rho = 0.7, beta = c(1.2, 0.4))
st  <- superficie(d2t, lambda = 0.2, alpha = 1, n = 25)
verifica("superficie: grade correta", all(dim(st$Z) == c(25, 25)) && !anyNA(st$Z))
tt  <- traj_2d(d2t, lambda = 0.2, alpha = 1, inicio = c(0, 0), n_passos = 60)
verifica("trajetoria: objetivo nao cresce", all(diff(tt$objetivo) <= 1e-12),
         max(diff(tt$objetivo)))
verifica("trajetoria e paralela aos eixos (uma coordenada por passo)",
         all(rowSums(abs(diff(as.matrix(tt[, c("b1", "b2")]))) > 1e-14) <= 1))
b_cd <- as.numeric(tt[nrow(tt), c("b1", "b2")])
b_en <- en_path(d2t$X, d2t$y, alpha = 1, lambda = 0.2, padronizar = FALSE)$beta[, 1]
verifica("traj_2d converge para o mesmo ponto que en_path",
         max(abs(b_cd - b_en)) < 1e-5, max(abs(b_cd - b_en)), "< 1e-5")
cw_t <- comparacao_warm(d2t, nlambda = 10)
pq <- cw_t$passos$passos[cw_t$passos$partida == "warm start (solucao anterior)"]
pf <- cw_t$passos$passos[cw_t$passos$partida == "partida fria (do zero)"]
verifica("warm start custa menos que a partida fria na quase totalidade",
         mean(pq <= pf) >= 0.9, c(total_quente = sum(pq), total_frio = sum(pf)))
vn <- vale_nao_unicidade(N = 80, lambda = 0.25, n = 25)
verifica("vale de nao unicidade: superficies geradas",
         !anyNA(vn$lasso$Z) && !anyNA(vn$enet$Z))

dcp <- decomposicao(d2t, lambda = 0.3, alpha = 1, n = 20)
verifica("decomposicao: soma = quadratico + penalidade",
         max(abs(dcp$soma$Z - dcp$quadratico$Z - dcp$penalidade$Z)) < 1e-12)
verifica("penalidade e minima na origem",
         which.min(dcp$penalidade$Z) == which.min(abs(dcp$penalidade$b1))
           + (which.min(abs(dcp$penalidade$b2)) - 1) * length(dcp$penalidade$b1))
ct <- corte_1d(d2t, lambdas = c(0.05, 0.6), alpha = 1)
verifica("corte 1d: o minimo da curva bate com a solucao do en_path",
         all(vapply(split(ct, ct$lambda), function(z)
           abs(z$beta[which.min(z$objetivo)] - z$otimo[1]) < 0.02, logical(1))))
pv <- perfil_vale(N = 80, lambda = 0.25, alphas = c(1, 0.5))
amp <- vapply(split(pv, pv$alpha), function(z) diff(range(z$objetivo)), numeric(1))
verifica("perfil do vale: plano no lasso, curvo na elastic net",
         amp[["1"]] < 1e-8 && amp[["0.5"]] > 1e-4, amp)

cat("\n12. Diagnostico de ordem\n")
dg_u <- diagnostico_ordem(d$X, d$y, alpha = 1, B = 6, nlambda = 20)
verifica("solucao unica -> Jaccard = 1", dg_u$jaccard_limiar > 0.999,
         dg_u$jaccard_limiar)
verifica("solucao unica -> objetivo constante entre reordenacoes",
         dg_u$amplitude_obj < 1e-9, dg_u$amplitude_obj, "< 1e-9")
Xd <- cbind(d$X, d$X[, 1], d$X[, 2])          # duas colunas duplicadas
dg_d <- diagnostico_ordem(Xd, d$y, alpha = 1, B = 6, nlambda = 20)
verifica("colunas duplicadas -> objetivo ainda constante",
         dg_d$amplitude_obj < 1e-8, dg_d$amplitude_obj, "< 1e-8")

## a grade do diagnostico deve TERMINAR no lambda alvo -- e isso que o torna
## barato o suficiente para rodar centenas de vezes
lam_alvo <- 0.3 * lambda_max(d$X, d$y, alpha = 1)
dg_a <- diagnostico_ordem(d$X, d$y, alpha = 1, lambda_alvo = lam_alvo,
                          B = 4, nlambda = 12)
verifica("diagnostico: a grade termina no lambda alvo",
         abs(dg_a$lambda - lam_alvo) / lam_alvo < 1e-10, dg_a$lambda)
verifica("lambda_max coincide com o primeiro lambda de en_path",
         abs(lambda_max(d$X, d$y, alpha = 1) -
             en_path(d$X, d$y, alpha = 1, nlambda = 5)$lambda[1]) < 1e-12)
verifica("guarda de estagnacao nao impede a convergencia normal",
         all(en_path(d$X, d$y, alpha = 1, nlambda = 20)$convergiu))

cat("\n13. Metricas\n")
m <- metricas_selecao(c(1, 0, 1, 0, 1), c(1, 1, 0, 0, 1))
verifica("TPR = 2/3", isTRUE(all.equal(unname(m["TPR"]), 2/3)), m["TPR"])
verifica("FPR = 1/2", isTRUE(all.equal(unname(m["FPR"]), 1/2)), m["FPR"])
verifica("Jaccard({1,2},{2,3}) = 1/3",
         isTRUE(all.equal(jaccard(c(1, 2), c(2, 3)), 1/3)))
verifica("estabilidade de conjuntos iguais = 1",
         isTRUE(all.equal(estabilidade(list(1:3, 1:3, 1:3)), 1)))

cat("\n14. Casos-limite\n")
verifica("p = 1", !anyNA(en_path(matrix(rnorm(50)), rnorm(50), nlambda = 5)$beta))
Xz <- d$X; Xz[, 3] <- 0
verifica("coluna constante", !anyNA(en_path(Xz, d$y, nlambda = 10)$beta))
verifica("coluna duplicada",
         !anyNA(en_path(cbind(d$X, d$X[, 1]), d$y, nlambda = 10)$beta))
verifica("p > N", !anyNA(en_path(matrix(rnorm(30 * 80), 30, 80), rnorm(30),
                                 nlambda = 10)$beta))
verifica("alpha = 0 (ridge puro)",
         !anyNA(en_path(d$X, d$y, alpha = 0, nlambda = 10)$beta))
verifica("y constante devolve solucao nula",
         all(en_path(d$X, rep(2, 80), nlambda = 5)$beta == 0))

cat("\n===========================================\n")
if (.falhas == 0L) cat(" TODAS AS VERIFICACOES PASSARAM.\n") else
  cat(sprintf(" %d VERIFICACAO(OES) FALHARAM.\n", .falhas))
cat("===========================================\n")
invisible(.falhas)
