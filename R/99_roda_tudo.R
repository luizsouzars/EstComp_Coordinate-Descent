# =============================================================================
#  99_roda_tudo.R
#
#  Roda todos os experimentos e guarda os resultados em resultados/*.rds.
#  O relatorio (relatorio.Rmd) le esse cache; se ele nao existir, o proprio
#  relatorio calcula na hora do knit. Ou seja: este script e um atalho, nao
#  uma dependencia.
#
#      source("R/99_roda_tudo.R")                    # modo rapido
#      RAPIDO <- FALSE; source("R/99_roda_tudo.R")   # estudo completo
# =============================================================================

if (!exists("RAPIDO")) RAPIDO <- TRUE
source("R/00_setup.R")

t0 <- Sys.time()
cat(sprintf("\nModo: %s | inicio: %s\n\n", if (RAPIDO) "RAPIDO" else "COMPLETO",
            format(t0, "%H:%M:%S")))

## --- verificacao antes de gastar tempo simulando -----------------------------
falhas <- tryCatch(source("R/00_testes.R")$value,
                   error = function(e) { message("ERRO nos testes: ",
                                                 conditionMessage(e)); NA })
if (!identical(falhas, 0L))
  cat("\n*** A verificacao nao passou limpa. Os experimentos vao rodar mesmo",
      "assim, mas confira os avisos acima. ***\n\n")

## --- experimentos (parametros vem de PARAM, definido em 00_setup.R) ----------
abl <- com_cache("ec1_ablacao",
  experimento_ablacao(R = PARAM$R_ABL, rhos = PARAM$RHO_ABL,
                      nlambda = PARAM$NL_ABL))

cri <- com_cache("ec2_criterios", experimento_criterios(R = PARAM$R_CRI))

rot <- com_cache("rotas_1996", experimento_rotas())

vel <- com_cache("co_velocidade",
  experimento_ordem_velocidade(R = PARAM$R_VEL))

sel <- com_cache("co_selecao",
  experimento_ordem_selecao(R = PARAM$R_SEL, B = PARAM$B_ORD))

cop <- com_cache("co_copias", exemplo_copias())

app <- com_cache("aplicacao",
  experimento_aplicacao(com_cache("dados_crime", carrega_dados()),
                        B_ordem = PARAM$B_ORD))

## --- resumo no console -------------------------------------------------------
cat("\n\n==================== RESUMO ====================\n")

cat("\nEC1 -- ablacao: custo em visitas de coordenada (media +/- EP)\n")
print(resume_mc(abl, "visitas", c("cenario", "rotulo")), row.names = FALSE)
cat(sprintf("\nDiferenca maxima entre os VETORES de coeficientes: %.3e\n",
            max(abl$difer)))
cat(sprintf("Diferenca maxima entre os VALORES DO OBJETIVO (relativa): %.3e\n",
            max(abl$dif_obj)))
cat(sprintf("Violacao KKT maxima: %.3e absoluta | %.3e relativa\n",
            max(abl$kkt), max(abl$kkt_rel)))

cat("\nEC2 -- selecao de lambda: erro contra E(y|x)\n")
print(resume_mc(cri, "eqm_mu", c("rho", "criterio")), row.names = FALSE)

cat("\nCO -- efeito da ordem sobre o custo\n")
print(resume_mc(vel, "visitas", c("rho", "ordem")), row.names = FALSE)

cat("\nCO -- efeito da ordem sobre a selecao (Jaccard: 1 = mesma lista sempre)\n")
print(resume_mc(sel, "jaccard", c("rho_g", "alpha")), row.names = FALSE)
cat("\nCO -- amplitude do objetivo (proxima de 0 => os pontos sao igualmente otimos)\n")
print(resume_mc(sel, "amplitude_obj", c("rho_g", "alpha")), row.names = FALSE)

cat("\nAplicacao -- desempenho fora da amostra\n")
print(app$desempenho, row.names = FALSE)
cat(sprintf("\nVariaveis intercambiaveis (lasso): %d de %d\n",
            app$diag_lasso$intercambiaveis, app$p))

if (length(.TEMPOS)) {
  cat("\nTempo por experimento (minutos):\n")
  print(round(unlist(.TEMPOS), 2))
}

cat(sprintf("\n\nTempo total: %.1f minutos\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Agora renderize: rmarkdown::render(\"relatorio.Rmd\")\n")
