## Verificação independente do experimento das rotas (cole no console).
## Não usa nada do cache: recalcula e confere as invariantes na sua frente.
source("R/00_setup.R")
r <- experimento_rotas()

cat("\n1) Nenhuma rota fica abaixo do ótimo de referência\n")
print(data.frame(eps = r$sensibilidade$eps,
                 excesso = signif(r$sensibilidade$excesso, 3)))
cat("   todos >= 0? ", all(r$sensibilidade$excesso >= -1e-9), "\n")

cat("\n2) Menos limiar, mais variáveis (deve ser monótono)\n")
print(r$sensibilidade[order(-r$sensibilidade$eps), c("eps", "nao_nulos")])

cat("\n3) O coordinate descent tem o menor excesso?\n")
i_cd <- which(is.na(r$sensibilidade$eps))
cat("   linha do CD é a de menor excesso? ",
    which.min(r$sensibilidade$excesso) == i_cd, "\n")

cat("\n4) A distância ao ótimo converge para zero quando eps -> 0?\n")
print(r$sensibilidade[, c("eps", "dist_ao_otimo")])
cat("   (esperado: NÃO converge; estabiliza num patamar)\n")

cat("\n5) A curva do coordinate descent é monótona decrescente?\n")
g <- subset(r$traj, rota == "descida coordenada (2010)")$gap
cat("   monótona? ", all(diff(g) <= 1e-12), " | gap final: ", signif(min(g), 3), "\n")
