# =============================================================================
#  07_aplicacao.R -- APLICACAO A DADOS REAIS
#
#  Communities and Crime (Redmond, 2002), UCI Machine Learning Repository,
#  DOI 10.24432/C53W3X, licenca CC BY 4.0.
#  N = 1994 comunidades dos EUA; 128 atributos (5 nao preditivos, 122
#  preditivos, 1 resposta); resposta = crimes violentos por 100 mil habitantes.
#
#  POR QUE ESTE CONJUNTO
#   - nao e o do artigo (expressao genica) nem o dos tutoriais do glmnet;
#   - tem blocos de preditores quase colineares -- rendas, alugueis, imigracao
#     recente em 3/5/8/10 anos, estruturas familiares. E exatamente a situacao
#     em que a contribuicao original deste trabalho morde: se a solucao nao for
#     unica, a lista de variaveis "relevantes" pode ser artefato de ordenacao;
#   - tem faltantes reais, que obrigam a decisoes explicitas;
#   - a leitura dos coeficientes tem consequencia social, o que desloca a
#     pergunta de "qual R2" para "esta lista e confiavel".
#
#  PARA TROCAR PELOS SEUS DADOS: substitua o corpo de carrega_dados() por
#  qualquer coisa que devolva list(X = matriz numerica, y = vetor numerico).
#  Nada mais no projeto precisa mudar.
# =============================================================================

## Nomes das colunas: o arquivo communities.data nao tem cabecalho. Vieram do
## cabecalho .arff publicado na pagina do conjunto no UCI.
NOMES_CRIME <- c(
  "state","county","community","communityname","fold","population","householdsize",
  "racepctblack","racePctWhite","racePctAsian","racePctHisp","agePct12t21",
  "agePct12t29","agePct16t24","agePct65up","numbUrban","pctUrban","medIncome",
  "pctWWage","pctWFarmSelf","pctWInvInc","pctWSocSec","pctWPubAsst","pctWRetire",
  "medFamInc","perCapInc","whitePerCap","blackPerCap","indianPerCap","AsianPerCap",
  "OtherPerCap","HispPerCap","NumUnderPov","PctPopUnderPov","PctLess9thGrade",
  "PctNotHSGrad","PctBSorMore","PctUnemployed","PctEmploy","PctEmplManu",
  "PctEmplProfServ","PctOccupManu","PctOccupMgmtProf","MalePctDivorce",
  "MalePctNevMarr","FemalePctDiv","TotalPctDiv","PersPerFam","PctFam2Par",
  "PctKids2Par","PctYoungKids2Par","PctTeen2Par","PctWorkMomYoungKids","PctWorkMom",
  "NumIlleg","PctIlleg","NumImmig","PctImmigRecent","PctImmigRec5","PctImmigRec8",
  "PctImmigRec10","PctRecentImmig","PctRecImmig5","PctRecImmig8","PctRecImmig10",
  "PctSpeakEnglOnly","PctNotSpeakEnglWell","PctLargHouseFam","PctLargHouseOccup",
  "PersPerOccupHous","PersPerOwnOccHous","PersPerRentOccHous","PctPersOwnOccup",
  "PctPersDenseHous","PctHousLess3BR","MedNumBR","HousVacant","PctHousOccup",
  "PctHousOwnOcc","PctVacantBoarded","PctVacMore6Mos","MedYrHousBuilt",
  "PctHousNoPhone","PctWOFullPlumb","OwnOccLowQuart","OwnOccMedVal","OwnOccHiQuart",
  "RentLowQ","RentMedian","RentHighQ","MedRent","MedRentPctHousInc",
  "MedOwnCostPctInc","MedOwnCostPctIncNoMtg","NumInShelters","NumStreet",
  "PctForeignBorn","PctBornSameState","PctSameHouse85","PctSameCity85",
  "PctSameState85","LemasSwornFT","LemasSwFTPerPop","LemasSwFTFieldOps",
  "LemasSwFTFieldPerPop","LemasTotalReq","LemasTotReqPerPop","PolicReqPerOffic",
  "PolicPerPop","RacialMatchCommPol","PctPolicWhite","PctPolicBlack","PctPolicHisp",
  "PctPolicAsian","PctPolicMinor","OfficAssgnDrugUnits","NumKindsDrugsSeiz",
  "PolicAveOTWorked","LandArea","PopDens","PctUsePubTrans","PolicCars",
  "PolicOperBudg","LemasPctPolicOnPatr","LemasGangUnitDeploy","LemasPctOfficDrugUn",
  "PolicBudgPerPop","ViolentCrimesPerPop")


## Baixa uma unica vez e guarda em dados/. Nas execucoes seguintes le do disco.
baixa_crime <- function(destino = file.path(DIR_DAT, "communities.data")) {
  if (file.exists(destino)) return(destino)
  urls <- c(
    "https://archive.ics.uci.edu/static/public/183/communities+and+crime.zip",
    "https://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.data")
  for (u in urls) {
    ok <- tryCatch({
      if (grepl("\\.zip$", u)) {
        tmp <- tempfile(fileext = ".zip")
        utils::download.file(u, tmp, mode = "wb", quiet = TRUE)
        dentro <- utils::unzip(tmp, list = TRUE)$Name
        alvo <- grep("communities\\.data$", dentro, value = TRUE)[1]
        if (is.na(alvo)) stop("arquivo nao encontrado no zip")
        utils::unzip(tmp, files = alvo, exdir = DIR_DAT, junkpaths = TRUE)
        file.rename(file.path(DIR_DAT, basename(alvo)), destino)
      } else {
        utils::download.file(u, destino, mode = "wb", quiet = TRUE)
      }
      file.exists(destino) && file.size(destino) > 1e5
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(destino)
  }
  stop("Nao foi possivel baixar os dados. Baixe manualmente de\n",
       "  https://archive.ics.uci.edu/dataset/183/communities+and+crime\n",
       "e salve communities.data em ", DIR_DAT)
}


## Carrega e limpa. As decisoes de pre-processamento estao no corpo da funcao,
## explicitas, porque cada uma delas e uma escolha do analista.
carrega_dados <- function() {
  bruto <- utils::read.csv(baixa_crime(), header = FALSE, na.strings = "?",
                           col.names = NOMES_CRIME, stringsAsFactors = FALSE)
  resposta <- "ViolentCrimesPerPop"

  ## (i) fora os 5 atributos declaradamente nao preditivos
  d <- bruto[, setdiff(names(bruto), c("state", "county", "community",
                                       "communityname", "fold"))]

  ## (ii) fora as variaveis com muitos faltantes. As 22 da pesquisa LEMAS
  ##      faltam em 1675 das 1994 comunidades (84%): imputar 84% seria inventar
  ##      dado, entao descartamos.
  prop_na <- colMeans(is.na(d))
  descartadas <- names(d)[prop_na > 0.5]
  d <- d[, prop_na <= 0.5]

  ## (iii) imputacao pela mediana no residual (OtherPerCap tem 1 faltante)
  for (j in names(d)) if (anyNA(d[[j]]))
    d[[j]][is.na(d[[j]])] <- stats::median(d[[j]], na.rm = TRUE)

  ## (iv) fora linhas sem resposta e colunas constantes
  d <- d[!is.na(d[[resposta]]), ]
  d <- d[, !vapply(d, function(z) stats::sd(z) == 0, logical(1))]

  y <- d[[resposta]]
  X <- as.matrix(d[, setdiff(names(d), resposta)])
  list(X = X, y = y, descartadas = descartadas, n_descartadas = length(descartadas))
}


## Perfil de colinearidade: para cada variavel, sua maior correlacao absoluta
## com qualquer outra. E o diagnostico que justifica a preocupacao com
## unicidade neste conjunto.
perfil_colinearidade <- function(X) {
  C <- stats::cor(X); diag(C) <- 0
  apply(abs(C), 1, max)
}


## Analise principal. Usa SOMENTE o que foi implementado neste trabalho:
## en_path, cv_en, os criterios analiticos, viol_kkt e o diagnostico de ordem.
## nlambda_diag menor que nlambda: o diagnostico refaz o caminho B vezes, tres
## vezes seguidas. E a parte mais cara de todo o trabalho, e uma grade fina nao
## acrescenta nada a ele -- o que importa e o lambda alvo, nao a resolucao do
## caminho que leva ate la.
experimento_aplicacao <- function(dados, prop_treino = 0.7, nfolds = 10L,
                                  nlambda = 50L, nlambda_diag = 18L,
                                  B_ordem = 20L, semente = 909) {
  set.seed(semente)
  X <- dados$X; y <- dados$y
  n <- nrow(X)
  tr <- sample(n, floor(prop_treino * n))
  Xtr <- X[tr, ]; ytr <- y[tr]
  Xte <- X[-tr, ]; yte <- y[-tr]

  message("    aplicacao: validacao cruzada ", nfolds, "-fold ...")
  cv  <- cv_en(Xtr, ytr, alpha = 1, nfolds = nfolds, nlambda = nlambda)
  fit <- cv$ajuste

  ## as mesmas cinco regras de selecao do EC2, agora em dados reais
  escolhas <- c("CV (minimo)" = cv$i_min, "CV (1 EP)" = cv$i_1se,
                "AIC" = escolhe_por_criterio(fit, "AIC"),
                "BIC" = escolhe_por_criterio(fit, "BIC"),
                "GCV" = escolhe_por_criterio(fit, "GCV"))

  desempenho <- do.call(rbind, lapply(names(escolhas), function(nm) {
    k <- escolhas[[nm]]
    eta <- as.vector(fit$a0[k] + Xte %*% fit$beta[, k])
    data.frame(regra = nm, lambda = fit$lambda[k], tamanho = fit$df[k],
               eqm_teste = mean((yte - eta)^2),
               ## TRUE indica que o criterio parou na ultima casa da grade, isto
               ## e, queria um lambda menor do que o caminho ofereceu
               no_limite = k == length(fit$lambda))
  }))

  ## referencia sem penalizacao (possivel aqui porque N > p)
  mqo <- stats::lm(ytr ~ Xtr)
  cf  <- stats::coef(mqo); cf[is.na(cf)] <- 0        # posto incompleto -> 0
  desempenho <- rbind(desempenho,
    data.frame(regra = "MQO (sem penalizacao)", lambda = 0, tamanho = ncol(X),
               eqm_teste = mean((yte - cbind(1, Xte) %*% cf)^2),
               no_limite = FALSE))

  ## A CONTRIBUICAO ORIGINAL, aplicada aos dados reais.
  ##
  ## Rodamos o diagnostico com DOIS criterios de parada, e a comparacao e o
  ## ponto: `solto` reproduz o criterio absoluto usual (equivalente a
  ## max (delta beta)^2 < 1e-7, o padrao do glmnet); `apertado` usa o criterio
  ## proporcional a lambda. Se a instabilidade sob reordenacao desaparecer ao
  ## apertar a tolerancia, ela era erro numerico e nao degenerescencia -- e
  ## qualquer leitura substantiva da lista de variaveis teria sido equivocada.
  message("    aplicacao: diagnostico de reordenacao (3 x ", B_ordem, " ajustes) ...")
  diag_solto <- diagnostico_ordem(Xtr, ytr, alpha = 1, lambda_alvo = cv$lambda_min,
                                  B = B_ordem, nlambda = nlambda_diag, semente = semente,
                                  nomes = colnames(X),
                                  tol = sqrt(1e-7), escala_tol = "absoluta")
  diag1 <- diagnostico_ordem(Xtr, ytr, alpha = 1, lambda_alvo = cv$lambda_min,
                             B = B_ordem, nlambda = nlambda_diag, semente = semente,
                             nomes = colnames(X))
  diag2 <- diagnostico_ordem(Xtr, ytr, alpha = 0.5, lambda_alvo = cv$lambda_min,
                             B = B_ordem, nlambda = nlambda_diag, semente = semente,
                             nomes = colnames(X))

  list(cv = cv, ajuste = fit, desempenho = desempenho, escolhas = escolhas,
       kkt = max(viol_kkt(fit)), kkt_rel = max(viol_kkt(fit, relativa = TRUE)),
       colinearidade = perfil_colinearidade(X),
       diag_lasso = diag1, diag_en = diag2, diag_solto = diag_solto,
       n_treino = length(tr), n_teste = n - length(tr), p = ncol(X),
       nomes = colnames(X))
}
