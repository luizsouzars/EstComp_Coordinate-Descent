# Regularization path por coordinate descent

Visite o repositório no [Github](https://github.com/luizsouzars/EstComp_Coordinate-Descent).

Implementação em R, do zero, do caminho de regularização do elastic net
(Friedman, Hastie & Tibshirani, 2010, *JSS* 33(1), 1–22), com estudo de Monte
Carlo, comparação contra a rota de Tibshirani (1996) e uma contribuição
original sobre o efeito da **ordem de visita das coordenadas**.

Nenhum pacote de modelagem é usado: `glmnet` **não** é dependência. A validação
é contra teoria (KKT, monotonicidade, casos-limite), não contra outro software.

---

## 1. Requisitos

| Item | Versão | Para quê |
|---|---|---|
| R | ≥ 4.1 | roda tudo |
| `ggplot2` | ≥ 3.4 (usa `linewidth`) | gráficos — **única dependência externa obrigatória** |
| `rmarkdown` + `knitr` | qualquer | só para renderizar o `relatorio.Rmd` |

```r
install.packages(c("ggplot2", "rmarkdown", "knitr"))
```

Para renderizar o relatório em HTML também é preciso ter **Pandoc**
(já vem embutido no RStudio; fora dele, instale pelo site do Pandoc).

---

## 2. Clonar e abrir

```bash
git clone git@github.com:luizsouzars/EstComp_Coordinate-Descent.git
cd EstComp_Coordinate-Descent
```

> **Importante:** todos os scripts devem ser rodados **a partir da raiz do
> projeto** (a pasta que contém `R/`). O `00_setup.R` verifica isso e aborta com
> mensagem clara se o diretório de trabalho estiver errado. Em R:
> `setwd("caminho/para/EstComp_Coordinate-Descent")`.

---

## 3. Uso rápido

Abra o R na raiz do projeto e rode, nesta ordem:

```r
# 1. Verificação da implementação (~20 s). Só testes, nada de simulação.
source("R/00_testes.R")

# 2. Todos os experimentos em modo rápido.
source("R/99_roda_tudo.R")

# 3. Relatório final em HTML.
rmarkdown::render("relatorio.Rmd")
```

O passo 2 grava os resultados em `resultados/*.rds`. O passo 3 lê esse cache —
por isso a renderização é rápida na segunda vez.

**Atalho:** o passo 2 é opcional. Se o cache não existir, o próprio
`relatorio.Rmd` calcula tudo durante o *knit*. Rodar `99_roda_tudo.R` antes só
serve para separar "calcular" de "escrever", e para ver o resumo no console.

### Modo rápido × modo completo

O tamanho do estudo é controlado por uma única variável, `RAPIDO`, definida em
`R/00_setup.R`:

```r
# modo rápido (padrão): poucas replicações, alguns minutos
source("R/99_roda_tudo.R")

# estudo completo: muitas replicações, dezenas de minutos
RAPIDO <- FALSE
source("R/99_roda_tudo.R")
```

Os números de replicações dos dois modos estão no objeto `PARAM`, em
`R/00_setup.R` — e em nenhum outro lugar, de forma que o script e o relatório
não podem divergir.

---

## 4. Só quero usar a função principal

Não precisa rodar experimento nenhum:

```r
source("R/00_setup.R")   # carrega toda a implementação

set.seed(1)
X <- matrix(rnorm(200 * 20), 200, 20)
y <- X %*% c(3, -2, rep(0, 18)) + rnorm(200)

fit <- en_path(X, y, alpha = 1)        # alpha = 1 -> lasso; 0 -> ridge
fit$lambda[1:5]                        # grade de lambda (decrescente)
fit$beta[, 30]                         # coeficientes no 30º lambda
fit$a0[30]                             # intercepto
fit$df                                 # nº de não nulos por lambda

predict.enpath(fit, newx = X)          # predições (uma coluna por lambda)
max(viol_kkt(fit))                     # violação das condições KKT ~ 0

# escolha de lambda
cv  <- cv_en(X, y, alpha = 1, nfolds = 10)
b   <- coef_em(cv, regra = "1se")
```

Argumentos úteis de `en_path()`:

| Argumento | Padrão | Efeito |
|---|---|---|
| `alpha` | `1` | 1 = lasso, 0 = ridge, entre os dois = elastic net |
| `nlambda`, `lambda_min_ratio` | `100`, `1e-4` (`0.01` se p > N) | grade de λ |
| `pf` | `rep(1, p)` | fator de penalidade por variável (0 = nunca penalizada) |
| `padronizar` | `TRUE` | padroniza X internamente e devolve na escala original |
| `warm` | `TRUE` | **ablação:** liga/desliga a inicialização quente |
| `ativo` | `TRUE` | **ablação:** liga/desliga o conjunto ativo |
| `ordem` | `"ciclica"` | ordem de visita das coordenadas (contribuição original) |

`warm` e `ativo` **não mudam a solução**, só o custo para chegar nela — é
exatamente isso que o experimento de ablação mede.

---

## 5. O que é cada arquivo

```
R/
  00_setup.R            infraestrutura: carrega tudo, cache, tema gráfico, PARAM
  00_testes.R           verificação da implementação (~20 s)
  01_caminho.R          >>> A IMPLEMENTAÇÃO PRINCIPAL: en_path(), soft-threshold,
                            ciclo de descida coordenada, KKT
  02_selecao_lambda.R   validação cruzada K-fold + critérios analíticos (AIC/BIC/GCV)
  03_referencia_1996.R  a rota de Tibshirani (1996): ridge iterado (MM)
  04_geradores.R        cenários de simulação (AR(1), equicorrelacionado, grupos)
                            e métricas (TPR, FPR, Jaccard)
  05_experimentos.R     EC1 (ablação) e EC2 (seleção de lambda); experimento de rotas
  06_contribuicao.R     contribuição original: efeito da ordem de visita
  07_aplicacao.R        aplicação a dados reais (Communities and Crime, UCI)
  08_visualizacoes.R    geometria do problema em p = 2 (só calcula; quem desenha é o .Rmd)
  99_roda_tudo.R        roda todos os experimentos e imprime o resumo

relatorio.Rmd           o relatório (lê o cache; nenhuma linha do método é digitada aqui)
verifica_rotas.R        conferência independente do experimento das rotas
dados/                  communities.data (baixado automaticamente se faltar)
resultados/             cache dos experimentos (*.rds)
figuras/                saída de figuras
```

Convenção do projeto: o `.Rmd` **nunca** redigita código do método. A função
`mostra("en_path", "01_caminho.R")` lê o arquivo `.R` e imprime o fonte, de modo
que relatório e implementação não podem divergir.

---

## 6. Dados

A aplicação usa o **Communities and Crime** (Redmond, 2002; UCI, DOI
10.24432/C53W3X, CC BY 4.0): 1994 comunidades, 122 preditores.

O arquivo já vem em `dados/communities.data`. Se ele não existir, o projeto
baixa sozinho do UCI na primeira execução — o que exige internet. Para usar
outros dados, basta substituir o corpo de `carrega_dados()` em
`R/07_aplicacao.R` por qualquer coisa que devolva
`list(X = <matriz numérica>, y = <vetor numérico>)`. Nada mais precisa mudar.

---

## 7. Cache: como forçar o recálculo

Cada experimento passa por `com_cache()`, que grava
`resultados/<nome>_<VERSAO_CACHE>.rds` e reaproveita o arquivo se ele já existir.
No console você vê `[cache] nome` (reaproveitou) ou `[calculando] nome ...`.

Para recalcular:

```r
# tudo
unlink("resultados", recursive = TRUE)

# um experimento só
file.remove("resultados/ec1_ablacao_v6.rds")
```

Ao mudar o **delineamento** de um experimento, incremente `VERSAO_CACHE` em
`R/00_setup.R` (`"v6"` → `"v7"`): os arquivos antigos deixam de ser encontrados e
tudo é recalculado. É a proteção contra reportar números de uma versão anterior
do código.

---

## 8. Verificação independente

Além de `R/00_testes.R`, o script `verifica_rotas.R` recalcula o experimento das
rotas do zero (sem tocar no cache) e confere as invariantes na tela:

```r
source("verifica_rotas.R")
```

---

## 9. Problemas comuns

| Sintoma | Causa / solução |
|---|---|
| `Rode a partir da raiz do projeto v2 (pasta que contem R/).` | diretório de trabalho errado — use `setwd()` para a raiz do repositório |
| `PASTA MISTURADA. Ha arquivos da versao anterior...` | há `.R` de uma versão antiga em `R/`. Extraia o projeto numa pasta **vazia** ou apague os arquivos listados na mensagem |
| `Faltam arquivos em R/` | clone/extração incompleta |
| `Instale ggplot2: install.packages('ggplot2')` | falta o único pacote externo obrigatório |
| `linewidth` / erro de tema no ggplot2 | `ggplot2` anterior à 3.4 — atualize |
| Falha ao baixar os dados | sem internet: baixe `communities.data` do UCI e coloque em `dados/` |
| `99_roda_tudo.R` demorando demais | você está em `RAPIDO <- FALSE`. Rode `RAPIDO <- TRUE` (ou reinicie a sessão) |
| Números do relatório não batem com o código | cache velho — apague `resultados/` ou incremente `VERSAO_CACHE` |

---

## 10. Referências

- Friedman, J., Hastie, T. & Tibshirani, R. (2010). *Regularization Paths for
  Generalized Linear Models via Coordinate Descent.* **JSS** 33(1), 1–22.
- Tibshirani, R. (1996). *Regression Shrinkage and Selection via the Lasso.*
  **JRSS-B** 58(1), 267–288.
- Tibshirani, R. J. (2013). *The Lasso Problem and Uniqueness.*
  **Electronic Journal of Statistics** 7, 1456–1490.
- Redmond, M. (2002). *Communities and Crime.* UCI Machine Learning Repository.
  DOI 10.24432/C53W3X.
