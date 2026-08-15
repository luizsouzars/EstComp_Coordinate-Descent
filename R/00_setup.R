# =============================================================================
#  00_setup.R -- infraestrutura comum
#
#  Fornece tres coisas ao relatorio:
#    mostra()     exibe o codigo-fonte de uma funcao LENDO O ARQUIVO .R.
#                 E o que garante que relatorio e implementacao nao divirjam:
#                 nenhuma linha de codigo do metodo e digitada no .Rmd.
#    com_cache()  guarda o resultado de um experimento em resultados/*.rds,
#                 para que reconhecer o documento nao recalcule tudo.
#    tema grafico em tons de cinza.
# =============================================================================

RAIZ <- if (file.exists("R/01_caminho.R")) "." else
        if (file.exists("../R/01_caminho.R")) ".." else
        stop("Rode a partir da raiz do projeto v2 (pasta que contem R/).")

DIR_R   <- file.path(RAIZ, "R")
DIR_RES <- file.path(RAIZ, "resultados")
DIR_DAT <- file.path(RAIZ, "dados")
for (d in c(DIR_RES, DIR_DAT)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## Modo rapido: TRUE = poucas replicacoes (alguns minutos). FALSE = estudo completo.
if (!exists("RAPIDO")) RAPIDO <- TRUE

## --- TAMANHO DO ESTUDO -------------------------------------------------------
## Definido AQUI e em nenhum outro lugar. Tanto R/99_roda_tudo.R quanto
## relatorio.Rmd leem estes valores, de modo que os dois nao podem divergir.
PARAM <- if (RAPIDO) list(
  R_ABL = 2L,  R_CRI = 20L,  R_VEL = 10L, R_SEL = 8L,
  RHO_ABL = c(0, 0.9), NL_ABL = 20L, B_ORD = 15L
) else list(
  R_ABL = 10L, R_CRI = 100L, R_VEL = 20L, R_SEL = 50L,
  RHO_ABL = c(0, 0.5, 0.9), NL_ABL = 50L, B_ORD = 40L
)

VERSAO_PROJETO <- "v2 -- recorte: regularization path por coordinate descent"

## VERSAO_CACHE entra no nome do arquivo. Ao mudar o delineamento de um
## experimento, basta incrementar esta string: os resultados antigos deixam de
## ser encontrados e tudo e recalculado. Evita o erro classico de reportar
## numeros de uma versao anterior do codigo.
VERSAO_CACHE <- "v6"


# ---------------------------------------------------------------------------
# GUARDA CONTRA PASTA MISTURADA
# ---------------------------------------------------------------------------
# Este projeto substituiu uma versao anterior, mais larga, que implementava
# tambem GLMs, covariance updates, Huber e ressincronizacao. Os dois tem apenas
# TRES nomes de arquivo em comum, de modo que extrair um por cima do outro NAO
# apaga o antigo: sobram 01_implementacao.R, 02_extensoes.R e companhia. O R
# acaba carregando funcoes erradas -- em silencio, porque muitos nomes
# coincidem. Um erro assim custa horas para ser rastreado; paramos aqui.

.esperados <- c("00_setup.R", "00_testes.R", "01_caminho.R", "02_selecao_lambda.R",
                "03_referencia_1996.R", "04_geradores.R", "05_experimentos.R",
                "06_contribuicao.R", "07_aplicacao.R", "08_visualizacoes.R",
                "99_roda_tudo.R")
.antigos <- c("01_implementacao.R", "02_extensoes.R", "03_monte_carlo.R",
              "04_aplicacao.R", "05_contribuicao.R", "06_figuras_metodologia.R")
.presentes <- list.files(DIR_R, pattern = "[.][Rr]$")

if (length(intersect(.antigos, .presentes)))
  stop("\n  PASTA MISTURADA. Ha arquivos da versao anterior do projeto em ", DIR_R, ":\n    ",
       paste(intersect(.antigos, .presentes), collapse = ", "),
       "\n  O R vai carregar funcoes erradas. Extraia o pacote numa pasta VAZIA",
       "\n  ou apague os arquivos listados acima.", call. = FALSE)

if (length(setdiff(.esperados, .presentes)))
  stop("\n  Faltam arquivos em ", DIR_R, ":\n    ",
       paste(setdiff(.esperados, .presentes), collapse = ", "),
       "\n  Extraia o pacote inteiro numa pasta vazia.", call. = FALSE)

if (length(setdiff(.presentes, .esperados)))
  warning("Arquivos nao reconhecidos em ", DIR_R, ": ",
          paste(setdiff(.presentes, .esperados), collapse = ", "),
          ". Se forem seus, ignore.", call. = FALSE)

message("projeto: ", VERSAO_PROJETO, "  |  cache: ", VERSAO_CACHE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Instale ggplot2: install.packages('ggplot2'). Nenhum outro pacote externo e usado.")
library(ggplot2)


# ---------------------------------------------------------------------------
# Exibicao de codigo a partir do arquivo-fonte
# ---------------------------------------------------------------------------
# Convencao seguida em todos os .R deste projeto: toda funcao e definida na
# forma
#       nome <- function(...) {
#         ...
#       }
# com a chave de fechamento sozinha na coluna 0. Isso torna a extracao trivial
# e confiavel.

mostra <- function(nome, arquivo, comentarios = TRUE) {
  linhas <- readLines(file.path(DIR_R, arquivo), warn = FALSE)
  ini <- grep(paste0("^", nome, "\\s*<-\\s*function"), linhas)
  if (!length(ini)) stop("funcao nao encontrada em ", arquivo, ": ", nome)
  ini <- ini[1]
  fim <- ini
  while (fim <= length(linhas) && linhas[fim] != "}") fim <- fim + 1L
  if (comentarios) {                      # inclui o cabecalho de comentarios
    j <- ini - 1L
    while (j >= 1L && grepl("^\\s*##", linhas[j])) j <- j - 1L
    ini <- j + 1L
  }
  cat("```r\n", paste(linhas[ini:fim], collapse = "\n"), "\n```\n", sep = "")
}

## Exibe um trecho delimitado por marcadores  ## <<marca  ...  ## marca>>
mostra_trecho <- function(marca, arquivo) {
  linhas <- readLines(file.path(DIR_R, arquivo), warn = FALSE)
  a <- grep(paste0("^\\s*## <<", marca, "\\s*$"), linhas)
  b <- grep(paste0("^\\s*## ", marca, ">>\\s*$"), linhas)
  if (!length(a) || !length(b)) stop("marcador nao encontrado: ", marca)
  cat("```r\n", paste(linhas[(a[1] + 1L):(b[1] - 1L)], collapse = "\n"), "\n```\n", sep = "")
}


# ---------------------------------------------------------------------------
# Cache de experimentos
# ---------------------------------------------------------------------------
.TEMPOS <- list()

com_cache <- function(nome, expr, forcar = FALSE) {
  arq <- file.path(DIR_RES, paste0(nome, "_", VERSAO_CACHE, ".rds"))
  if (!forcar && file.exists(arq)) {
    message("[cache] ", nome)
    return(readRDS(arq))
  }
  message("[calculando] ", nome, " ...")
  t0 <- Sys.time()
  valor <- expr                       # avaliacao preguicosa: so acontece aqui
  gasto <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("[pronto]   %s em %.1f min", nome, gasto))
  ## guarda os tempos para que se possa ver DEPOIS onde o tempo foi gasto
  .TEMPOS[[nome]] <<- gasto
  saveRDS(valor, arq)
  valor
}


# ---------------------------------------------------------------------------
# Tema grafico: tons de cinza
# ---------------------------------------------------------------------------
tema_cinza <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.border     = element_rect(colour = "grey40", fill = NA),
    strip.background = element_rect(fill = "grey92", colour = "grey40"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.key       = element_blank(),
    plot.title       = element_text(face = "bold", size = 11.5),
    plot.subtitle    = element_text(colour = "grey30", size = 9.5),
    plot.caption     = element_text(colour = "grey40", size = 8, hjust = 0)
  )
theme_set(tema_cinza)

## Paleta das superficies 3D. "Grays" mantem o padrao do restante do trabalho;
## troque por "Spectral", "Zissou 1" ou "Viridis" para versoes coloridas.
PALETA <- "Grays"

cinza_cor  <- function(...) scale_colour_grey(start = 0.05, end = 0.78, ...)
cinza_fill <- function(...) scale_fill_grey(start = 0.30, end = 0.95, ...)


## Aviso de progresso dentro de experimentos longos. Sem isso nao ha como saber
## se o processo esta trabalhando ou travado -- e a diferenca importa quando um
## experimento leva dezenas de minutos.
progresso <- function(i, n, o_que = "") {
  passo <- max(1L, round(n / 20))
  if (i == 1L || i == n || i %% passo == 0L)
    message(sprintf("      %s %d/%d (%.0f%%)  %s", o_que, i, n, 100 * i / n,
                    format(Sys.time(), "%H:%M:%S")))
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Carrega a implementacao
# ---------------------------------------------------------------------------
for (f in c("01_caminho.R", "02_selecao_lambda.R", "03_referencia_1996.R",
            "04_geradores.R", "05_experimentos.R", "06_contribuicao.R",
            "07_aplicacao.R", "08_visualizacoes.R"))
  source(file.path(DIR_R, f))
