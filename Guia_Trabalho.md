# Guia do trabalho

**O regularization path por coordinate descent** — implementação, medição e extensão de
Friedman, Hastie & Tibshirani (2010), *Regularization Paths for Generalized Linear
Models via Coordinate Descent*, JSS 33(1), 1–22.

Este arquivo é o guia de leitura do `relatorio.html`, o material de estudo para a
apresentação, e o lugar onde as decisões do trabalho estão justificadas. Ele não
contém código: para isso, o próprio relatório exibe as funções lendo os arquivos
`.R`.

> **Nota sobre as referências.** A numeração das seções de **2010** foi conferida
> diretamente no PDF. A de **1996** vem de conhecimento consolidado, mas o PDF
> disponível é um *scan* sem camada de texto e não pôde ser verificado por
> máquina — os pontos marcados com ⚠️ devem ser confirmados na sua cópia antes da
> apresentação.

---

## Sumário

1. [A tese em uma frase](#1-a-tese-em-uma-frase)
2. [O problema, e por que ele é computacional](#2-o-problema-e-por-que-ele-é-computacional)
3. [O recorte: o que implementamos e por quê](#3-o-recorte-o-que-implementamos-e-por-quê)
4. [Mapa: cada seção do HTML contra os artigos](#4-mapa-cada-seção-do-html-contra-os-artigos)
5. [Os conceitos, com intuição](#5-os-conceitos-com-intuição)
6. [Os experimentos: pergunta, desenho, achado](#6-os-experimentos-pergunta-desenho-achado)
7. [A contribuição original](#7-a-contribuição-original)
8. [Os dois erros que encontramos no caminho](#8-os-dois-erros-que-encontramos-no-caminho)
9. [Números para levar decorados](#9-números-para-levar-decorados)
10. [Como verificar cada resultado](#10-como-verificar-cada-resultado)
11. [Perguntas prováveis da banca](#11-perguntas-prováveis-da-banca)
12. [Fronteiras e limitações](#12-fronteiras-e-limitações)
13. [Glossário](#13-glossário)
14. [Roteiro de apresentação](#14-roteiro-de-apresentação)

---

## 1. A tese em uma frase

> **1996 entregou um estimador. 2010 entregou um algoritmo — e o algoritmo
> redefiniu o que o estimador podia ser usado para fazer.**

O lasso de 2010 é **exatamente** o de 1996: mesma solução, mesmas propriedades
estatísticas. Não há estimador novo. O que mudou foi o que se consegue calcular,
em que escala, com que confiabilidade — e, como consequência, **qual prática
estatística passou a ser recomendada**.

O trabalho não apenas narra essa transição: ela vira hipótese testável no
experimento EC2.

---

## 2. O problema, e por que ele é computacional

Regressão com muitos preditores impõe dois problemas ao mesmo tempo:

- **acurácia preditiva** — mínimos quadrados tem viés baixo e variância alta;
- **interpretação** — queremos um subconjunto, não 100 coeficientes.

O lasso resolve os dois de uma vez. Essa é a contribuição de 1996. Mas *resolver
o lasso na prática*, para tamanhos de problema que interessam, é uma questão
separada — e é a contribuição de 2010.

> **Intuição — por que isso muda o modo de avaliar o artigo.** Um artigo de
> **estimador** é julgado por viés, variância e consistência. Um artigo de
> **algoritmo** precisa ser julgado por **correção, custo e robustez**. É por
> isso que o estudo de Monte Carlo deste trabalho não mede viés de coeficiente:
> mede violação das condições de otimalidade, visitas de coordenada e
> estabilidade de seleção.

---

## 3. O recorte: o que implementamos e por quê

O artigo de 2010 tem sete componentes. Implementamos **um**:

> **O *regularization path* percorrido por *coordinate descent* cíclico com
> *warm starts*** — a função `en_path()`, no caso gaussiano.

### As três razões (e é isso que a banca vai querer ouvir)

**1. Não é o *coordinate descent* que o artigo reivindica como novidade.**
Os próprios autores o atribuem a Fu (1998), Daubechies et al. (2004), Shevade &
Keerthi (2003), Van der Kooij (2007) e Wu & Lange (2008), e observam que o método
foi redescoberto várias vezes sem que seu poder fosse percebido. O que
reivindicam é tê-lo aplicado ao **caminho inteiro**. O precursor direto deles
chama-se, literalmente, *Pathwise Coordinate Optimization* (Friedman et al.,
2007), e o título de 2010 põe **Paths** antes de *Coordinate Descent*.

**2. É a peça que faltava em 1996.** O gargalo do lasso original não era a
qualidade do estimador. Era que **cada valor do parâmetro de regularização exigia
uma execução completa e independente de um otimizador caro** — nada era herdado
de uma solução para a próxima. Escolher o parâmetro por validação cruzada
multiplicava esse custo pelo número de *folds*.

**3. As demais partes dependem dela.** As *covariance updates* (2.2) são uma
otimização interna ao caminho. A extensão a GLMs (Seção 3) funciona porque os
autores **não implementam verificação de divergência do Newton** — o argumento é
que os *warm starts* mantêm cada problema tão próximo do anterior que a
aproximação quadrática nunca chega a quebrar. Retirar o caminho derruba as
extensões; retirar as extensões não derruba o caminho.

### O que ficou fora, e uma observação que reforça a coerência

Não implementamos: *covariance updates* (2.2), *sparse updates* (2.3),
*weighted updates* (2.4), regressão logística (3), multinomial (4), resposta
agrupada (4.2).

> **Repare:** a Seção 2.4 (atualizações ponderadas) é a base das Seções 3 e 4.
> Cortando 2.4 caem 3 e 4 junto. Não cortamos pedaços avulsos — cortamos um
> bloco inteiro e conectado do artigo.

---

## 4. Mapa: cada seção do HTML contra os artigos

### 4.1. O que é implementação do artigo de 2010

| Seção do HTML | Artigo de 2010 | O que faz |
| :-- | :-- | :-- || **O algoritmo** (caixa com o pseudocódigo) | 2.1 + 2.5 + 2.6 reunidas | Enuncia o método completo num lugar só |
| **A implementação → O passo elementar** | **2.1 Naive updates** | O operador de *soft-thresholding* e a atualização coordenada |
| **A implementação → Um ciclo** | **2.1 Naive updates** | Uma varredura de $j=1,\dots,p$, com atualização do resíduo em $O(N)$ só quando o coeficiente muda |
| **A implementação → Convergência para um λ** | **2.6 Other details** | A estratégia de *active set* |
| **A implementação → O caminho** | **2.5 Pathwise coordinate descent** + **2.6** | *Warm starts*, $\lambda_{\max}$ em forma fechada, grade log-espaçada, padronização, intercepto, fatores de penalidade |

Essas quatro funções **são** o artigo. Todo o resto do documento gira em torno delas.

### 4.2. Uma parte do artigo que aparece só pelos efeitos

A **Seção 6 (Selecting the tuning parameters)** está implementada em
`R/02_selecao_lambda.R` — validação cruzada $K$-fold e a regra de um erro-padrão
— mas o **código não é exibido** no HTML. Ela aparece por:

- a figura da curva de validação cruzada, na aplicação;
- as linhas `CV (mínimo)` e `CV (1 EP)` nas tabelas;
- as linhas tracejada e pontilhada nos gráficos de caminho.

Se perguntarem "onde está a Seção 6?", a resposta é: existe, roda, e optamos por
não imprimi-la no relatório.

### 4.3. O que é do artigo de 1996

| Seção do HTML | Artigo de 1996 | O que faz |
| :-- | :-- | :-- |
| **1996 → 2010: O lasso como Tibshirani o definiu** | Definição do lasso na forma **restrita**, $\sum_j\lvert\beta_j\rvert\le t$ | O ponto de partida |
| **Validação → caso ortonormal** | ⚠️ eq. (3), caso de desenho ortonormal | A fórmula fechada que nosso `S` deve reproduzir à precisão de máquina |
| **O contraste com 1996, medido** | ⚠️ Seção 4, algoritmos; a aproximação por ridge iterado | A rota de 1996 implementada, para contraste medido e não afirmado |
| **Monte Carlo II — a hipótese** | ⚠️ Seção 3, *prediction error and estimation of t*: GCV e estimador de risco tipo Stein | Os critérios analíticos que existiam **porque** a CV era cara |
| **Geradores** (`gera_ar1`) | ⚠️ Exemplo 1 das simulações: desenho AR(1) | O cenário de referência |

> ⚠️ Confirme os números de seção e equação na sua cópia do PDF de 1996. O
> conteúdo está certo; a numeração é o que não pude verificar.

### 4.4. O que está no HTML e **não** é de nenhum dos dois artigos

Esta é a fronteira que a banca testa. Tenha na ponta da língua:

| No HTML | De onde vem |
| :-- | :-- |
| **A auditoria de otimalidade** (`viol_kkt`) | Nossa. Baseada nas condições de otimalidade padrão de otimização convexa. O artigo de 2010 **não** propõe verificação de KKT — quem faz isso é o artigo das *strong rules*, de 2012 |
| **A geometria do problema** (superfícies 3D, corte 1D, zigue-zague, *warm start*) | Nossas. **Mas** a trajetória desenhada é *coordinate descent* de verdade, a mesma atualização da Seção 2.1 rodando com $p=2$ — não é esquema ilustrativo |
| **AIC, BIC e GCV** no Monte Carlo II | Não estão em 2010. O GCV vem de 1996; AIC e BIC, da literatura geral |
| **Os interruptores da ablação** (`warm`, `ativo`) | O delineamento é nosso; os ingredientes ligados e desligados são do artigo |
| **Visitas de coordenada** como unidade de custo | Nossa |
| **O critério de parada proporcional a λ** | Nosso, e é uma correção a uma falha real do critério usual |
| **Contribuição original: ordem de visita** | Nossa, com base teórica em R. J. Tibshirani (2013) |
| **O diagnóstico de reordenação** | Nosso |

---

## 5. Os conceitos, com intuição

### 5.1. Restrita contra penalizada: a palavra é *separável*

1996 define o lasso por uma **restrição**: $\sum_j\lvert\beta_j\rvert\le t$. Ela
equivale a $2^p$ desigualdades lineares, uma para cada combinação de sinais.

2010 escreve o mesmo problema na forma **penalizada**, com $\lambda P_\alpha(\beta)$.
Por dualidade lagrangiana os dois são equivalentes: para cada $t$ existe um
$\lambda$ que dá a mesma solução. **Estatisticamente nada muda. Estruturalmente
muda tudo.**

> **Intuição.** Numa restrição $\sum_j\lvert\beta_j\rvert\le t$, mexer em
> $\beta_3$ consome orçamento dos outros: as coordenadas estão amarradas e não se
> pode otimizar uma sem renegociar todas. Na forma penalizada cada coordenada
> paga o **próprio pedágio** $\lambda\lvert\beta_j\rvert$, independentemente das
> demais. Fixados os outros coeficientes, sobra um problema de uma variável só —
> e esse tem solução analítica.

A cadeia inteira sai daí:

```text
separável → forma fechada por coordenada → ciclo barato → caminho barato → warm start
```

E o *warm start* devolve o favor: cada problema do caminho é fácil **porque** a
solução anterior já é quase a resposta.

### 5.2. De onde vêm os zeros

Esta é a intuição mais importante do trabalho, e a figura de três painéis existe
para ela.

O termo de erro quadrático é uma **tigela elíptica**, lisa em toda parte. O termo
$\lambda\sum_j\lvert\beta_j\rvert$ acrescenta um **cone com quinas** ao longo dos
eixos — é onde $\lvert\beta_j\rvert$ deixa de ser derivável. A soma herda os
vincos.

> **Intuição.** Num ponto liso o mínimo só pode estar onde a derivada é zero, e a
> chance de isso cair exatamente em $\beta_j=0$ é nula. Numa **quina** não existe
> *uma* derivada, existe um **intervalo** de subgradientes: basta que o gradiente
> do termo quadrático caiba dentro desse intervalo para o mínimo ficar preso ali
> — e ele fica preso para uma **faixa inteira** de valores de $\lambda$.
>
> **É a quina, e não a inclinação da penalidade, que produz zeros exatos.**
> O ridge encolhe porque inclina a tigela; o lasso zera porque a amassa.

O corte unidimensional mostra isso com mais nitidez ainda: o objetivo em função
de um coeficiente tem um **bico** em zero, onde a inclinação salta de $-\lambda$
para $+\lambda$. Enquanto a inclinação do termo quadrático em zero for menor que
$\lambda$ em módulo, **nenhuma direção desce** — e o mínimo é exatamente zero.
Essa desigualdade é literalmente a condição KKT, e é o que o operador
$S(\cdot,\lambda\alpha)$ implementa.

### 5.3. Por que o caminho pode começar

$\lambda_{\max}$ é o menor $\lambda$ que zera **todos** os coeficientes
penalizados. Acima dele o *soft-thresholding* zera toda coordenada.

> **Intuição.** O caminho parte de um ponto **exato e conhecido de antemão**
> ($\beta=0$), não de um chute. É isso que faz o *warm start* funcionar já na
> primeira iteração — não existe fase de aquecimento.

### 5.4. Warm start

> **Intuição.** As soluções ao longo do caminho formam uma **curva contínua**:
> quando $\lambda$ cai um pouco, a solução se move pouco. Partir da solução
> anterior significa começar a poucos passos do destino; partir da origem
> significa atravessar o espaço inteiro toda vez.

A figura do *warm start* mostra um detalhe que só ela revela: as soluções ficam
**sobre os eixos** em boa parte do caminho — isto é, o *active set* é pequeno —
e é por isso que o segundo ingrediente tem o que economizar.

### 5.5. Active set

> **Intuição.** Numa solução com 12 variáveis ativas entre 1000, 988 coordenadas
> por ciclo são apenas **consultadas**: o produto interno é calculado, o
> coeficiente não muda, e o resíduo não é tocado. A atualização do resíduo, que
> custa $O(N)$, só acontece quando há mudança. **A esparsidade que o método
> produz é também o que o torna rápido.**

Detalhe de correção que vale citar: a convergência só é declarada após um ciclo
**completo** sem mudanças. Parar depois de convergir apenas dentro do *active
set* seria um erro sutil — o modelo pareceria convergido e estaria errado, porque
uma variável de fora poderia ter direito de entrar.

### 5.6. Por que a correlação atrapalha

> **Intuição.** Sem correlação as curvas de nível são quase circulares e dois
> passos bastam: ajusta $\beta_1$, ajusta $\beta_2$, chegou. Com correlação alta
> as curvas viram elipses finas e **inclinadas** — e a direção de descida útil é
> a diagonal, que o algoritmo **não pode tomar**, porque só anda paralelamente
> aos eixos. Ele escorrega pela parede do vale em zigue-zague, dando muitos
> passos curtos.

Essa é a razão teórica para esperar que a correlação encareça o método. Os
autores dizem ter esperado isso e **não observado**. O EC1 resolve a
controvérsia, e a resposta é interessante (§6.1).

### 5.7. KKT: a diferença entre "parou" e "está certo"

As condições de Karush–Kuhn–Tucker dizem, para o lasso: para coeficientes não
nulos, o gradiente do termo quadrático tem de igualar $\lambda\alpha$ com o sinal
do coeficiente; para coeficientes nulos, tem de ficar **dentro** de
$[-\lambda\alpha, +\lambda\alpha]$.

> **Intuição.** O critério de parada do algoritmo pergunta "os coeficientes
> pararam de mudar?". As KKT perguntam "esta é a resposta certa?". **São
> perguntas diferentes**, e é perfeitamente possível que a primeira responda sim
> e a segunda não.

Por que isso importa metodologicamente: é o que permite afirmar que a
implementação está correta **sem comparar com o `glmnet`**. Comparar duas
implementações não diz qual das duas está certa. Verificar as condições de
otimalidade diz.

### 5.8. A tolerância e a equivariância de escala (nossa correção)

Duas armadilhas, ambas encontradas rodando o trabalho:

**Primeira.** O critério olha o **passo**, não a **distância ao ótimo**. Sob
correlação alta o método converge linearmente com razão próxima de 1, e a
distância ao ótimo é da ordem do passo dividido por $(1-\text{razão})$. Passo
pequeno não é sinônimo de estar perto.

**Segunda, e mais sutil.** Um limiar **absoluto** não é equivariante em escala.

> **Intuição.** As KKT comparam o gradiente com $\lambda\alpha$, e $\lambda$
> percorre quatro ordens de grandeza ao longo do caminho. Um erro de $10^{-5}$ é
> **desprezível** quando $\lambda\alpha=1$ e **catastrófico** quando
> $\lambda\alpha=10^{-5}$. O mesmo número absoluto significa coisas opostas em
> pontos diferentes do mesmo caminho.

Nossa correção: exigir $\lvert\Delta\beta_j\rvert < \text{tol}\cdot\lambda_k$. O
critério fica proporcional à escala local, e a violação KKT relativa passa a ser
aproximadamente constante ao longo do caminho.

### 5.9. Não unicidade

> **Intuição.** Com duas colunas **idênticas**, o ajuste só enxerga a soma
> $\beta_1+\beta_2$ — e a penalidade $\ell_1$ também só enxerga a soma, desde que
> os sinais sejam iguais. Duas funções que dependem apenas da soma têm o **mesmo
> valor** ao longo de toda a reta $\beta_1+\beta_2=c$: o fundo do vale é
> perfeitamente plano. O conjunto de ótimos é um **segmento**, não um ponto.
>
> O termo $\ell_2$ quebra o empate porque $\beta_1^2+\beta_2^2$ **não** depende
> só da soma — é mínimo quando os dois são iguais. Daí a divisão igualitária.

O que muda entre as ordens de visita é a **atribuição**, não a soma. O modelo
prediz igual e conta uma história diferente.

---

## 6. Os experimentos: pergunta, desenho, achado

### 6.1. Monte Carlo I — ablação dos ingredientes

**A pergunta.** Os dois ingredientes que fazem do caminho algo barato — *warm
start* e *active set* — funcionam? Quanto cada um economiza?

**Por que ablação em vez de comparar com outro software.** Comparar nossa
implementação com o `glmnet` confundiria o efeito do algoritmo com diferenças de
linguagem, compilador e estilo de código. Ligando e desligando os ingredientes
**dentro da mesma implementação**, com a mesma amostra e a mesma grade de
$\lambda$, sobra só o efeito que interessa.

**O desenho.** Pares $(N,p)$ explícitos cobrindo os dois regimes ($N>p$ duas
vezes, $p>N$ duas vezes) $\times$ $\rho$ $\times$ 4 configurações, comparação
pareada.

**Primeiro achado: as quatro configurações resolvem o mesmo problema.** Isso é a
**premissa** de todo o resto — se não valesse, comparar custos seria sem sentido.

> **Por que duas medidas de concordância, e não uma.** Comparar os vetores
> $\beta$ pode enganar: quando a solução não é única, dois pontos distintos podem
> ter exatamente o mesmo valor do objetivo. A diferença entre os vetores mede
> degenerescência do **problema**; a diferença entre os objetivos mede erro do
> **algoritmo**. Reportamos as duas.

**Segundo achado: o *warm start* domina.** Desligá-lo custa muito mais que
desligar o *active set*.

**Terceiro achado, e o mais interessante: há interação, não multiplicação.**
Num dos cenários, desligar o *warm start* **mantendo** o *active set* custa
**mais** do que desligar os dois.

> **Intuição.** Numa partida fria o conjunto ativo começa vazio e cresce aos
> poucos. O laço interno converge com esmero sobre um conjunto que ainda está
> **errado**, e todo esse trabalho é refeito quando o ciclo completo seguinte
> admite variáveis novas. **O *active set* isoladamente pode ser prejuízo — ele
> só paga quando o *warm start* já entrega um conjunto quase certo.**

**Quarto achado: a correlação encarece, e as duas posições se reconciliam.**
A teoria de otimização está certa sobre o *coordinate descent* cíclico; os
autores estão certos sobre o **algoritmo deles**, porque o *warm start* absorve
boa parte do problema.

### 6.2. Monte Carlo II — o que o caminho barato permitiu

**A hipótese, que é a tese do trabalho virando teste.** Em 1996, critérios
analíticos (GCV, risco tipo Stein) eram atraentes **porque a validação cruzada
era cara**: $K$ *folds* sobre uma grade de $G$ valores custavam $K \times G$
execuções completas e independentes do otimizador. Depois de 2010, a mesma
validação cruzada custa $K+1$ **caminhos**, e um caminho inteiro custa da ordem de
um único ajuste.

> **A comparação que vale levar ao slide.** Com 100 valores de $\lambda$ e 10
> *folds*: mil execuções do otimizador em 1996, onze caminhos depois de 2010. **A
> grade cresceu de dezenas para centenas de valores e o custo caiu.** É isso que
> explica por que os critérios analíticos saíram de cena — e é uma afirmação
> segura, porque não depende de nenhum número específico citado do artigo.

Os critérios analíticos deixaram silenciosamente de ser recomendados. **A troca
compensou?**

**O desenho.** Cinco regras — CV no mínimo, regra de 1 EP, AIC, BIC, GCV —
aplicadas ao **mesmo caminho**, na **mesma réplica**. Comparação pareada: as
diferenças não carregam variabilidade amostral entre métodos.

**Uma decisão de desenho que vale explicar.** Os seis efeitos verdadeiros estão
em **escada**, de 2 a 0,25.

> Na primeira versão do experimento eram todos fortes, e o resultado foi que
> **todos os critérios acertaram todos**: a taxa de verdadeiros positivos ficou
> presa em 1 nos quinze cenários, sem poder de discriminação nenhum. Um
> delineamento em que a tarefa é fácil demais não compara nada.

**Duas conclusões que vão além de comparar médias.**

*A CV é a única regra definida em todo o caminho.* AIC, BIC e GCV usam graus de
liberdade estimados por $|A|$ e deixam de fazer sentido
quando $\lvertA\rvert\ge N$ — corriqueiro no início do caminho quando
$p>N$. O GCV é o caso mais visível: o fator $(1-\text{gl}/N)^{-2}$ explode.
**Quando a viabilidade computacional mudou, o método que ficou barato também era
o mais geral.**

*A justificativa teórica veio depois do uso.* Usar $|A|$ como
graus de liberdade do lasso só foi justificado por Zou, Hastie e Tibshirani
(2007) — **depois** de o critério já estar em uso.

**A coluna "proporção no limite da grade"** merece atenção: uma proporção alta
não é uma escolha, é uma **reclamação**. O critério queria ir além do que o
caminho ofereceu, e o que se reporta ali é a fronteira da grade, não o ótimo do
critério.

### 6.3. Aplicação — Communities and Crime

**Os dados.** UCI, DOI 10.24432/C53W3X, CC BY 4.0. 1994 comunidades dos EUA,
~100 preditores após limpeza, resposta = crimes violentos por 100 mil habitantes.

**Por que este conjunto, e não o do artigo.** Três razões, nenhuma delas
conveniência:

1. tem a patologia que motiva a discussão — blocos de preditores quase
   colineares (rendas, aluguéis, imigração recente em 3/5/8/10 anos, estruturas
   familiares);
2. tem **faltantes reais**: 22 variáveis da pesquisa LEMAS faltam em 1675 das
   1994 comunidades (84%). Imputar 84% seria inventar dado — descartamos, e a
   decisão está explícita no código;
3. a leitura dos coeficientes tem **consequência social**, o que desloca a
   pergunta de "qual $R^2$" para "esta lista de variáveis é confiável".

**Uma leitura honesta que convém não maquiar.** Com $N/p \approx 14$, a
penalização compra **pouco** em predição — fica na mesma casa do MQO. O que ela
compra é parcimônia (14 variáveis contra 100, na regra de 1 EP). Em conjuntos
assim o argumento a favor da regularização é **interpretativo e de estabilidade**,
não preditivo. Dizer isso é mais forte do que fingir um ganho que não existe.

---

## 7. A contribuição original

### 7.1. A lacuna

O artigo varre as coordenadas na ordem $j=1,\dots,p$ e **nunca justifica essa
escolha**. Ela aparece como detalhe de implementação.

### 7.2. O que a teoria diz, e o que os dados disseram

Com colunas **exatamente** linearmente dependentes, o conjunto de soluções ótimas
é uma **face** do poliedro, não um ponto — e o *coordinate descent* devolve um
ponto dessa face, escolhido por quem foi visitado primeiro.

Quando as colunas são apenas **muito correlacionadas**, mas distintas, a solução
é única com probabilidade 1 (R. J. Tibshirani, 2013).

**O desenho, em uma frase, para quando perguntarem como foi medido.** $N=100$,
$p=30$ — três grupos de cinco colunas correlacionadas mais quinze de ruído.
Cruzamento de $\rho_g\in\{0{,}90;\,0{,}99;\,1\}$ contra
$\alpha\in\{1;\,0{,}95;\,0{,}80;\,0{,}50\}$, com $\rho_g=1$ significando colunas
idênticas. Para cada célula, o caminho é reajustado $B$ vezes com ordens de
visita sorteadas e comparado num único $\lambda=0{,}05\,\lambda_{\max}$. **Os
dados, o $\lambda$ e o $\alpha$ são sempre os mesmos; só a ordem muda.**

> **O achado, e ele contrariou a expectativa inicial.** A primeira versão do
> experimento usava correlação intragrupo 0,99 e encontrou Jaccard **igual a 1**:
> nenhuma instabilidade. Só com colunas **idênticas** o fenômeno aparece.
> **Correlação alta não gera degenerescência; dependência exata gera.**

**Não pare no Jaccard — e este é o ponto que a banca pode cobrar.** Conjunto que
varia sob reordenação é compatível com **duas** explicações: degenerescência
verdadeira ou tolerância frouxa (Seção 7.3). O que separa uma da outra é a
**amplitude do objetivo**, e na simulação ela é nula. Conjuntos diferentes com o
mesmo valor da função objetivo é a assinatura de não unicidade. Se fosse
tolerância, o objetivo variaria junto.

Se você tiver só uma frase para a contribuição, use esta: *os reajustes chegam a
listas diferentes de variáveis com exatamente o mesmo valor do objetivo.*

Isso não é curiosidade acadêmica. Dependência exata é comum em dados aplicados —
indicadoras redundantes, variáveis derivadas de outras, totais que são soma de
partes, e **todo** caso com $p>N$ — e é justamente onde ninguém procura, porque a
saída do algoritmo parece perfeitamente determinada.

**E o remédio é barato, com um detalhe que vale citar.** Sob dependência exata,
$\alpha=0{,}95$ já devolve Jaccard 1, e o erro de predição sobe pouco. Sob
correlação 0,99, o mesmo $\alpha=0{,}95$ **reduz** o erro em relação ao lasso. O
custo do remédio muda de sinal conforme o cenário, e no cenário mais comum na
prática ele é negativo. Não é preciso abandonar a parcimônia do lasso: basta não
estar **exatamente** em $\alpha=1$.

**Duas ressalvas que você mesmo deve levantar.** A grade salta de $\rho_g=0{,}99$
para 1 e de $\alpha=0{,}95$ para 1: sabemos que há descontinuidade em cada
intervalo, não onde está. E **seleção estável não é seleção correta** — a
recuperação do grupo verdadeiro é pior em $\rho_g=0{,}99$, onde o Jaccard é 1, do
que sob dependência exata. Em 0,99 o algoritmo elege uma variável de forma
perfeitamente reprodutível, inclusive quando elege a errada.

### 7.3. O segundo mecanismo, que se confunde com o primeiro

Existe outra fonte de variação sob reordenação, de natureza completamente
diferente: com tolerância finita, coeficientes numericamente minúsculos **oscilam
entre zero e não zero** conforme a ordem de visita.

> Para quem lê a saída, isso é **indistinguível** de degenerescência verdadeira —
> e o remédio é o oposto: apertar a tolerância, não trocar de estimador.

**Ordem de exposição, se for apresentar.** A tabela de decisão abaixo já foi
usada na Seção 7.2, mesmo aparecendo formalmente aqui. Vale antecipá-la no slide:
mostre os dois sinais **antes** de mostrar a simulação, e a conclusão de 7.2
deixa de parecer um salto. O trabalho aplica o mesmo par de sinais duas vezes e
chega a conclusões **opostas** — degenerescência na simulação, tolerância nos
dados reais. Esse contraste é o argumento mais forte de que o diagnóstico
discrimina de fato, em vez de sempre dizer a mesma coisa.

### 7.4. A proposta

Reajustar $B$ vezes com ordens de visita sorteadas e olhar **duas** coisas:

| Conjunto selecionado | Valor do objetivo | Diagnóstico |
| :-- | :-- | :-- |
| varia | constante | **não unicidade verdadeira** — os pontos são igualmente ótimos |
| varia | varia | **convergência incompleta** — aperte a tolerância antes de concluir |
| constante | constante | seleção estável |

> **Por que o *bootstrap* não serve para isso.** O *bootstrap* muda os dados, e
> portanto mistura degenerescência, erro numérico e incerteza amostral numa
> medida só. Aqui os dados e o $\lambda$ são **sempre os mesmos** — a única coisa
> que muda é a ordem de visita. Qualquer variação observada tem só duas causas
> possíveis, e o valor do objetivo diz qual.

O **limiar de magnitude** cumpre o papel simétrico: sem ele, um coeficiente de
$10^{-12}$ conta como variável do modelo, e o diagnóstico mede arredondamento em
vez de estrutura.

### 7.5. O que é original, e o que não é

**Não é original:** a existência de soluções não únicas no lasso é conhecida, e a
*elastic net* foi proposta em parte por causa disso.

**É original:** (a) mostrar que a **ordem de varredura** — apresentada como
detalhe de implementação — é o mecanismo concreto pelo qual a não unicidade se
manifesta na saída; (b) mostrar que o fenômeno é de **dependência exata**, não de
correlação alta; (c) transformar isso num **diagnóstico executável** que separa
degenerescência de erro numérico.

### 7.6. O desfecho na aplicação, que é o melhor argumento do trabalho

O mesmo diagnóstico rodou nos dados reais com **dois critérios de parada**, tudo
o mais idêntico. Com o critério absoluto usual, aparecem variáveis
intercambiáveis e o Jaccard cai. Com o critério proporcional a $\lambda$, a
instabilidade **desaparece**.

> **A lição.** Se tivéssemos parado na primeira execução, teríamos escrito que a
> lista de fatores associados à criminalidade violenta é instável e que parte
> dela é artefato de ordenação — uma afirmação substantiva, publicável, e
> **errada**. O que havia era tolerância frouxa.
>
> É por isso que o diagnóstico precisa das duas peças: **reordenar expõe o
> sintoma; o valor do objetivo e o aperto da tolerância dizem qual é a doença.**

**Recomendação prática, para levar ao slide de impacto.** Ao reportar variáveis
selecionadas por lasso em dados colineares: (i) rode com ordens de visita
sorteadas e publique a frequência de seleção; (ii) antes de interpretar qualquer
instabilidade, repita com a tolerância apertada. O primeiro passo custa $B$
ajustes; o segundo custa um.

---

## 8. Os dois erros que encontramos no caminho

Material de primeira para a declaração de uso de IA e para a arguição — mostram
que o trabalho foi de fato auditado.

**1. O critério de parada não era equivariante em escala.** Na segunda execução,
com limiar absoluto, a violação KKT **relativa** chegou a 4,6 no fim do caminho —
o gradiente estava mais distante da condição de otimalidade do que o próprio
$\lambda$. O diagnóstico só apareceu **porque medimos a violação relativa**; em
valor absoluto ela era $2{,}7\times10^{-4}$ e passaria por aceitável.

**2. O experimento da contribuição testava a hipótese errada.** Media correlação
0,99 e não achava nada, porque ali a solução é única com probabilidade 1. Faltava
o caso de dependência exata. O resultado nulo virou o achado, depois de entender
*por que* era nulo.

**3. Um erro de desenho que a auditoria pegou.** Ao encolher o EC1 para ganhar
tempo, cruzei $N\in\{120,320\}$ com $p\in\{30,90\}$ — e todas as quatro células
ficaram com $N>p$. O cenário $p>N$ desapareceu da tabela **sem aviso**. A
correção foi trocar o cruzamento por pares explícitos.

---

## 9. Números para levar decorados

Separe os que são **estruturais** (não mudam de execução para execução) dos que
são **numéricos** (mudam com semente e delineamento).

### Estruturais — pode afirmar sem consultar

- O artigo de 2010 tem **sete** componentes; implementamos **um**.
- Em 1996, validação cruzada de $K$ *folds* sobre grade de $G$ valores custava
  $K \times G$ **execuções completas e independentes** do otimizador, mais o
  ajuste final. Depois de 2010: $K+1$ **caminhos**. ⚠️ Não cite "75 execuções"
  como se estivesse no artigo — é aritmética ilustrativa ($5\times15$), e a grade
  que Tibshirani usa não foi verificada.
- A restrição $\ell_1$ equivale a $2^p$ desigualdades lineares.
- Em 1996, na prática: entre $0{,}5p$ e $0{,}75p$ problemas de mínimos quadrados restritos.
- LEMAS: 22 variáveis faltando em **1675 de 1994** comunidades (84%).
- Communities and Crime: **1994** observações, ~**100** preditores após limpeza.
- Com colunas idênticas, a **soma** dos coeficientes é estável; a **divisão** não.
- Custo por ciclo: $O(pN)$ na versão *naive*.

### Numéricos — confira na execução que você vai apresentar

- Violação KKT relativa nas validações (deve ficar na ordem de `tol`).
- Custo relativo de cada configuração na ablação.
- Fator de encarecimento com a correlação.
- Jaccard sob reordenação: **1,000** para correlação 0,90 e 0,99; **abaixo de 1**
  só com colunas idênticas.
- Amplitude do objetivo entre reordenações: **nula** no caso degenerado — é o que
  prova que é não unicidade e não erro numérico. Sem este número, o Jaccard
  sozinho não conclui nada.
- Nº médio de variáveis **intercambiáveis** no caso degenerado (entram em algumas
  ordens de visita e não em outras).
- Erro contra $E(y\mid x)$ em $\alpha=1$ **contra** $\alpha=0{,}95$, nos dois
  cenários: sob dependência exata o remédio custa um pouco, sob correlação 0,99
  ele **melhora** o erro.
- Na aplicação: número de intercambiáveis com critério solto **contra** critério
  apertado. **Este é o par de números mais importante do trabalho.**

---

## 10. Como verificar cada resultado

O princípio: em vez de olhar a tabela renderizada, confira as **propriedades que
precisam valer**.

| Saída | Invariantes que precisam valer |
| :-- | :-- |
| Validação, caso ortonormal | Discrepância exatamente **zero** — é fórmula fechada |
| Validação, KKT | Violação relativa da ordem de `tol`. Um erro de implementação daria ordem 1 |
| Contraste com 1996 | Nenhuma linha do ridge iterado com objetivo **abaixo** do ótimo de referência; nº de não nulos **cresce** quando o limiar diminui; distância ao ótimo **não** vai a zero quando $\epsilon\to0$ |
| Ablação | As quatro configurações no **mesmo objetivo**; os ajustes que não convergiram reportados **separadamente** |
| Ablação, custo | *Warm start* dominante; ganhos **sub**-multiplicativos |
| EC2 | TPR **não** presa em 1 (senão o desenho é fácil demais); proporção no limite da grade baixa |
| Cópias exatas | Soma estável entre ordens; divisão instável com $\alpha=1$; **igual** com $\alpha<1$ |
| Contribuição | Jaccard = 1 para correlação < 1; < 1 só com dependência exata; amplitude do objetivo ≈ 0 **em todas as células** — se alguma célula tiver amplitude apreciável, o diagnóstico ali é convergência incompleta e a leitura muda |
| Contribuição, remédio | Jaccard volta a 1 em $\alpha<1$; a soma dos coeficientes do grupo é estável entre ordens mesmo quando a divisão não é |
| Aplicação | Intercambiáveis **caem** ao apertar a tolerância; amplitude do objetivo vai a zero |
| Qualquer experimento | O tempo por experimento impresso no fim do `99_roda_tudo.R` diz onde o custo está |

Há um script pronto, `verifica_rotas.R`, que faz isso para o experimento das
rotas sem depender do cache nem do HTML. O mesmo padrão vale para os outros.

---

## 11. Perguntas prováveis da banca

**"Por que não implementou o artigo todo?"**
Porque amplitude produz um trabalho raso. Escolhemos o componente que o artigo de
fato reivindica, que era o que faltava em 1996, e do qual todas as outras partes
dependem. E cortamos um bloco conectado — a Seção 2.4 é a base de 3 e 4 —, não
pedaços avulsos.

**"Como sabe que sua implementação está correta?"**
Por duas verificações contra **teoria**, não contra outro software: o caso
ortonormal reproduz a fórmula fechada de 1996 à precisão de máquina, e as
condições KKT são auditadas ao longo de todo o caminho. Comparar com o `glmnet`
não diria qual das duas implementações está certa.

**"Por que não comparou com o `glmnet`?"**
Justamente por isso. E, na parte de custo, comparar implementações confundiria o
efeito do algoritmo com diferenças de linguagem e compilador. A ablação isola o
efeito dentro da mesma implementação.

**"Seus tempos são muito piores que os do artigo."**
São. Os laços são interpretados em R, os do artigo são Fortran. Por isso o custo
é reportado em **visitas de coordenada**, que é invariante à máquina, e todas as
comparações são relativas.

**"A não unicidade do lasso já é conhecida."**
Sim. O que é nosso é mostrar que a **ordem de varredura** é o mecanismo concreto,
que o fenômeno exige dependência **exata** e não correlação alta, e transformar
isso num diagnóstico que separa degenerescência de erro numérico — coisa que o
*bootstrap* não faz.

**"Como você sabe que a instabilidade da simulação não era só tolerância
frouxa?"**
Porque a **amplitude do objetivo** era nula. Conjuntos diferentes com o mesmo
valor da função objetivo só acontece se os pontos forem igualmente ótimos. Se
fosse tolerância, os reajustes teriam parado em valores diferentes do objetivo —
foi exatamente isso que encontramos nos dados reais, e lá a conclusão foi a
oposta. O mesmo instrumento, duas respostas.

**"Por que a grade pula de 0,99 direto para 1?"**
Porque o desenho foi feito para contrastar dependência exata com correlação alta,
não para localizar a transição. A consequência é que sabemos que existe uma
descontinuidade entre 0,99 e 1, e não onde ela está. Acrescentar $\rho_g=0{,}999$
responderia, e é trabalho de uma execução.

**"Jaccard 1 quer dizer que a seleção está certa?"**
Não, e os próprios dados mostram isso. Em $\rho_g=0{,}99$ o Jaccard é 1 e a
recuperação do grupo verdadeiro é a **pior** das três condições. Estabilidade
significa que a lista não depende da ordem de visita; a lista pode ser
reprodutivelmente errada. São perguntas diferentes e o diagnóstico responde só a
primeira.

**"Por que a auditoria KKT não é circular?"**
Porque o critério de parada é sobre **mudança de coeficiente**, e a auditoria é
sobre **gradiente**. São quantidades diferentes. Se usássemos KKT como critério de
parada, aí sim seria circular — e isso está registrado nas limitações como o
próximo passo, com a ressalva.

**"Alguns ajustes não convergiram."**
Correto, e está reportado. São as configurações sem *warm start*, no fim do
caminho. Consequência: o custo dessas configurações é uma **subestimativa** —
elas pararam antes de resolver. Isso reforça a conclusão em vez de enfraquecê-la.

**"Por que o α = 0,5 na aplicação seleciona mais variáveis?"**
Porque a parte $\ell_2$ não zera nada por si; ela distribui o coeficiente entre
correlacionadas em vez de eleger uma. Modelo maior, seleção mais estável — é o
*trade-off* que a *elastic net* oferece.

---

## 12. Fronteiras e limitações

Diga você antes que perguntem:

- **Escalas absolutas de tempo não são comparáveis às do artigo** (R interpretado
  contra Fortran). Por isso o custo é reportado em visitas de coordenada.
- **Não implementamos** GLMs, *covariance updates*, *sparse updates*,
  *weighted updates*, multinomial nem resposta agrupada.
- **O experimento de velocidade da ordem de varredura usa preditores
  equicorrelacionados**, delineamento favorável ao efeito medido. Correlação
  AR(1) com muitas variáveis ativas é o teste apropriado, e fica pendente.
- **Alguns ajustes sem *warm start* atingem o limite de iterações** no fim do
  caminho. Registrado no objeto e reportado.
- **O diagnóstico detecta a não unicidade, mas não a resolve.** O remédio testado
  ($\alpha<1$) muda o estimador.
- **A grade do experimento de seleção salta de $\rho_g=0{,}99$ para 1 e de
  $\alpha=0{,}95$ para 1.** Há uma descontinuidade em cada intervalo e não
  sabemos onde. Custo de fechar: uma execução.
- **No mesmo experimento, a semente depende de $\rho_g$**, então cada valor de
  correlação usa dados distintos. Comparar linhas de $\rho_g$ mistura
  delineamento com ruído amostral; comparar $\alpha$ dentro de uma linha, não.
- **Estabilidade não é correção.** O diagnóstico responde se a lista depende da
  ordem de visita, não se a lista está certa — e há uma célula em que ela é
  estável e errada.
- **A constante `tol` continua arbitrária.** Um critério que devolvesse um
  certificado explícito — parar quando a violação KKT relativa cai abaixo de um
  limiar — seria mais honesto, ao custo de tornar a auditoria parte do algoritmo
  em vez de verificação independente.
- **Nada aqui trata de inferência após seleção**, que está fora do alcance dos
  dois artigos.

---

## 13. Glossário

Termos mantidos em inglês no relatório, por serem os nomes originais:

| Termo | O que é |
| :-- | :-- |
| *coordinate descent* | Otimizar uma coordenada por vez, com as outras fixas |
| *regularization path* | A sequência de soluções para uma grade decrescente de $\lambda$ |
| *warm start* | Iniciar cada problema na solução do $\lambda$ anterior |
| *active set* | $A=\{j:\beta_j\ne0\}$; ciclos restritos a ele são baratos |
| *soft-thresholding* | $S(z,\gamma)=\operatorname{sign}(z)(\lvert z\rvert-\gamma)^+$ |
| *naive updates* | Manter o resíduo e atualizá-lo em $O(N)$ (Seção 2.1) |
| *covariance updates* | Manter $X^\top r$ e atualizá-lo em $O(p)$ (Seção 2.2) — **não implementado** |
| *cross-validation* | Escolha de $\lambda$ por reamostragem (Seção 6) |
| Visitas de coordenada | Quantas vezes uma coordenada foi avaliada; cada uma custa $O(N)$ |
| Jaccard | $\lvert A\cap B\rvert / \lvert A\cup B\rvert$; = 1 quando o mesmo conjunto é sempre escolhido |
| Regra de 1 EP | Maior $\lambda$ cujo erro de CV fica a menos de um erro-padrão do mínimo |
| *fold* / dobra | Cada uma das $K$ partições da amostra na validação cruzada. Os dois termos são equivalentes; use *fold* por consistência com os demais nomes originais |

---

## 14. Roteiro de apresentação

Os oito itens exigidos, com peso sugerido. Total de referência: 20 minutos.

| Tempo | Item | O que mostrar | A frase que fica |
| --: | :-- | :-- | :-- |
| 1 min | **a) Título** | Artigo, recorte | "Um componente, levado até o fim" |
| 3 min | **b) Introdução** | 1996 e seus três gargalos; o que mudou desde 2010 | "1996 entregou um estimador; 2010 entregou um algoritmo" |
| 1 min | **b′) O recorte** | As três razões | "É a peça que os autores de fato reivindicam" |
| 4 min | **c) Metodologia** | Os três painéis 3D e o corte 1D; o pseudocódigo | "É a quina, não a inclinação, que produz zeros" |
| 1 min | **c′) Validação** | Ortonormal + KKT | "Validamos contra teoria, não contra outro software" |
| 4 min | **d) Monte Carlo** | Ablação (a interação!) e os cinco critérios | "O *active set* sozinho pode ser prejuízo" |
| 2 min | **e) Aplicação** | Colinearidade, caminho, CV | "O ganho aqui é interpretativo, não preditivo" |
| 3 min | **f) Contribuição** | Cópias exatas, perfil do vale, a tabela dos dois sinais, os dois critérios de parada | "Reordenar expõe o sintoma; o objetivo diz a doença" |
| 0,5 min | **g) Uso de IA** | A declaração, com os erros reais | "A IA não executou nada; os números são meus" |
| 1 min | **h) Impacto** | A recomendação prática | "Publique a frequência de seleção sob reordenação" |

**Se o tempo apertar, corte nesta ordem:** o contraste com 1996 (é bonito mas
dispensável), depois o experimento de velocidade da ordem, depois metade do
Monte Carlo II.

**Nunca corte:** os três painéis 3D, a tabela das cópias exatas, e o par de
números do diagnóstico com as duas tolerâncias. São eles que sustentam,
respectivamente, o entendimento do método, a contribuição, e a honestidade do
trabalho.

**E nunca mostre o Jaccard sem a amplitude do objetivo ao lado.** Sozinho, o
Jaccard é ambíguo entre degenerescência e tolerância frouxa — é o par que
conclui, e é o par que a banca vai cobrar. Se o slide couber só um número, que
seja: *listas diferentes, mesmo objetivo.*
