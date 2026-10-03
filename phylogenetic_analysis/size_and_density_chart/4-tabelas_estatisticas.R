################################################################################
# Gera as tabelas numéricas para a correlação entre tamanho do genoma e densidade de TEs: n, r de Pearson, IC 95%, p, p ajustado, R² e a verificação dos pressupostos (normalidade, homocedasticidade, linearidade, pontos influentes), além do rho de Spearman e do IC 95% por bootstrap.
#
# Utilizar tanto para o Excel das FAMÍLIAS como para o Excel dos ELEMENTOS:
#
# Saídas (em PASTA_SAIDA), com o PREFIXO escolhido abaixo:
#   <prefixo>_tabelas.xlsx           todas as tabelas + resultados completos
#   <prefixo>_tabelas.docx           as mesmas tabelas, prontas para copiar no Word
#   <prefixo>_texto_resultados.txt   frases com os números reais (revisar e colar)
#   <prefixo>_diagnosticos.png       resíduos x ajustados e Q-Q, por grupo
################################################################################

# ---- CONFIGURAÇÃO ------------------------------------------------------------
ARQUIVO_EXCEL <- "D:/Mestrado/Phylogenetic_analysis/FAMILY - size and density chart.xlsx"
COL_TAMANHO   <- "Genome Size (bp)"
COL_DENSIDADE <- "% Quantity of Families"   # % Quantity of Families/Elements
NOME_ANALISE  <- "famílias de TEs"           # famílias/elementos
PREFIXO       <- "familias"                 # familias/elementos
PASTA_SAIDA   <- "D:/Resultados/tabelas_estatisticas"

HOSPEDEIROS   <- c("Fish", "Human", "Bovine")  # ordem das linhas nas tabelas
ALFA          <- 0.05
N_BOOT        <- 10000   # reamostragens do bootstrap
SEMENTE       <- 2026    # fixa o bootstrap (resultado reprodutível)

# Análise de sensibilidade (opcional): nomes, como aparecem na coluna Strain, das cepas a retirar (p.ex. derivadas isogênicas de uma mesma cepa-mãe, que não são observações independentes). Com c() vazio, essa análise não é gerada.
CEPAS_EXCLUIR <- c()
# ------------------------------------------------------------------------------

pacotes  <- c("readxl", "dplyr", "lmtest", "writexl")
ausentes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(ausentes) > 0) install.packages(ausentes)
invisible(lapply(pacotes, library, character.only = TRUE))

dir.create(PASTA_SAIDA, showWarnings = FALSE, recursive = TRUE)

# ---- Funções de formatação (vírgula decimal, padrão da dissertação) -----------
fmt    <- function(x, d = 3) ifelse(is.na(x), "-", sub("\\.", ",", formatC(x, format = "f", digits = d)))
fmt_p  <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "< 0,001", fmt(p, 3)))
fmt_pe <- function(p, lab = "p") ifelse(p < 0.001, paste0(lab, " < 0,001"), paste0(lab, " = ", fmt(p, 3)))
fmt_ic <- function(lo, hi, d = 3) ifelse(is.na(lo) | is.na(hi), "-",
                                         paste0("[", fmt(lo, d), "; ", fmt(hi, d), "]"))

# ---- Leitura e limpeza (mesmas regras dos scripts 2 e 3) ----------------------
para_num <- function(v) as.numeric(gsub(",", ".", as.character(v)))

bruto <- read_excel(ARQUIVO_EXCEL, sheet = 1)
faltam_colunas <- setdiff(c("Host", COL_TAMANHO, COL_DENSIDADE), names(bruto))
if (length(faltam_colunas) > 0)
  stop("Colunas não encontradas na planilha: ", paste(faltam_colunas, collapse = ", "),
       "\nColunas disponíveis: ", paste(names(bruto), collapse = " | "))

dados <- bruto %>%
  mutate(Host = trimws(as.character(Host))) %>%
  filter(!is.na(Host), !grepl("^all", Host, ignore.case = TRUE)) %>%   # remove linhas-resumo
  mutate(genome_size_bp = para_num(.data[[COL_TAMANHO]]),
         percent_is     = para_num(.data[[COL_DENSIDADE]]))

fora <- setdiff(unique(dados$Host), HOSPEDEIROS)
if (length(fora) > 0)
  warning("Hospedeiros ignorados (não estão em HOSPEDEIROS): ", paste(fora, collapse = ", "))
dados <- filter(dados, Host %in% HOSPEDEIROS)

# ---- Contagem de cepas incluídas / excluídas (transparência) ------------------
amostra <- dados %>%
  group_by(Host) %>%
  summarise(n_planilha       = n(),
            n_sem_valor      = sum(is.na(genome_size_bp) | is.na(percent_is)),
            n_densidade_zero = sum(!is.na(genome_size_bp) & !is.na(percent_is) & percent_is == 0),
            n_analisadas     = sum(!is.na(genome_size_bp) & !is.na(percent_is) & percent_is > 0),
            .groups = "drop") %>%
  arrange(match(Host, HOSPEDEIROS))
total <- amostra %>% select(-Host) %>% summarise(across(everything(), sum)) %>%
  mutate(Host = "Todas as cepas", .before = 1)
amostra <- bind_rows(amostra, total)

# ---- Dados da análise: densidade > 0 e transformação log10 -------------------
analise <- dados %>%
  filter(!is.na(genome_size_bp), !is.na(percent_is), percent_is > 0) %>%
  mutate(log10_tam = log10(genome_size_bp), log10_dens = log10(percent_is))

# ---- Análise de um grupo ------------------------------------------------------
analisar <- function(d, rotulo) {
  n <- nrow(d)
  if (n < 4) return(tibble(Grupo = rotulo, n = n))   # IC de Fisher exige n > 3
  x <- d$log10_tam
  y <- d$log10_dens
  
  out <- tryCatch({
    ct  <- cor.test(x, y, method = "pearson", conf.level = 1 - ALFA)
    mod <- lm(y ~ x)
    sm  <- summary(mod)
    sp  <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    cb  <- confint(mod, level = 1 - ALFA)["x", ]
    
    set.seed(SEMENTE)
    rb <- replicate(N_BOOT, {
      i <- sample.int(n, replace = TRUE)
      suppressWarnings(cor(x[i], y[i]))
    })
    ic_boot <- quantile(rb, c(ALFA / 2, 1 - ALFA / 2), na.rm = TRUE, names = FALSE)
    
    cook   <- cooks.distance(mod)
    infl   <- which(cook > 4 / n)
    cepas  <- if ("Strain" %in% names(d) && length(infl) > 0)
                paste(d$Strain[infl], collapse = "; ") else ""
    sw     <- function(v) tryCatch(shapiro.test(v)$p.value, error = function(e) NA_real_)
    bp     <- tryCatch(bptest(mod)$p.value, error = function(e) NA_real_)
    rs     <- tryCatch(resettest(mod, power = 2:3, type = "fitted")$p.value,
                       error = function(e) NA_real_)
    
    tibble(Grupo = rotulo, n = n,
           r = unname(ct$estimate), ic_inf = ct$conf.int[1], ic_sup = ct$conf.int[2],
           t = unname(ct$statistic), gl = unname(ct$parameter), p = ct$p.value,
           r2 = sm$r.squared, r2_aj = sm$adj.r.squared,
           beta = unname(coef(mod)["x"]), beta_inf = cb[[1]], beta_sup = cb[[2]],
           rho = unname(sp$estimate), p_rho = sp$p.value,
           boot_inf = ic_boot[1], boot_sup = ic_boot[2],
           sw_res = sw(residuals(mod)), bp = bp, reset = rs,
           cook_max = max(cook), n_infl = length(infl), cepas_infl = cepas,
           sw_tam_bruto = sw(d$genome_size_bp), sw_tam_log = sw(x),
           sw_dens_bruto = sw(d$percent_is),    sw_dens_log = sw(y))
  }, error = function(e) {
    message("Grupo '", rotulo, "': análise não realizada (", conditionMessage(e), ")")
    tibble(Grupo = rotulo, n = n)
  })
  out
}

calcular <- function(df) {
  grupos <- c(setNames(lapply(HOSPEDEIROS, function(h) filter(df, Host == h)), HOSPEDEIROS),
              list("Todas as cepas" = df))
  res <- bind_rows(Map(analisar, grupos, names(grupos)))
  res$p_holm <- p.adjust(res$p, method = "holm")   # correção para as 4 correlações da tabela
  list(grupos = grupos, res = res)
}

# ---- Tabelas formatadas -------------------------------------------------------
montar_tabelas <- function(res, amostra = NULL) {
  tabs <- list(
    Tab1_correlacao = res %>% transmute(
      Grupo, n,
      r = fmt(r), `IC 95% de r` = fmt_ic(ic_inf, ic_sup), `R²` = fmt(r2),
      p = fmt_p(p), `p ajustado (Holm)` = fmt_p(p_holm),
      `ρ de Spearman` = fmt(rho), `p (Spearman)` = fmt_p(p_rho)),
    Tab2_pressupostos = res %>% transmute(
      Grupo, n,
      `Shapiro-Wilk (p)` = fmt_p(sw_res),
      `Breusch-Pagan (p)` = fmt_p(bp),
      `RESET (p)` = fmt_p(reset),
      `D de Cook máx.` = fmt(cook_max),
      `Cepas com D > 4/n` = ifelse(is.na(n_infl), "-", as.character(n_infl)),
      `IC 95% (bootstrap)` = fmt_ic(boot_inf, boot_sup)),
    Tab3_normalidade = res %>% transmute(
      Grupo, n,
      `Tamanho (pb)` = fmt_p(sw_tam_bruto), `Tamanho (log10)` = fmt_p(sw_tam_log),
      `Densidade (%)` = fmt_p(sw_dens_bruto), `Densidade (log10)` = fmt_p(sw_dens_log)))
  if (!is.null(amostra))
    tabs$Tab4_amostra <- amostra %>% transmute(
      Grupo = Host, `Cepas na planilha` = n_planilha, `Sem valor (NA)` = n_sem_valor,
      `Densidade = 0 (excluídas)` = n_densidade_zero, `Cepas analisadas` = n_analisadas)
  tabs
}

principal <- calcular(analise)
res       <- principal$res
tabs      <- montar_tabelas(res, amostra)

sens <- NULL
if (length(CEPAS_EXCLUIR) > 0 && !("Strain" %in% names(analise)))
  warning("CEPAS_EXCLUIR preenchido, mas a planilha não tem a coluna 'Strain': sensibilidade não gerada.")
if (length(CEPAS_EXCLUIR) > 0 && "Strain" %in% names(analise)) {
  sens <- calcular(filter(analise, !(Strain %in% CEPAS_EXCLUIR)))
  tabs$Tab5_sensibilidade <- montar_tabelas(sens$res)$Tab1_correlacao
}

# ---- Texto com os números reais -------------------------------------------------
sentenca <- function(r) {
  if (is.na(r$r)) return(paste0(r$Grupo, " (n = ", r$n, "): amostra insuficiente para o teste."))
  paste0(r$Grupo, " (n = ", r$n, "): r = ", fmt(r$r), " (IC 95% ", fmt_ic(r$ic_inf, r$ic_sup),
         "); t(", r$gl, ") = ", fmt(r$t, 2), "; ", fmt_pe(r$p), "; R² = ", fmt(r$r2), "; ",
         fmt_pe(r$p_holm, "p ajustado (Holm)"), " - ",
         ifelse(r$p_holm < ALFA,
                paste0("correlação ", ifelse(r$r > 0, "positiva", "negativa"), " significativa após a correção"),
                "sem significância estatística após a correção"),
         ". ρ de Spearman = ", fmt(r$rho), " (", fmt_pe(r$p_rho), ").")
}
avisos <- vapply(seq_len(nrow(res)), function(i) {
  r <- res[i, ]
  if (is.na(r$r)) return(paste0(r$Grupo, ": amostra insuficiente."))
  f <- c()
  if (!is.na(r$sw_res) && r$sw_res < ALFA)
    f <- c(f, paste0("resíduos fora da normalidade (Shapiro-Wilk, ", fmt_pe(r$sw_res), ")"))
  if (!is.na(r$bp) && r$bp < ALFA)
    f <- c(f, paste0("variância dos resíduos não constante (Breusch-Pagan, ", fmt_pe(r$bp), ")"))
  if (!is.na(r$reset) && r$reset < ALFA)
    f <- c(f, paste0("indício de não linearidade (RESET, ", fmt_pe(r$reset), ")"))
  if (!is.na(r$rho) && (sign(r$r) != sign(r$rho) || (r$p < ALFA) != (r$p_rho < ALFA)))
    f <- c(f, "Pearson e Spearman divergem (sinal ou significância): investigar pontos influentes")
  if (r$n_infl > 0)
    f <- c(f, paste0(r$n_infl, " ponto(s) com D de Cook > 4/n a conferir",
                     ifelse(nzchar(r$cepas_infl), paste0(": ", r$cepas_infl), "")))
  paste0(r$Grupo, ": ", ifelse(length(f) > 0, paste(f, collapse = "; "),
                               "nenhum desvio detectado nos testes"))
}, character(1))

texto <- c(
  paste0("RESULTADOS - ", NOME_ANALISE),
  paste0("Correlação entre log10(tamanho do genoma, pb) e log10(densidade de ", NOME_ANALISE, ", %)."),
  paste0("Cepas com densidade = 0 foram excluídas (ver Tab4_amostra). IC 95% de r: transformação z de Fisher."),
  paste0("Bootstrap percentil: ", N_BOOT, " reamostragens, semente ", SEMENTE, "."),
  "",
  vapply(seq_len(nrow(res)), function(i) sentenca(res[i, ]), character(1)),
  "",
  "PRESSUPOSTOS (o que discutir no texto):",
  avisos)
if (!is.null(sens)) {
  texto <- c(texto, "", paste0("SENSIBILIDADE (sem ", length(CEPAS_EXCLUIR), " cepa(s) excluída(s)):"),
             vapply(seq_len(nrow(sens$res)), function(i) sentenca(sens$res[i, ]), character(1)))
}
con <- file(file.path(PASTA_SAIDA, paste0(PREFIXO, "_texto_resultados.txt")),
            open = "w", encoding = "UTF-8")
writeLines(texto, con); close(con)

# ---- Excel ----------------------------------------------------------------------
write_xlsx(c(tabs, list(Resultados_completos = as.data.frame(res))),
           file.path(PASTA_SAIDA, paste0(PREFIXO, "_tabelas.xlsx")))

# ---- Figura de diagnóstico ----------------------------------------------------
grupos <- principal$grupos
png(file.path(PASTA_SAIDA, paste0(PREFIXO, "_diagnosticos.png")),
    width = 3200, height = 1600, res = 300)
par(mfcol = c(2, length(grupos)), mar = c(4, 4, 2.6, 1), cex.main = 0.9, cex.lab = 0.85)
for (g in names(grupos)) {
  d <- grupos[[g]]
  if (nrow(d) < 4) { plot.new(); plot.new(); next }
  m <- lm(log10_dens ~ log10_tam, data = d)
  plot(fitted(m), resid(m), pch = 19, col = "#00000088", xlab = "Valores ajustados",
       ylab = "Resíduos", main = paste0(g, " (n = ", nrow(d), ")"))
  abline(h = 0, lty = 2)
  qqnorm(resid(m), pch = 19, col = "#00000088", main = "",
         xlab = "Quantis teóricos", ylab = "Quantis dos resíduos")
  qqline(resid(m))
}
invisible(dev.off())

# ---- Word (três linhas, sem bordas laterais: padrão ABNT/IBGE) -----------------
if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  library(flextable); library(officer)
  ft_abnt <- function(df) {
    flextable(df) %>% theme_booktabs() %>%
      font(fontname = "Arial", part = "all") %>% fontsize(size = 9, part = "all") %>%
      align(j = seq_len(ncol(df))[-1], align = "center", part = "all") %>%
      autofit() %>% fit_to_width(6.2)
  }
  doc <- read_docx()
  add <- function(doc, titulo, df, nota) {
    doc <- body_add_par(doc, titulo, style = "Normal")
    doc <- body_add_flextable(doc, ft_abnt(df))
    doc <- body_add_par(doc, nota, style = "Normal")
    body_add_par(doc, "", style = "Normal")
  }
  fonte <- " Fonte: elaborado pelo autor."
  doc <- add(doc,
    paste0("Tabela X – Correlação de Pearson entre o log10 do tamanho do genoma (pb) e o log10 da densidade de ",
           NOME_ANALISE, " (%), por hospedeiro e no conjunto das cepas"),
    tabs$Tab1_correlacao,
    paste0("Nota: n = cepas analisadas (densidade > 0); IC 95% de r pela transformação z de Fisher; ",
           "p ajustado pelo método de Holm para as quatro correlações da tabela; ",
           "ρ = correlação de postos de Spearman.", fonte))
  doc <- add(doc,
    "Tabela Y – Verificação dos pressupostos do modelo linear (log10 da densidade em função do log10 do tamanho do genoma)",
    tabs$Tab2_pressupostos,
    paste0("Nota: Shapiro-Wilk aplicado aos resíduos; Breusch-Pagan (homocedasticidade); RESET (forma funcional); ",
           "D de Cook > 4/n usado apenas para triagem de pontos influentes; ",
           "IC 95% de r por bootstrap percentil (", N_BOOT, " reamostragens).", fonte))
  doc <- add(doc,
    "Tabela Z – Teste de normalidade de Shapiro-Wilk (valor de p) das variáveis antes e depois da transformação log10",
    tabs$Tab3_normalidade, paste0("Nota: valores de p < 0,05 indicam desvio da normalidade.", fonte))
  doc <- add(doc,
    "Tabela W – Cepas incluídas e excluídas da análise de correlação, por hospedeiro",
    tabs$Tab4_amostra,
    paste0("Nota: cepas sem elementos (densidade = 0) foram excluídas porque log10(0) é indefinido.", fonte))
  if (!is.null(sens))
    doc <- add(doc,
      paste0("Tabela V – Análise de sensibilidade: correlação de Pearson sem ", length(CEPAS_EXCLUIR),
             " cepa(s) potencialmente não independente(s)"),
      tabs$Tab5_sensibilidade, paste0("Nota: mesmas convenções da Tabela X.", fonte))
  print(doc, target = file.path(PASTA_SAIDA, paste0(PREFIXO, "_tabelas.docx")))
}

# ---- Resumo no console ------------------------------------------------------------
cat("\n", paste(texto, collapse = "\n"), "\n\nArquivos salvos em: ", PASTA_SAIDA, "\n", sep = "")
