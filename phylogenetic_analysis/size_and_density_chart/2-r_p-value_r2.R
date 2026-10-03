# Utilizar este script .R tanto para o arquivo excel das Familias como o arquivo dos elementos de TEs, so modificar o path do arquivo excel e fazer algumas pequenas alteraações durante o código.

install.packages(c("readxl", "dplyr", "ggplot2", "ggpubr"))

library(readxl) # Para ler arquivos Excel
library(dplyr) # Para manipulação de dados
library(ggplot2) # Para visualização de dados
library(ggpubr) # Para funções adicionais de visualização

dados_bacterias <- read_excel("D:/Mestrado/Phylogenetic_analysis/FAMILY - size and density chart.xlsx", sheet = 1)

head(dados_bacterias)
str(dados_bacterias)

# Renomear colunas para facilitar o uso, removendo espaços e caracteres especiais
## Colunas -> "Genome Size (bp)" e "% Quantity ..."
dados_bacterias <- dados_bacterias %>%
  rename(
    genome_size_bp = `Genome Size (bp)`,
    percent_is = `% Quantity of Families` # Families or Elements, depende do excel
  )

## Converter ambas as colunas para numérico ANTES da transformação log
dados_bacterias$genome_size_bp <- as.numeric(dados_bacterias$genome_size_bp)
dados_bacterias$percent_is <- as.numeric(dados_bacterias$percent_is)

### Verificar se a conversão funcionou
str(dados_bacterias$genome_size_bp)
str(dados_bacterias$percent_is)

# Aplicar a transformação logarítmica de base 10
dados_transformados <- dados_bacterias %>%
  mutate(
    log10_genome_size = log10(genome_size_bp),
    log10_percent_is = log10(percent_is)
  )

## Verificar os dados transformados
head(dados_transformados)
str(dados_transformados$genome_size_bp)
str(dados_transformados$percent_is)

################################################################################
#* Se o arquivo excel contiver linhas com % Quantity of Elements/percent_is igual a zero ou NA, essas linhas devem ser removidas antes de realizar a análise de correlação, pois o log10(0) é indefinido. *# 
#* Alem disso só faz sentido biologicamente analisar apenas linhagens que possuem elementos transponiveis (IS e TNs).*#
dados_filtrados <- dados_transformados %>%
  filter(percent_is > 0)
#* Não esquecer de ajustar o código, mudando a variavel dados_transformados para dados_filtrados nos proximos passos do codigo 

################################################################################
# 2. Correlação de Pearson e coeficiente de determinação de todas as cepas
################################################################################

# Interromper a análise com uma mensagem clara se não houver observações
# suficientes para calcular a correlação.
if (nrow(dados_filtrados) < 3) {
  stop("São necessárias pelo menos 3 observações válidas para a correlação de Pearson.")
}

# Calcular a correlação de Pearson entre as variáveis transformadas em log10.
teste_correlacao <- cor.test(
  ~ log10_percent_is + log10_genome_size,
  data = dados_filtrados,
  method = "pearson"
)

# Exibir os resultados completos do teste
print(teste_correlacao)

# Ajustar o modelo de regressão linear: percentual de TEs em função do tamanho
# do genoma, com ambas as variáveis transformadas em log10.
modelo_linear <- lm(
  log10_percent_is ~ log10_genome_size,
  data = dados_filtrados
)

# Obter e exibir o resumo do modelo, incluindo o R²
resumo_modelo <- summary(modelo_linear)
print(resumo_modelo)

# Extrair os principais resultados da análise geral
r_geral <- unname(teste_correlacao$estimate)
p_geral <- teste_correlacao$p.value
r2_geral <- resumo_modelo$r.squared
n_geral <- nrow(dados_filtrados)

cat("\n====================================================\n")
cat("RESULTADOS GERAIS\n")
cat("====================================================\n")
cat("Número de observações (n):", n_geral, "\n")
cat("Correlação de Pearson (r):", round(r_geral, 4), "\n")
cat("p-value:", format.pval(p_geral, digits = 4, eps = 0.0001), "\n")
cat("Coeficiente de determinação (R²):", round(r2_geral, 4), "\n")
cat("====================================================\n")

################################################################################
# 3. Correlação de Pearson e R² por hospedeiro
################################################################################

analisar_hospedeiro <- function(nome_hospedeiro) {
  cat("\n====================================================\n")
  cat("Hospedeiro:", as.character(nome_hospedeiro), "\n")
  cat("====================================================\n")
  
  # Filtrar somente as observações do hospedeiro atual.
  dados_subgrupo <- dados_filtrados %>%
    filter(Host == nome_hospedeiro)
  
  # Remover pares não finitos, caso ainda exista algum valor problemático.
  dados_subgrupo <- dados_subgrupo %>%
    filter(
      is.finite(log10_genome_size),
      is.finite(log10_percent_is)
    )
  
  numero_observacoes <- nrow(dados_subgrupo)
  
  # A correlação de Pearson exige pelo menos três pares completos. Também é
  # necessário que as duas variáveis apresentem variação dentro do subgrupo.
  if (
    numero_observacoes < 3 ||
    dplyr::n_distinct(dados_subgrupo$log10_genome_size) < 2 ||
    dplyr::n_distinct(dados_subgrupo$log10_percent_is) < 2
  ) {
    cat("Análise não realizada: dados insuficientes ou sem variação.\n")
    
    return(data.frame(
      Host = as.character(nome_hospedeiro),
      n = numero_observacoes,
      r = NA_real_,
      p_value = NA_real_,
      R2 = NA_real_,
      status = "Dados insuficientes ou sem variação",
      stringsAsFactors = FALSE
    ))
  }
  
  teste_subgrupo <- cor.test(
    ~ log10_percent_is + log10_genome_size,
    data = dados_subgrupo,
    method = "pearson"
  )
  
  modelo_subgrupo <- lm(
    log10_percent_is ~ log10_genome_size,
    data = dados_subgrupo
  )
  
  r_subgrupo <- unname(teste_subgrupo$estimate)
  p_subgrupo <- teste_subgrupo$p.value
  r2_subgrupo <- summary(modelo_subgrupo)$r.squared
  
  cat("Número de observações (n):", numero_observacoes, "\n")
  cat("Correlação de Pearson (r):", round(r_subgrupo, 4), "\n")
  cat("p-value:", format.pval(p_subgrupo, digits = 4, eps = 0.0001), "\n")
  cat("Coeficiente de determinação (R²):", round(r2_subgrupo, 4), "\n")
  
  data.frame(
    Host = as.character(nome_hospedeiro),
    n = numero_observacoes,
    r = r_subgrupo,
    p_value = p_subgrupo,
    R2 = r2_subgrupo,
    status = "Análise realizada",
    stringsAsFactors = FALSE
  )
}

# Obter a lista de hospedeiros únicos, descartando valores ausentes e vazios.
hospedeiros_unicos <- unique(dados_filtrados$Host)
hospedeiros_unicos <- hospedeiros_unicos[
  !is.na(hospedeiros_unicos) & trimws(as.character(hospedeiros_unicos)) != ""
]

# Aplicar a função a cada hospedeiro e consolidar os resultados em um dataframe.
resultados_por_hospedeiro <- lapply(
  hospedeiros_unicos,
  analisar_hospedeiro
) %>%
  bind_rows()

# Exibir a tabela consolidada
print(resultados_por_hospedeiro)

################################################################################
# 4. Gráfico geral da relação analisada
################################################################################

grafico_geral <- ggplot(
  dados_filtrados,
  aes(x = log10_genome_size, y = log10_percent_is)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Relação entre o tamanho do genoma e o percentual de TEs",
    subtitle = paste0(
      "Pearson r = ", round(r_geral, 3),
      "; p = ", format.pval(p_geral, digits = 3, eps = 0.001),
      "; R² = ", round(r2_geral, 3),
      "; n = ", n_geral
    ),
    x = "log10 do tamanho do genoma (bp)",
    y = "log10 do percentual de TEs"
  ) +
  theme_classic()

print(grafico_geral)
