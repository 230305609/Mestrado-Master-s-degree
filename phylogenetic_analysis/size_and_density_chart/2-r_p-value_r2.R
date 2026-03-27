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
    percent_is = `% Quantity of Families`
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#* Se o arquivo excel contiver linhas com % Quantity of Elements/percent_is igual a zero ou NA, essas linhas devem ser removidas antes de realizar a análise de correlação, pois o log10(0) é indefinido. *# 
#* Alem disso só faz sentido biologicamente analisar apenas linhagens que possuem elementos IS.*#
dados_filtrados <- dados_transformados %>%
  filter(percent_is > 0)
#* Não esquecer de ajustar o código, mudando a variavel dados_transformados para dados_filtrados nos proximos passos do codigo 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

################################################################################
# Calculando a correlação de Pearson (r e p-value) e o Coeficiente de Determinação (R²) de todas as cepas 
################################################################################

# Realizar o teste de correlação de Pearson
## A fórmula ~ y + x especifica as variáveis a serem correlacionadas.
teste_correlacao <- cor.test(
  ~ log10_percent_is + log10_genome_size,
  data = dados_transformados,
  method = "pearson"
)
### Exibir os resultados completos do teste
print(teste_correlacao)

# Cálculo do Coeficiente de Determinação (R²)
## Ajustar um modelo de regressão linear
### A fórmula y ~ x indica que estamos modelando y como uma função de x.
modelo_linear <- lm(log10_percent_is ~ log10_genome_size, data = dados_transformados)
## Obter o resumo do modelo, que inclui o R-quadrado
resumo_modelo <- summary(modelo_linear)
## Exibir o resumo completo
print(resumo_modelo)
### Somente o R-quadrado
r_quadrado <- resumo_modelo$r.squared
cat("O coeficiente de determinação (R-quadrado) é:", r_quadrado, "\n")

################################################################################
# Calculando a correlação de Pearson (r e p-value) e o Coeficiente de Determinação (R²) por hospedeiro 
################################################################################

# Filtrar os dados por hospedeiro e realizar análises separadas

# Função para realizar e imprimir a análise para um subgrupo de hospedeiro
analisar_hospedeiro <- function(nome_hospedeiro) {
  cat("====================================================\n")
  cat("Análise para o Hospedeiro:", nome_hospedeiro, "\n")
  cat("====================================================\n")
  
  # Filtrar o dataframe
  dados_subgrupo <- dados_transformados %>%
    filter(Host == nome_hospedeiro)
  
  # Verificar se há dados suficientes para a análise
  if (nrow(dados_subgrupo) < 3) {
    cat("Não há dados suficientes para a análise.\n\n")
    return(NULL)
  }
  
  # Realizar o teste de correlação
  teste_cor_subgrupo <- cor.test(
    ~ log10_percent_is + log10_genome_size,
    data = dados_subgrupo
  )
  
  # Ajustar o modelo linear
  modelo_linear_subgrupo <- lm(
    log10_percent_is ~ log10_genome_size,
    data = dados_subgrupo
  )
  
  # Imprimir os resultados
  print(teste_cor_subgrupo)
  print(summary(modelo_linear_subgrupo))
  
  # Retornar os resultados para preencher a Tabela 2
  return(list(
    r = teste_cor_subgrupo$estimate,
    p_value = teste_cor_subgrupo$p.value,
    r_squared = summary(modelo_linear_subgrupo)$r.squared,
    n = nrow(dados_subgrupo)
  ))
}

# Executar a análise para cada hospedeiro presente nos dados
## Obter a lista de hospedeiros únicos
hospedeiros_unicos <- unique(dados_transformados$Host)

# Aplicar a função a cada hospedeiro
## Armazenar os resultados em uma lista
resultados_subgrupos <- lapply(hospedeiros_unicos, analisar_hospedeiro)
## Nomear a lista com os nomes dos hospedeiros
names(resultados_subgrupos) <- hospedeiros_unicos

### Exibir os resultados para cada hospedeiro
resultados_subgrupos$Fish
resultados_subgrupos$Human 
resultados_subgrupos$Bovine
