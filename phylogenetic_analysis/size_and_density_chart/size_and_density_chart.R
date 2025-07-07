install.packages(c("readr","ggplot2", "dplyr", "gridExtra", "stringr", "RColorBrewer","ggpubr","scales"))

library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(stringr)
library(RColorBrewer)
library(ggpubr)
library(scales)

# Histogramas ##################################################################

data <- read.csv("E:/ApresentaçãoTCC/Phylogenetic_analysis/size and density chart.csv", sep=";", stringsAsFactors = FALSE) # stringsAsFactors = FALSE evita que as strings sejam convertidas em fatores

colnames(data) # Verificar os nomes das colunas

# Converter vírgulas em pontos para números decimais
data$Genome.Size..mb. <- as.numeric(gsub(",", ".", data$Genome.Size..mb.))
data$X..Quantity.of.ISs <- as.numeric(gsub(",", ".", data$X..Quantity.of.ISs))

# Remover linhas com dados de resumo -> as que contêm "all"
data <- data[!grepl("all", data$Host, ignore.case = TRUE), ]

# Definir a ordem dos hospedeiros
host_order <- c("Fish", "Human", "Bovine")
data$Host <- factor(data$Host, levels = host_order)

# Ordenar primeiro por hospedeiro, depois por tamanho do genoma (crescente)
data <- data %>%
  arrange(Host, Genome.Size..mb.)

# Criar índice de ordenação para manter a ordem no gráfico
data$order_index <- 1:nrow(data)

# Definir cores para cada hospedeiro
host_colors <- c("Fish" = "#3498db", "Human" = "#e74c3c", "Bovine" = "#2ecc71")

# Criar o primeiro histograma (Tamanho do Genoma) ##############################
p1 <- ggplot(data, aes(x = reorder(Strain, order_index), y = Genome.Size..mb., fill = Host)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = Serotype), 
            angle = 90, hjust = -0.1, vjust = 0.5, size = 2.5, color = "black") +
  scale_fill_manual(values = host_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 20, r = 20, b = 60, l = 20, unit = "pt")
  ) +
  labs(
    y = "Genome Size (Mb)",
    title = "Genome Size Distribution by Host",
    fill = "Host"
  ) +
  coord_cartesian(ylim = c(0, max(data$Genome.Size..mb.) * 1.15))

print(p1)

# Criar o segundo histograma (% de TEs) ########################################

p2 <- ggplot(data, aes(x = reorder(Strain, order_index), y = X..Quantity.of.ISs, fill = Host)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = TEs), 
            angle = 90, hjust = 1.1, vjust = 0.5, size = 2, color = "black") +
  scale_fill_manual(values = host_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Strains (ordered by genome size)",
    y = "% Quantity of ISs",
    title = "Transposable Elements Distribution",
    fill = "Host"
  ) +
  coord_cartesian(ylim = c(0, max(data$X..Quantity.of.ISs) * 1.15))

print(p2)

# Combinar os dois gráficos
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))

# Exibir o gráfico combinado
print(combined_plot)

png("E:/histogramas.png", width = 16, height = 20, units = "in", res = 360)
grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))
dev.off()

# Criar o arranjo usando arrangeGrob (em vez de grid.arrange)
combined_plot <- arrangeGrob(p1, p2, ncol = 1, heights = c(1, 1))
# Salvar usando ggsave
ggsave("E:/combined_plot.png", combined_plot, width = 16, height = 20, dpi = 360)

################################################################################

# Ler os dados do CSV
csv <- read.csv("E:/ApresentaçãoTCC/Phylogenetic_analysis/size and density chart.csv", sep=";", stringsAsFactors = FALSE)

colnames(csv) # Verificar os nomes das colunas

# Preparar os dados
dados_clean <- csv %>%
  # Remover linhas de resumo (que começam com "all")
  filter(!grepl("^all", Host)) %>%
  # Renomear colunas para facilitar o uso
  rename(
    genome_size_bp = `Genome.Size..bp.`,
    percent_is = `X..Quantity.of.ISs`
  ) %>%
  # Converter vírgulas para pontos nos dados numéricos
  mutate(
    percent_is = as.numeric(gsub(",", ".", percent_is))
  ) %>%
  # Remover linhas com dados faltantes
  filter(!is.na(genome_size_bp) & !is.na(percent_is) & percent_is > 0)

# Aplicar transformação logarítmica
dados_transformados <- dados_clean %>%
  mutate(
    log10_genome_size = log10(genome_size_bp),
    log10_percent_is = log10(percent_is)
  )

# Função para calcular estatísticas
calcular_stats <- function(data) {
  if(nrow(data) < 3) return(NULL)
  
  # Teste de correlação
  cor_test <- cor.test(data$log10_genome_size, data$log10_percent_is, 
                       method = "pearson")
  
  # Modelo linear para R²
  modelo <- lm(log10_percent_is ~ log10_genome_size, data = data)
  r_squared <- summary(modelo)$r.squared
  
  # Formatação dos valores para exibição (garantindo que não sejam NULL)
  r_formatted <- round(as.numeric(cor_test$estimate), 3)
  p_formatted <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                        round(cor_test$p.value, 3))
  r2_formatted <- round(r_squared, 3)
  
  # Debug: imprimir valores para verificação
  cat("Hospedeiro com", nrow(data), "amostras:\n")
  cat("r =", r_formatted, "p =", p_formatted, "R² =", r2_formatted, "\n")
  
  return(list(
    r = cor_test$estimate,
    p_value = cor_test$p.value,
    r_squared = r_squared,
    label = paste0("r = ", r_formatted, "\n",
                   "p = ", p_formatted, "\n",
                   "R² = ", r2_formatted),
    n = nrow(data)
  ))
}

# Função para criar gráfico individual por hospedeiro
criar_grafico_hospedeiro <- function(host_name, data) {
  dados_host <- data %>% filter(Host == host_name)
  
  if(nrow(dados_host) < 3) {
    return(NULL)
  }
  
  # Calcular estatísticas
  stats <- calcular_stats(dados_host)
  
  # Definir cores por hospedeiro
  cores_hospedeiro <- c("Fish" = "#3498db", "Human" = "#e74c3c", "Bovine" = "#2ecc71")
  cor_host <- cores_hospedeiro[host_name]
  
  # Determinar limites dos eixos baseados nos dados
  x_min <- min(dados_host$log10_genome_size)
  x_max <- max(dados_host$log10_genome_size)
  y_min <- min(dados_host$log10_percent_is)
  y_max <- max(dados_host$log10_percent_is)
  
  # Calcular posição para o texto das estatísticas (canto superior esquerdo)
  x_text <- x_min + 0.02 * (x_max - x_min)
  y_text <- y_max - 0.05 * (y_max - y_min)
  
  # Criar o gráfico
  p <- ggplot(dados_host, aes(x = log10_genome_size, y = log10_percent_is)) +
    geom_point(color = cor_host, size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = cor_host, 
                fill = cor_host, alpha = 0.2, size = 1) +
    labs(
      title = paste("Host:", host_name),
      x = expression(log[10] ~ "Genome size (bp)"),
      y = expression(log[10] ~ "IS density (%)")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  # Adicionar estatísticas como anotação
  if(!is.null(stats)) {
    p <- p + annotate("text", 
                      x = x_text,
                      y = y_text,
                      label = stats$label,
                      hjust = 0, vjust = 1, size = 3.2, 
                      fontface = "plain", color = "black")
  }
  
  # Definir limites dos eixos
  p <- p + coord_cartesian(xlim = c(x_min - 0.01 * (x_max - x_min), 
                                    x_max + 0.01 * (x_max - x_min)),
                           ylim = c(y_min - 0.15 * (y_max - y_min), 
                                    y_max + 0.15 * (y_max - y_min)))
  
  return(p)
}

# Obter hospedeiros únicos
hospedeiros <- unique(dados_transformados$Host)
hospedeiros <- hospedeiros[hospedeiros %in% c("Fish", "Human", "Bovine")]

# Criar gráficos individuais
graficos_lista <- list()
for(host in hospedeiros) {
  grafico <- criar_grafico_hospedeiro(host, dados_transformados)
  if(!is.null(grafico)) {
    graficos_lista[[host]] <- grafico
  }
}

# Organizar os gráficos
if(length(graficos_lista) >= 3) {
  # Gráfico A: Fish
  grafico_A <- graficos_lista[["Fish"]]
  
  # Gráfico B: Human  
  grafico_B <- graficos_lista[["Human"]]
  
  # Gráfico C: Bovine
  grafico_C <- graficos_lista[["Bovine"]]
  
  # Combinar os três gráficos
  grafico_combinado <- grid.arrange(grafico_A, grafico_B, grafico_C, 
                                    ncol = 1, 
                                    top = "Genome Size vs IS Density by Host")
  
  # Salvar o gráfico combinado
  ggsave("E:/line_graphs_by_host.png", grafico_combinado, 
         width = 10, height = 20, dpi = 360)
  
  # Exibir estatísticas detalhadas
  cat("=== ESTATÍSTICAS POR HOSPEDEIRO ===\n\n")
  for(host in names(graficos_lista)) {
    dados_host <- dados_transformados %>% filter(Host == host)
    stats <- calcular_stats(dados_host)
    
    cat("Hospedeiro:", host, "\n")
    cat("Número de amostras:", stats$n, "\n")
    cat("Coeficiente de correlação (r):", round(stats$r, 4), "\n")
    cat("Valor p:", ifelse(stats$p_value < 0.001, "< 0.001", round(stats$p_value, 4)), "\n")
    cat("Coeficiente de determinação (R²):", round(stats$r_squared, 4), "\n")
    cat("---\n")
  }
}

