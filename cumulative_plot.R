# Obter todas as categorias possíveis de clade_name
all_clades <- node_annotation %>%
  arrange(desc(root)) %>%
  dplyr:: select(clade_name) %>%
  unique()

# Filtrar e calcular cumulativo
cumulative_data <- node_annotation %>%
  dplyr::select(-Signature) %>%
  unique() %>%
  arrange(desc(root)) %>%
  group_by(-root, clade_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(cumulative_sum = cumsum(count)) %>%
  rename("root" = "-root")

cumulative_data <- all_clades %>%
  left_join(cumulative_data)

for (i in 2:nrow(cumulative_data)) {
  if (is.na(cumulative_data$cumulative_sum[i])) {
    cumulative_data$cumulative_sum[i] <- cumulative_data$cumulative_sum[i - 1]
  }
}

# Selecionar nomes do eixo x
clade_name <- unique(all_clades$clade_name)

# Calcular 90% do total
total_genes <- sum(na.omit(cumulative_data$count))
threshold_90 <- 0.9 * total_genes
root_90 <- cumulative_data$clade_name[which(cumulative_data$cumulative_sum >= threshold_90)][1]
cumulative_90 <- cumulative_data$cumulative_sum[which(cumulative_data$cumulative_sum >= threshold_90)][1]

# Criar gráfico
ggplot(cumulative_data, aes(x = factor(clade_name, levels = clade_name), y = cumulative_sum)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = threshold_90, linetype = "dashed", color = "red") +
  #geom_text(aes(label=cumulative_sum), vjust = 2, size = 3) +
  annotate("text", x = as.numeric(factor(root_90, levels = rev(all_clades))), y = cumulative_90, 
  label = paste("90% of genes present at", root_90), hjust = 1, vjust = 7, color = "gray") +
  labs(x = "Clade", y = "Cumulative Genes") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(1, 1, 0, 1))
