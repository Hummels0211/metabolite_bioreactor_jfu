library(readr)
library(ggplot2)
library(igraph)
library(ggraph)
library(RColorBrewer)

# Load the correlation result table
data <- read_csv('/home/jiayu/Documents/Pepe2Jia_10_02_2021/network_spearman_correlation_data_12032022.csv')

# Load the 16s-seq mapping result table
sc_mapping <- read_csv('/home/jiayu/Documents/SynCom_Info/Syncom_mapping_2019.csv')

# Reduce the complexity of the mapping table
sc_mapping <- sc_mapping[, c('ID', 'family', 'genus')]

# Adding the taxonomic information into the correlation table
data$Sample1_family <- sc_mapping$family[match(data$Sample1, sc_mapping$ID)]
data$Sample2_family <- sc_mapping$family[match(data$Sample2, sc_mapping$ID)]

data$Sample1_genus <- sc_mapping$genus[match(data$Sample1, sc_mapping$ID)]
data$Sample2_genus <- sc_mapping$genus[match(data$Sample2, sc_mapping$ID)]

# Create a hierarchy vector
hierarchy <- sc_mapping[, c('family', 'ID')]
# Add the original point for the brach structure
df_ext <- data.frame(family = 'Origin', ID = c(unique(sc_mapping$family)))
hierarchy <- rbind(hierarchy, df_ext)

# Create a vertical indice list 
vertices <- data.frame(
	name = unique(c(as.character(hierarchy$family), as.character(hierarchy$ID)))
)
vertices$Family <- hierarchy$family[match(vertices$name, hierarchy$ID)]
vertices$id <- NA
all_leaves <- which(is.na(match(vertices$name, hierarchy$from)))
n_leaves <- length(all_leaves)
vertices$id[all_leaves] <- seq(1:n_leaves)
# Modify the angles for the labelling text
# vertices$angle <- 360 * vertices$id / n_leaves

# Create a list for edge colouring
Correlation <- rep(NA, 339)  # Not clear why the length is 339
for(i in 1:length(data$r)){	# Replace the 'key' position with correlation results
	Correlation[5*i - 6] <- (data$r)[i]
}
# Two odd positions need to be manually changed
Correlation[1] <- data$r[1]
Correlation[6] <- data$r[2]

# Create a canvas object
canvas <- graph_from_data_frame(hierarchy, vertices = vertices)

# Create a connection object
from <- match(data$Sample1, vertices$name)
to <- match(data$Sample2, vertices$name)

# Create the figure
png(file = 'Classified_correlation_chord_diagram.png', width = 600, height = 500)

fig <- ggraph(canvas, layout = 'dendrogram', circular = TRUE) +
	labs(title = 'The complete Spearman correlation co-occurrence network') +
	geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.65, tension = 0.45, aes(colour = Correlation), width = 2.5) +
	scale_edge_colour_distiller(palette = 'RdYlBu', guide = 'edge_colourbar') +
	scale_fill_continuous(limits = c(-1.0, 1.0), breaks = c(-1.0, -0.5, 0.0, 0.5, 1.0), guide = guide_colourbar(reverse = TRUE)) +
	geom_node_point(aes(filter = leaf, x = x * 1.12, y = y * 1.12, colour = Family), size = 7.5, alpha = 1.0) +
	geom_node_text(aes(x = x * 1.3, y = y * 1.25, filter = leaf, label = name, colour = Family, angle = 0, hjust = 0.5, vjust = 0.5), size = 5.5, alpha = 1.0) +
	scale_colour_manual(values = rep(brewer.pal(12, 'Paired'), 30)) +
	theme_void() +
       	theme(
		legend.key.size = unit(1.1, 'cm'),
		legend.title = element_text(size = 20, family='AvantGarde', face = 'bold'),
		legend.text = element_text(size=12, face = 'italic', family = 'Helvetica'),
		plot.margin = unit(c(0.3,0.3,0.3,0.3),"cm"),
	        plot.title = element_text(size = 25, family = 'AvantGarde', face = 'bold')
	) +
	expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))


# Save the plot as a .png file
dev.off()
