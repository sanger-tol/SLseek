librerary(ggplot2)
library(scales)
library(optparse)

# getting data from terminal
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input CSV file with k-mer histogram data", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="kmer_histogram.png", 
              help="Output plot file name [default= %default]", metavar="character"),
  make_option(c("-a", "--abundance_col"), type="character", default="abundance", 
              help="Column name for k-mer abundance [default= %default]", metavar="character"),
  make_option(c("-c", "--count_col"), type="character", default="X.of.kmers", 
              help="Column name for number of k-mers [default= %default]", metavar="character"),
  make_option(c("-x", "--xlim"), type="character", default=NULL, 
              help="X-axis limits as 'min,max' (optional)", metavar="character"),
  make_option(c("-y", "--ylim"), type="character", default=NULL, 
              help="Y-axis limits as 'min,max' (optional)", metavar="character"),
  make_option(c("-t", "--title"), type="character", default="K-mer Abundance Histogram", 
              help="Plot title [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be supplied", call.=FALSE)
}
# Read the data
data <- read.csv(opt$input)
# Parse xlim and ylim if provided
xlim <- if (!is.null(opt$xlim)) as.numeric(unlist(strsplit(opt$xlim, ","))) else NULL
ylim <- if (!is.null(opt$ylim)) as.numeric(unlist(strsplit(opt$ylim, ","))) else NULL


plot_kmer_histogram_gg <- function(data,
                                   abundance_col = "abundance",
                                   count_col = "X.of.kmers",
                                   xlim = NULL,
                                   ylim = NULL,
                                   title = "K-mer Abundance Histogram") {
  # Base plot
  p <- ggplot(data, aes_string(x = abundance_col, y = count_col)) +
    geom_col(fill = "steelblue") +
    labs(x = "Abundance",
         y = "Number of k-mers",
         title = title) +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(labels = comma) +
    theme_minimal()
  
  # Add x/y limits only if provided
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  # print(p)
  # Return the plot object
    return(p)
}

# Generate and save the plot
kmer_plot <- plot_kmer_histogram_gg(data,
                                    abundance_col = opt$abundance_col,
                                    count_col = opt$count_col,
                                    xlim = xlim,
                                    ylim = ylim,
                                    title = opt$title)
ggsave(opt$output, plot = kmer_plot, width = 10, height = 6)