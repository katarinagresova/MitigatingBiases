import scipy.stats as stats

def get_descriptives(background_counts, motif_counts):
    background = background_counts.describe()
    motif = motif_counts.describe()
    return background, motif

def get_comparatives(background_counts, motif_counts): # Wilcoxon Test
    for column in background_counts.columns:
        yield stats.ranksums(background_counts[column], motif_counts[column])

def output_statistics(background_counts, motif_counts):
    with open("statistics.txt", "w") as f:
         print(get_descriptives(background_counts, motif_counts), file=f)
         print(list(get_comparatives(background_counts, motif_counts)), file=f)