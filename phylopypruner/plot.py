"Module for generating plots."

import os
from phylopypruner import report

try:
    import matplotlib as mpl
    if "DISPLAY" not in os.environ:
        mpl.use("agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'figure.max_open_warning': 0})
    MATPLOTLIB = True
except ImportError:
    report.tip("install Matplotlib (https://matplotlib.org/) to generate \
plots")
    MATPLOTLIB = False

OCCUPANCY_PLOT_FILE = "occupancy_matrix.png"
FREQ_PLOT_FILE = "paralogy_freq_plot.png"
PLOTS_DIR = "/plots"
PPI = 300
DEFAULT_AXIS_SIZE = 4
FONT_SIZE = 2.5
MAX_PLOT_SIZE = 500


def set_size(labels):
    """
    Takes a list of x- or y labels as an input. Returns the size of the
    provided list times 0.05, if it is larger than 100. Otherwise, return the
    default axis size.
    """
    if len(labels) > 100:
        return 0.05 * len(labels)
    return DEFAULT_AXIS_SIZE


def flag_outliers(axes, genes_below, otus_below):
    "Highlight genes and OTUs below their defined threshold in red."
    # Highlight gene partitions below the allowed threshold in red.
    for index in range(genes_below)[::-1]:
        axes.get_xticklabels()[-index - 1].set_color("red")

    # Highlight OTUs below the allowed threshold in red.
    for index in range(otus_below)[::-1]:
        axes.get_yticklabels()[-index - 1].set_color("red")

    return axes


def get_outliers(total, current_end, below_threshold):
    """
    Return the number of genes that should be flagged in a subplot, given the
    total number of genes, the last subplot index, and the number of genes
    below the threshold in total.
    """
    otus_below, genes_below = below_threshold

    in_this_round = genes_below - (total - current_end)

    if current_end == total and genes_below > MAX_PLOT_SIZE:
        in_this_round = current_end % MAX_PLOT_SIZE
    elif in_this_round > MAX_PLOT_SIZE:
        in_this_round = MAX_PLOT_SIZE

    return (otus_below, in_this_round) if in_this_round >= 0 \
        else (otus_below, 0)

def chunks(items, size):
    """
    Takes a list and an integer as an input. Yield subsets of the list of the
    same size as the provided integer.
    """
    for index in range(0, len(items), size):
        yield items[index:index + size]


def occupancy_subplot(matrix, xlabels, ylabels, from_and_to=None, dir_out=None,
                      below_threshold=None):
    """
    Takes a matrix as a list, a list of X-labels and a list of Y-labels as an
    input and generates an occupancy subplot.
    """
    if from_and_to:
        mat_sub = [None] * len(matrix)
        begin, end = from_and_to
        filename = "occupancy_subplot_{}-{}".format(begin, end)

        for index, _ in enumerate(matrix):
            mat_sub[index] = matrix[index][begin - 1:end - 1]
    else:
        filename = OCCUPANCY_PLOT_FILE
        mat_sub = matrix

    fig = plt.figure()
    plotting_dir = dir_out + PLOTS_DIR + "/"
    axes = fig.add_subplot(111)
    plot = axes.matshow(mat_sub, cmap="BuGn", interpolation="nearest")
    fig.colorbar(plot)

    axes.set_title("Occupancy Matrix")
    axes.xaxis.set_ticks_position("bottom")
    axes.set_xticks(list(range(len(xlabels))))
    axes.set_yticks(list(range(len(ylabels))))
    axes.set_xticklabels(xlabels, rotation="vertical",
                         fontsize=FONT_SIZE, stretch="expanded")
    axes.set_yticklabels(ylabels, fontsize=FONT_SIZE)

    if below_threshold:
        otus_below, genes_below = below_threshold
        axes = flag_outliers(axes, genes_below, otus_below)

    fig.set_size_inches(set_size(xlabels), set_size(ylabels))
    plt.xlabel("Gene partitions")
    plt.ylabel("OTUs")
    # Pad margins so that markers don't get clipped by the axes.
    plt.margins(0.2)
    # Tweak spacing to prevent clipping of tick-labels.
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(plotting_dir + filename, dpi=PPI)


def occupancy_matrix(matrix, xlabels, ylabels, dir_out, below_threshold):
    """
    Generate a heatmap that displays the number of positions covered for each
    gene and OTU.
    """
    report.progress_bar(
        "plotting (disable with '--no-plot')", replace=False)

    if len(xlabels) > MAX_PLOT_SIZE:
        for index, chunk in enumerate(chunks(xlabels, MAX_PLOT_SIZE)):
            begin = (1 + (index * MAX_PLOT_SIZE))

            if len(chunk) < MAX_PLOT_SIZE:
                end = MAX_PLOT_SIZE + ((index - 1) * MAX_PLOT_SIZE) + \
                        len(chunk)
            else:
                end = MAX_PLOT_SIZE + (index * MAX_PLOT_SIZE)

            below = get_outliers(len(xlabels), end, below_threshold)

            occupancy_subplot(
                matrix, chunk, ylabels,
                (begin, end), dir_out, below)
    else:
        occupancy_subplot(matrix=matrix, xlabels=xlabels, ylabels=ylabels,
                          from_and_to=None, dir_out=dir_out,
                          below_threshold=below_threshold)


def paralogy_frequency(indices, frequencies, otus, threshold, dir_out):
    "Generate a paralogy frequency plot."
    try:
        plt.barh(y=indices, color="c0", width=frequencies, alpha=0.5)
    except TypeError:
        report.error("plotting function is lacking indices, try updating \
Matplotlib or disable plotting with '--no-plot'")
        return

    plotting_dir = dir_out + PLOTS_DIR + "/"
    freq_file = plotting_dir + FREQ_PLOT_FILE

    if not os.path.isdir(plotting_dir):
        os.mkdir(plotting_dir)

    if os.path.isfile(freq_file):
        os.remove(freq_file)

    plt.yticks(list(indices), otus)
    plt.ylabel("OTU")
    plt.xlabel("number of paralogs / number of alignments OTU is in")
    plt.title("Paralogy Frequency")

    if threshold:
        plt.axvline(x=threshold,
                    color="red",
                    label="cutoff = {}".format(round(threshold, 3)),
                    linestyle="--")
        plt.legend(loc='upper right', fontsize=8)

    fig = plt.gcf()
    fig.set_size_inches(20.0, 0.15 * len(otus))
    plt.savefig(freq_file, dpi=PPI)
