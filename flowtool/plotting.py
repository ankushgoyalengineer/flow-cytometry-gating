import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_scatter(df_sample, x_name, y_name):
    """Plain scatter plot of two channels."""
    plt.scatter(df_sample[x_name], df_sample[y_name], s=1, alpha=0.5)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title(f"Scatter plot: {x_name} vs {y_name}")
    plt.show()


def plot_density(df_sample, x_name, y_name, value_range=None):
    """2D density plot with a log colour scale. value_range, e.g. [[0, 1023], [0, 1023]];
    None uses the full range of the data."""
    plt.hist2d(df_sample[x_name], df_sample[y_name], range=value_range, bins=256, norm=LogNorm())
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title(f"Density plot: {x_name} vs {y_name}")
    plt.show()


def plot_histogram(df_sample, num_col):
    """Histogram of one channel: 256 bins over 0-1023, log-scaled counts."""
    sns.histplot(data=df_sample, x=num_col, bins=256, binrange=(0, 1024), kde=False)
    plt.title(f'Distribution of {num_col}')
    plt.yscale("log")
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.show()