import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns
from os import path

from plotting.plotting import melt_dataframe, create_box_plots, create_distribution_plot, create_pairplot


ALL_VARS_PLOTTING_FUNCTIONS = [create_box_plots, create_pairplot]
SUM_PLOTTING_FUNCTION = [create_distribution_plot]


def plot_functions(background, motif, output_directory, format='pdf', size=(20, 16),
                   all_vars_plotting_functions=None,
                   sum_plotting_functions=None) -> bool:
    """
    A function to automate creating multiple functions
    :param sum_plotting_functions: the plotting functions to use that relate to the sum of the genetic sequence
    :param all_vars_plotting_functions: the functions that plot each column individually
    :param background: the background dataset
    :param motif: the motif dataset
    :param output_directory: the output directory
    :param format: the format
    :return: True IFF the function executed successfully
    """
    if sum_plotting_functions is None:
        sum_plotting_functions = SUM_PLOTTING_FUNCTION
    if all_vars_plotting_functions is None:
        all_vars_plotting_functions = ALL_VARS_PLOTTING_FUNCTIONS
    background['source'] = 'background'
    motif['source'] = 'motif'
    dataset = pd.concat([background, motif])
    melted_dataset = melt_dataframe(dataset)

    for function in all_vars_plotting_functions:
        fig, ax = plt.subplots(figsize=size)
        _ = function(melted_dataset, ax=ax)
        fig.savefig(path.join(output_directory, f"{function.__name__.replace('create_', '')}.{format}"), format=format,
                    dpi=1200)

    for function in sum_plotting_functions:
        fig, ax = plt.subplots(figsize=size)
        _ = function(background, motif, ax=ax)
        fig.savefig(path.join(output_directory, f"{function.__name__.replace('create_', '')}.{format}"), format=format,
                    dpi=1200)

    return True