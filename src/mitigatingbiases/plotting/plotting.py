import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

sns.set_palette("pastel")


def create_distribution_plot(background: pd.DataFrame, motif: pd.DataFrame,
                             variable='sum',
                             title='Background length, vs Motif length', ax=None):
    """
    A function to create the Distribution plots for the background
    :param variable: the variable
    :param title: the title
    :param background: the background dataframe
    :param motif: the motif dataframe
    :return: the axis plot
    """
    ax = sns.kdeplot(background[variable], label='Background', fill=True, ax=ax)
    ax = sns.kdeplot(motif[variable], label='Background', fill=True, ax=ax)
    plt.title(title)
    return ax


def create_box_plots(dataset,
                     variable_name='nucleotide',
                     value_name='frequency',
                     source_variable='source',
                     title='Box plot comparing the frequency of nucleotides, stratified by background & motif',
                     ax=None):
    """
    A function to create the box plots
    :param ax:
    :param dataset: the dataset
    :param variable_name: the name to use on the X axis
    :param value_name: the value name to use n the Y axis
    :param source_variable: the source variable
    :param title: the title to use the plot
    :return:
    """
    ax = sns.boxplot(data=dataset, x=variable_name, y=value_name, hue=source_variable, ax=ax)
    plt.title(title)
    return ax


def melt_dataframe(dataset,
                   variables=None,
                   variable_name='nucleotide',
                   value_name='frequency',
                   source_variable='source',
                   ):
    """
    A function to melt the dataframe in 3 columns (source, nucleotide, frequency)
    :param dataset:
    :param value_name:
    :param variable_name:
    :param variables:
    :param source_variable:
    :return: the melted dataframe
    """
    if variables is None:
        variables = ['A', 'C', 'G', 'T', ]
    return pd.melt(dataset, id_vars=source_variable, value_vars=variables, var_name=variable_name,
                   value_name=value_name)


def create_pairplot(dataset: pd.DataFrame, title='Bar plot comparing the frequency', *_, **__):
    """
    A function to create the bar plot
    :param dataset: the dataset
    :param title: the title
    :return: the axis
    """
    ax = sns.pairplot(dataset, hue='source',)
    plt.title(title)
    return ax