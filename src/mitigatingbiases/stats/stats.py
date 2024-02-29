from scipy.stats import fisher_exact, ranksums
import numpy as np
from statsmodels.stats.multitest import fdrcorrection

def stats_composition_comparison(pos_df, neg_df, end_position, p_value_thresh=0.01):

    # get columns names
    bases = pos_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    p_values = {}
    for base in bases:

        pos_freq = pos_df[base][:end_position] / pos_df['sum'][:end_position]
        neg_freq = neg_df[base][:end_position] / neg_df['sum'][:end_position]

        _, p_value = ranksums(pos_freq, neg_freq)
        p_values[base] = p_value

    # Correcting for FDR
    _, new_p_values = fdrcorrection(list(p_values.values()))
    p_values = dict(zip(p_values.keys(), new_p_values))

    # If there is no significant p-value, this test passed
    passed = np.all(np.array(new_p_values) > p_value_thresh)
    
    return (p_values, passed)

def stats_per_base_sequence_comparison(pos_df, neg_df, end_position, p_value_thresh=0.01):

    # get columns names
    bases = pos_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    stats = {}
    passed = True
    for base in bases:

        pos_freq = pos_df[base][:end_position] / pos_df['sum'][:end_position]
        neg_freq = neg_df[base][:end_position] / neg_df['sum'][:end_position]

        stats[base] = []
        for i in range(len(pos_freq)):
            table=[[pos_freq[i] * 100, (1 - pos_freq[i]) * 100],
                [neg_freq[i] * 100, (1 - neg_freq[i]) * 100]]

            _, p_value = fisher_exact(table=table) 
            stats[base].append(p_value)

        # Correcting for FDR
        _, stats[base] = fdrcorrection(stats[base])
        # If there is no significant p-value, this test passed
        passed = passed and np.all(np.array(stats[base]) > p_value_thresh)

    return (stats, passed)

def stats_lenght_comparison(pos_df, neg_df, p_value_thresh=0.01):

    _, p_value = ranksums(pos_df, neg_df)
    # If there is no significant p-value, this test passed
    passed = p_value > p_value_thresh

    return (p_value, passed)