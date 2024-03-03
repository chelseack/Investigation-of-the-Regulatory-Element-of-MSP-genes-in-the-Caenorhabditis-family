# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 22:46:17 2024

@author: Chelsea
"""

# plot div time vs GC, div time vs length and GC vs length and div time vs similarity
# statistical tests for all four relationships

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress
import statsmodels.api as sm

def read_data(file_path):
    """Read the TSV file containing the data."""
    return pd.read_csv(file_path)

def scatter_plot(x, y, color, xlabel, ylabel):
    """Create a scatter plot."""
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, color=color, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()

def pearson_correlation(x, y):
    """Compute Pearson correlation coefficient and p-value."""
    corr_coef, p_value = pearsonr(x, y)
    print(f"Pearson Correlation Coefficient: {corr_coef}")
    print(f"P-value: {p_value}")

def linear_regression(x, y):
    """Perform linear regression and print the summary."""
    x_with_const = sm.add_constant(x)
    model = sm.OLS(y, x_with_const)
    results = model.fit()
    print(results.summary())

# Read the new TSV file containing the data
# new_data_file = r"C:\Users\kenne\Downloads\EEB498\div_time_GC_length_diff.tsv"
# df_sim_div = r"C:\Users\kenne\Downloads\EEB498\df_sim_div.tsv"
# GC_length=r"C:\Users\kenne\Downloads\EEB498\GC_length.tsv"
# avg_data=r"C:\Users\kenne\Downloads\EEB498\avg_l_GC_no_motifs.tsv"
sd_data=r"C:\Users\kenne\Downloads\EEB498\standard deviation GC and length.tsv"

# # Read data
# new_data_df = read_data(new_data_file)
# sim_div = read_data(df_sim_div)
# GC_length = read_data(GC_length)
# avg_data=read_data(avg_data)
sd_data=read_data(sd_data)

# # Prepare the data for plotting
# divergence_time = new_data_df['Divergence Time']
# diff_length = new_data_df['length']
# diff_GC = new_data_df['GC']
# length = GC_length['length']
# GC=GC_length["GC content"]
# simdiv_time=sim_div['Divergence Time']
# similarity=sim_div['Similarity']
# avg_GC=avg_data['GC-content']
# avg_length=avg_data["Average length"]
# num_of_motif=avg_data["Number of motifs"]
num_of_motif=sd_data["Number of motifs"]
sd_length=sd_data["S.d. length"]
sd_GC=sd_data["S.d. GC-content"]

# # Create scatter plots
# scatter_plot(divergence_time, diff_length, 'blue', 'Divergence Time', 'diff Length')
# scatter_plot(divergence_time, diff_GC, 'green', 'Divergence Time', 'diff GC')
# scatter_plot(diff_GC, diff_length, 'red', 'diff GC', 'diff Length')
# scatter_plot(GC,length,'orange','GC','length')
# scatter_plot(simdiv_time,similarity,'purple','divergence time','similarity')
# scatter_plot(avg_GC,avg_length,'blue','average GC content','average length')
# scatter_plot(num_of_motif,avg_GC,'red','number of motifs found','average GC-content')
# scatter_plot(num_of_motif,avg_length,'blue','number of motifs found','average length')
scatter_plot(num_of_motif,sd_length,'blue','number of motifs found','Sd of length')
scatter_plot(num_of_motif,sd_GC,'red','number of motifs found','Sd of GC-content')
scatter_plot(sd_GC,sd_length,'green','Sd of GC-length','Sd of length')

# # Perform statistical tests
# # Pearson Correlation
# pearson_correlation(divergence_time, diff_length)
# pearson_correlation(divergence_time, diff_GC)
# pearson_correlation(diff_GC, diff_length)
# pearson_correlation(GC,length)
# pearson_correlation(simdiv_time,similarity)

# # Linear Regression
# linear_regression(diff_length, divergence_time)
# linear_regression(diff_GC, divergence_time)
# linear_regression(diff_length, diff_GC)
# linear_regression(GC,length)
# linear_regression(simdiv_time, similarity)



