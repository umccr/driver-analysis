"""
Contains functions to adjust OncodriveCLUSTL models
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

def calculate_ks(file, pvalue, random_set=1000, alpha=0.1):
    """
    Calculate Kolmogorov-Smirnov statistic (KS) for raw-pvalues fitness to the uniform distribution (expected).
    KS is returned as positive if the difference between observed and expected is greater than 0 in subset of
    p-values defined by alpha cutoff.

    Args:
        file (str): path to file with elements results
        pvalue (str): p-value column name in dataframe to calculate KS
        random_set (int): random set of genes to calculate KS
        alpha (float): cutoff to calculate KS sign

    Returns:
        ks_statistic (float): KS statistic
        (float): ks_pvalue
        (int): number of p-values used to calculate KS
    """

    df = pd.read_csv(file, sep='\t', header=0)
    df = df[np.isfinite(pd.to_numeric(df[pvalue]))].copy()
    df.sort_values(by=[pvalue, 'SCORE', 'CGC'], ascending=[True, False, False], inplace=True)

    if len(df) != 0:
        significant_pvalues = df.loc[df[pvalue] > alpha][pvalue]

        # Subset p-values > 0.1
        if len(significant_pvalues) > 1000:
            pvalues = np.random.choice(significant_pvalues, size=random_set, replace=False)
        else:
            pvalues = significant_pvalues

        # Calculate KS
        ks = stats.kstest(pvalues, 'uniform')

        # Calculate KS sign
        observed = len(df.loc[df[pvalue] >= alpha])
        expected = (1 - alpha) * len(df)
        if observed < expected:
            ks_statistic = ks.statistic
        else:
            ks_statistic = - ks.statistic

        return ks_statistic, ks.pvalue, len(pvalues)
    else:
        return float('nan'), float('nan'), float('nan')
    
def get_weight(i, weight):
    """
    Weight contributions according to ranking

    Args:
        i (int): the ith ranking
        weight (str): the type of weighting [log, normal]

    Returns:
        (float): the weight of the ith position

    """
    if weight == "log":
        return 1.0 / np.log2(i + 2)
    if weight == "normal":
        return 1.0 / i


def calculate_percentage_cgc(ranking, cgc_genes):
    """
    Calculate the percentage of CGC genes in the input list

    Args:
        ranking (list): the input list of the ranked genes
        cgc_genes (set): set of genes included in the Cancer Gene Census (CGC)

    Returns:
        (float): percentage of cgc genes in the list

    """

    n = float(sum([1.0 if gene in cgc_genes else 0.0 for gene in ranking]))
    return n / len(ranking)


def evaluate_enrichment_method(file, cgc_genes, pvalue, weight="log", ranking_limit=100):
    """

    Args:
        file (str): path to file with elements results
        cgc_genes (set): set of genes in the CGC
        pvalue (str): p-value column name in dataframe to calculate KS
        weight (str): normalization of the weight. [log,normal] default: log
        ranking_limit (int): limit to calculate the area under the curve. Default: 40

    Returns:
        (foat): The weighted area under the curve of the CGC enrichment

    """

    df = pd.read_csv(file, sep='\t', header=0)
    df = df[np.isfinite(pd.to_numeric(df[pvalue]))].copy()
    df.sort_values(by=[pvalue, 'SCORE', 'CGC'], ascending=[True, False, False], inplace=True)

    ranking = df['SYMBOL'].tolist()

    if len(ranking) > ranking_limit:
        ranking = ranking[0:ranking_limit]
    xticks = range(len(ranking))
    area = 0.0
    for i in xticks:
        weight_i = get_weight(i, weight)
        x_i = calculate_percentage_cgc(ranking[0:i + 1], cgc_genes)
        area += x_i * weight_i

    return area

# Specify the root directory containing subdirectories with "elements_results.txt" files
root_directory = "/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/results/"

# Specify the output file path
output_file_path = "/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/results/calculate_ks.txt"

# Specify the COSMIC genes file path
cosmic_path = '/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/cosmic_grch38_v96.txt'

# Read cosmic gene list as set
with open(cosmic_path, 'r') as file:
    # Read lines from the file and create a set
    cosmic_set = set(file.read().splitlines())

# Specify the p-value column name
pvalue_column_name = "P_ANALYTICAL"

# Iterate over subdirectories and their files
results = []
with open(output_file_path, "w") as output_file:
    output_file.write("Directory\tKS Statistic\tKS P-value\tNumber of P-values\tENRICH_LOG\n")  # Header

    for root, dirs, files in os.walk(root_directory):
        if "elements_results.txt" in files:
            file_path = os.path.join(root, "elements_results.txt_updated")

            # Calculate KS for each file
            ks_statistic, ks_pvalue, num_pvalues = calculate_ks(file_path, pvalue_column_name)

            # Calculate area under the curve for CGC enrichment
            area = evaluate_enrichment_method(file_path, cosmic_set, pvalue_column_name)

	        # Extract directory name from root
            dir_name = os.path.basename(root)

            # Write results to the output file
            output_file.write(f"{dir_name}\t{abs(ks_statistic)}\t{abs(ks_pvalue)}\t{num_pvalues}\t{abs(area)}\n")

            # Append results to the list
            results.append({"Directory": dir_name, "KS Statistic": ks_statistic, "KS P-value": ks_pvalue, "Number of P-values": num_pvalues, "ENRICH_LOG": area})


print("Results written to:", output_file_path)

# Create a DataFrame from the results
df_results = pd.DataFrame(results)

if len(df_results) != 0:
    df_results.sort_values(by=['KS Statistic'], ascending=[True], inplace=True)
    print(min(df_results['KS Statistic'].tolist()))
    df_results = df_results.reset_index(drop=True)
    selected_rows = df_results.loc[df_results['KS Statistic'] <= min(df_results['KS Statistic'].tolist()) + 0.1].copy()

    # Sort by CGC enrichment and get top 1
    selected_rows.sort_values(by=['ENRICH_LOG'], ascending=[False], inplace=True)
    selected_rows = selected_rows.reset_index(drop=True)
    selected_rows = selected_rows.iloc[0]
    

# Print or use the selected rows
print(selected_rows[["Directory", "KS Statistic", "KS P-value", "Number of P-values", "ENRICH_LOG"]])