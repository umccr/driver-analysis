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

# Specify the root directory containing subdirectories with "elements_results.txt" files
root_directory = "/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/results/"

# Specify the output file path
output_file_path = "/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/results/calculate_ks.txt"

# Specify the p-value column name
pvalue_column_name = "P_ANALYTICAL"

# Iterate over subdirectories and their files
results = []
with open(output_file_path, "w") as output_file:
    output_file.write("Directory\tKS Statistic\tKS P-value\tNumber of P-values\n")  # Header

    for root, dirs, files in os.walk(root_directory):
        if "elements_results.txt" in files:
            file_path = os.path.join(root, "elements_results.txt")

            # Calculate KS for each file
            ks_statistic, ks_pvalue, num_pvalues = calculate_ks(file_path, pvalue_column_name)

	        # Extract directory name from root
            dir_name = os.path.basename(root)

            # Write results to the output file
            output_file.write(f"{dir_name}\t{abs(ks_statistic)}\t{abs(ks_pvalue)}\t{num_pvalues}\n")

            # Append results to the list
            results.append({"Directory": root, "KS Statistic": ks_statistic, "KS P-value": ks_pvalue, "Number of P-values": num_pvalues})


print("Results written to:", output_file_path)

# Create a DataFrame from the results
df_results = pd.DataFrame(results)

if len(df_results) != 0:
    df_results.sort_values(by=['KS Statistic'], ascending=[True], inplace=True)
    print(min(df_results['KS Statistic'].tolist()))
    df_results = df_results.reset_index(drop=True)
    selected_rows = df_results.loc[df_results['KS Statistic'] <= min(df_results['KS Statistic'].tolist()) + 0.1].copy()

# Print or use the selected rows
print(selected_rows[["Directory", "KS Statistic", "KS P-value", "Number of P-values"]])