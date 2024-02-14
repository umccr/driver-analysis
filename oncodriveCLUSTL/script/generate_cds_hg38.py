import os, sys
#os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
from functools import partial

import warnings
warnings.simplefilter("ignore") # Change the filter in this process
os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

import pandas as pd
import pybedtools
from pybedtools import bedtool


def generate_gene_merged_coordinates(group):
    """Merge transcripts of a gene

    :param group: pandas.Dataframe

    """

    # Collapse overlapping transcript of the same gene
    group_bed = pybedtools.BedTool.from_dataframe(group)
    group_genemerged_bed = group_bed.merge(c=[5, 4, 6, 7], o=['distinct', 'distinct', 'collapse', 'distinct'])

    # Reformat
    group_genemerged_df = pd.read_csv(group_genemerged_bed.fn, sep='\t', header=None)
    group_genemerged_df = group_genemerged_df[[0, 1, 2, 4, 3, 5, 6]].copy()
    group_genemerged_df.sort_values(by=[0, 1, 2], axis=0, ascending=True, inplace=True)
    group_genemerged_df.columns = ['CHROMOSOME', 'START', 'END', 'STRAND', 'GENE_ID', 'TRANSCRIPT_ID', 'SYMBOL']

    return group_genemerged_df


# Read gencode annotation
gencode_annotations = '/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/gencode.v39.annotation.gtf.gz'
gencode = pd.read_csv(
    gencode_annotations,
    compression='gzip',
    sep="\t",
    header=None,
    skiprows=5,
    usecols=[0, 2, 3, 4, 6, 8],
    names=['CHROMOSOME', 'TYPE', 'START', 'END', 'STRAND', 'INFOS']
)

# Remove the `chr` prefix from the CHROMOSOME column
gencode['CHROMOSOME'] = gencode['CHROMOSOME'].map(lambda x: x[3:])

# Filter chromosomes (1-23, X, Y and M)
# Filtering the original DataFrame (gencode) to include only those rows where the 'CHROMOSOME_FILTER' column has the value 'PASS'.
chromosomes = set(map(str, list(range(1, 23)) + ['X', 'Y', 'M']))
gencode['CHROMOSOME_FILTER'] = gencode.apply(lambda x: 'PASS' if x['CHROMOSOME'] in chromosomes else 'FAIL', axis=1)
gencode = gencode.loc[gencode['CHROMOSOME_FILTER'] == 'PASS'].copy()
gencode.drop(['CHROMOSOME_FILTER'], axis=1)

# Parse the INFOS column
gencode['INFOS'] = gencode['INFOS'].map(
    lambda line: {
        i.split()[0]: i.split()[1].replace('"', '') for i in line.split(";") if len(i.split()) > 1
    }
)

# Get gene ID
gencode['GENE_ID'] = gencode['INFOS'].map(lambda x: x.get('gene_id', 'NaN.'))
# Remove gene ID version (eg 'ENSG00000223972.5_2' -> 'ENSG00000223972')
# If 'PAR_Y' (gene in pseudoautosomal chrY region) in gene ID, keep it so that it can be distinguished from chrX gene
gencode['GENE_ID'] = gencode.apply(lambda x: x['GENE_ID'].split('.')[0] if 'PAR_Y' not in x['GENE_ID'] else f"{x['GENE_ID'].split('.')[0]}_PAR_Y", axis=1)

# Get transcript ID
gencode['TRANSCRIPT_ID'] = gencode['INFOS'].map(lambda x: x.get('transcript_id', 'NaN.'))
# Remove transcript ID version (eg 'ENST00000223972.5_2' -> 'ENST00000223972')
# If 'PAR_Y' (gene in pseudoautosomal chrY region) in gene ID, keep it so that it can be distinguished from chrX gene
gencode['TRANSCRIPT_ID'] = gencode.apply(lambda x: x['TRANSCRIPT_ID'].split('.')[0] if 'PAR_Y' not in x['TRANSCRIPT_ID'] else f"{x['TRANSCRIPT_ID'].split('.')[0]}_PAR_Y", axis=1)

# Get gene type
gencode['GENE_TYPE'] = gencode['INFOS'].map(lambda x: x.get('gene_type', 'NaN'))

# Get symbol
gencode['SYMBOL'] = gencode['INFOS'].map(lambda x: x.get('gene_name', 'NaN'))

# Get transcript type
gencode['TRANSCRIPT_TYPE'] = gencode['INFOS'].map(lambda x: x.get('transcript_type', 'NaN'))

# Subset protein coding entries
gencode = gencode[
    (gencode['GENE_TYPE'] == 'protein_coding') &
    (gencode['TRANSCRIPT_TYPE'] == 'protein_coding')
].copy()

# Extract CDS
cds_df = gencode[(gencode['TYPE'] == 'CDS')].copy()
cds_df = cds_df[['CHROMOSOME', 'START', 'END', 'STRAND', 'GENE_ID', 'TRANSCRIPT_ID', 'SYMBOL']].copy()
cds_df.sort_values(by=['CHROMOSOME', 'START', 'END'], axis=0, ascending=True, inplace=True)

# Merge overlaping exons of the same gene
output_file = '/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/output_regions.tsv.gz'
transcripts_by_gene = cds_df.groupby('GENE_ID')
#import subprocess
#subprocess.run(['which', 'mergeBed'])
transcripts_gene_merged = transcripts_by_gene.apply(partial(generate_gene_merged_coordinates))
transcripts_gene_merged.reset_index(drop=True, inplace=True)
transcripts_gene_merged.sort_values(by=['CHROMOSOME', 'START', 'END'], axis=0, ascending=True, inplace=True)
transcripts_gene_merged.to_csv(output_file, sep="\t", compression='gzip', index=None, header=True)
