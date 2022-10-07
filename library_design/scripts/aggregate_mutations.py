"""Implements ``snakemake`` rule to aggregate all mutations across primer sets."""


import re

import Bio.Data.IUPACData

import pandas as pd


# get variables from snakemake
paired_positiveSel_primers = snakemake.input.paired_positiveSel_primers
usher_primers = snakemake.input.usher_primers
positiveSel_primers = snakemake.input.positiveSel_primers
gisaid_primers = snakemake.input.gisaid_primers
output_csv = snakemake.output.csv

amino_acids = Bio.Data.IUPACData.protein_letters

# read paired positive selection sites
paired_positive = (
    pd.read_csv(paired_positiveSel_primers, header=None, names=['name', 'primer'])
    ['name'].str
    .extract(
        '^delta_PS_paired\-(?:for|rev)\-mut_epi\-(?P<site1>\d+)'
        '(?:\.0)?\-(?P<site2>\d+)(?:\.0)?(?:_NN[CG]){1,2}$'
    )
    .melt(value_vars=['site1', 'site2'], value_name='site')
    [['site']]
    .assign(
        mutation_type='paired positive selection',
        site=lambda x: x['site'].astype(int),
        amino_acid=amino_acids,
    )
    .drop_duplicates()
    .assign(amino_acid=lambda x: x['amino_acid'].map(list))
    .explode('amino_acid')
    [['site', 'amino_acid', 'mutation_type']]
)
assert paired_positive.notnull().all().all()

# read positive selection sites
positive = (
    pd.read_csv(positiveSel_primers)
    ['Primer name'].str
    .extract('^positiveSelectionSpike_NN[CG]\-(?:for|rev)\-mut(?P<site>\d+)$')
    .assign(
        site=lambda x: x['site'].astype(int),
        mutation_type='positive selection',
        amino_acid=amino_acids,
    )
    .drop_duplicates()
    .assign(amino_acid=lambda x: x['amino_acid'].map(list))
    .explode('amino_acid')
    [['site', 'amino_acid', 'mutation_type']]
)
assert positive.notnull().all().all()

# read UsHER recurrent mutations
recurrent = (
    pd.read_csv(usher_primers, header=None, names=['name', 'primer'])
    ['name'].str
    .extract('^variant_usher\-(?:for|rev)\-mut(?P<site>\d+)(?P<amino_acid>[A-Z])$')
    .assign(
        mutation_type='recurrent mutation',
        site=lambda x: x['site'].astype(int),
    )
    .drop_duplicates()
    [['site', 'amino_acid', 'mutation_type']]
)
assert recurrent.notnull().all().all()

# read GISAID observed mutations
observed = (
    pd.read_csv(gisaid_primers)
    ['primer_name'].str
    .extract('^variantGISAID\-(?:for|rev)\-mut(?P<site>\d+)(?P<amino_acid>[A-Z\-])$')
    .assign(
        mutation_type='observed mutation',
        site=lambda x: x['site'].astype(int),
    )
    .drop_duplicates()
    [['site', 'amino_acid', 'mutation_type']]
)
assert observed.notnull().all().all()

# concatenate everything and write to file
(pd.concat([paired_positive, positive, recurrent, observed])
 .sort_values(['mutation_type', 'site', 'amino_acid'])
 .rename(columns={"site": "sequential_site"})
 .to_csv(output_csv, index=False)
 )
