import numpy

import pandas as pd


print(f"Reading genomic selected sites from {snakemake.input.csv}")

p_cutoff = snakemake.params.p_cutoff
min_windows = snakemake.params.min_windows
print(f"Finding sites with P-value < {p_cutoff} in at least {min_windows} time windows")
df = (
    pd.read_csv(snakemake.input.csv)
    .query('gene == "S"')
    .query('beta > alpha')  # positively selected
    .query('p < @p_cutoff')
    .groupby('site', as_index=False)
    .aggregate(p=pd.NamedAgg('p', 'max'),  # get max selection at any timepoint
               n_time_windows=pd.NamedAgg('from', 'count'),
               )
    .query('n_time_windows > @min_windows')
    .sort_values('p')
    [['site', 'p', 'n_time_windows']]
    )
print(f"Found {len(df)} positively selected sites")

print(f"Writing sites to {snakemake.output.csv}")
df.to_csv(snakemake.output.csv, index=False)
