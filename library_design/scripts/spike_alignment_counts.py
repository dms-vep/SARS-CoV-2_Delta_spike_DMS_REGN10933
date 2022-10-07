import collections

import Bio.SeqIO

import pandas as pd


print(f"Reading alignment from {snakemake.input.alignment}")
counts = collections.defaultdict(lambda: collections.defaultdict(int))
n = 0
for seq in Bio.SeqIO.parse(snakemake.input.alignment, 'fasta'):
    n += 1
    for site, aa in enumerate(seq.upper(), 1):
        if aa != 'X':
            counts[site][aa] += 1
print(f"Counted amino acids at each site for {n} sequences")

df = (pd.DataFrame.from_dict(counts, orient='index')
      .rename_axis('site')
      .reset_index()
      .melt(id_vars='site',
            var_name='amino_acid',
            value_name='alignment_counts',
            )
      .query('alignment_counts.notnull()')
      .assign(alignment_counts=lambda x: x['alignment_counts'].astype(int))
      )

print(f"Writing counts for all {len(df)} amino-acids with non-zero counts "
      f"to {snakemake.output.alignment_counts}")
df.to_csv(snakemake.output.alignment_counts, index=False)
