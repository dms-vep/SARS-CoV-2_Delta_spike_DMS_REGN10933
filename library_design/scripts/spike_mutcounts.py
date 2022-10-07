import Bio.SeqIO
import Bio.Seq

import pandas as pd


# get variables from `snakemake`
refseq_file = snakemake.input.refseq
genomic_mutcounts_tsv = snakemake.input.genomic_mutcounts
spike_mutcounts_csv = snakemake.output.spike_mutcounts

# read reference sequence and get spike
refseq = Bio.SeqIO.read(refseq_file, 'gb')
spike_feature = [f for f in refseq.features
                 if f.type == 'CDS' and f.qualifiers['gene'] == ['S']]
if len(spike_feature) != 1:
    raise ValueError(f"found {len(spike_feature)} spikes:\n{spike_feature}")
else:
    spike_feature = spike_feature[0]
spike = spike_feature.extract(refseq).seq
print(f"Extracted spike of {len(spike)} nucleotides.")
spikeprot = spike.translate(cds=True)
print(f"Translated to {len(spikeprot)} amino acids.")
spike_start = spike_feature.location.start + 1
spike_end = spike_start + len(spikeprot) * 3 - 1
print(f"Spike starts at nt {spike_start} and goes to {spike_end}")
spike = str(spike)
spikeprot = str(spikeprot)

def aa_mut(mutation):
    wt = mutation[0]
    site = int(mutation[1: -1])
    mut = mutation[-1]
    if spike_start <= site <= spike_end:
        spike_site = site - spike_start
        spike_icodon = spike_site // 3
        spike_codon = spike[3 * spike_icodon: 3 * spike_icodon + 3]
        codon_position = spike_site % 3
        wt_codon = spike_codon[: codon_position] + wt + spike_codon[codon_position + 1:]
        mut_codon = spike_codon[: codon_position] + mut + spike_codon[codon_position + 1:]
        wt_aa = str(Bio.Seq.Seq(wt_codon).translate())
        mut_aa = str(Bio.Seq.Seq(mut_codon).translate())
        if wt_aa != mut_aa:
            return f"{wt_aa}{spike_icodon + 1}{mut_aa}"
        else:
            return None
    else:
        return None

print(f"Reading genomic mutcounts from {genomic_mutcounts_tsv}")
genomic_mutcounts = pd.read_csv(genomic_mutcounts_tsv, sep='\t')
print(f"Read counts for {len(genomic_mutcounts)} mutations")

spike_mutcounts = (
    genomic_mutcounts
    .assign(spike_mutation=lambda x: x['ID'].map(aa_mut))
    .query('spike_mutation.notnull()')
    .assign(site=lambda x: x['spike_mutation'].str[1: -1].astype(int),
            amino_acid=lambda x: x['spike_mutation'].str[-1],
            )
    .query('amino_acid != "*"')
    .rename(columns={'occurrence': 'n_mutations_to'})
    [['site', 'amino_acid', 'n_mutations_to']]
    .sort_values('n_mutations_to', ascending=False)
    )

print(f"These genomic mutations correspond to {len(spike_mutcounts)} spike mutations")
print(f"Writing to {spike_mutcounts_csv}")
spike_mutcounts.to_csv(spike_mutcounts_csv, index=False)
