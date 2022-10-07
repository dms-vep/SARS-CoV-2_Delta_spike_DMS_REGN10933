# Reference sequences
This subdirectory contains contains SARS-CoV-2 spike reference sequences and lookup tables required to renumber positions between variants.

## Files
- `pH2rU3_ForInd_sinobiological_617.2_Spiked21_CMV_ZsGT2APurR.gb` genbank file for plasmid used to generate spike DMS library

- `B16172_extended_ends_for_primers.txt` Delta variant spike sequence containing flanking plasmid sequences, copied from `pH2rU3_ForInd_sinobiological_617.2_Spiked21_CMV_ZsGT2APurR.gb` plasmid. Used in generating DMS primers.

- `B16172_extended_ends_for_paired_primers.fasta` fasta file of Delta variant spike sequence containing plasmid sequences, copied from `pH2rU3_ForInd_sinobiological_617.2_Spiked21_CMV_ZsGT2APurR.gb` plasmid. Used in generating paired DMS primers.

- `homo_codon_freq_del.csv` Human codon frequencies used in DMS primer design. Frequencies were obtained from [here](https://www.kazusa.or.jp/codon/)

- `reference_sequence_position_lookup.csv` table used to look up SARS-CoV-2 variant spike sequence positions in relation to Wuhan-1 spike sequence. Columns are as follows:
	- `parent_seq` Wuhan-1 spike amino acid sequence
	- `variant_seq` Delta variant spike amino acid sequence
	- `parent_pos` Wuhan-1 spike amino acid position
	- `variant_pos` Delta variant spike amino acid position
	- `variant_sig` Signiture amino acid changes in Delta variant. `Yes` for amino acid changes that are unique to Delta variant and `No` for amino acids that are the the same in relation to Wuhan-1 sequence. 



