import io
import subprocess
import tempfile

import Bio.SeqIO


# get variables from snakemake
allprots_file = snakemake.input.allprots
alignment_file = snakemake.output.alignment
ref_name = snakemake.params.ref_name
length_range = snakemake.params.length_range
max_ambiguous = snakemake.params.max_ambiguous
threads = snakemake.threads

print(f"Reading proteins from {allprots_file}")
prots = list(Bio.SeqIO.parse(allprots_file, 'fasta'))
print(f"Read {len(prots)=} proteins")

print(f"Finding protein with name that includes {ref_name}")
refprot = [p for p in prots if ref_name in p.id]
assert len(refprot) == 1, 'did not find unique reference'
refprot = refprot[0]
refprot_length = len(refprot)
print(f"Reference protein has length {refprot_length}")

# filter sequences
retained = []
n_invalid_length = n_too_ambiguous = 0
for prot in prots:
    if not (refprot_length - length_range <= len(prot) <= refprot_length + length_range):
        n_invalid_length += 1
    elif prot.seq.count('X') + prot.seq.count('x') > max_ambiguous:
        n_too_ambiguous += 1
    else:
        retained.append(prot)
prots = retained
print(f"Retained {len(prots)=} sequences")
print(f"Removed {n_invalid_length=} for having lengths >{length_range} different from reference")
print(f"Removed {n_too_ambiguous=} for having >{max_ambiguous} ambiguous residues.")

# align sequences
chunksize = 250000  # align in chunks of this many sequences
alignment = []
for chunk_start in range(0, len(prots), chunksize):
    chunk_end = min(len(prots), chunk_start + chunksize)
    print(f"Aligning proteins from {chunk_start=} to {chunk_end=}")
    with tempfile.NamedTemporaryFile(mode='w') as chunk_file, tempfile.NamedTemporaryFile(mode='w') as refprot_file:
        Bio.SeqIO.write([refprot], refprot_file, 'fasta')
        refprot_file.flush()
        Bio.SeqIO.write(prots[chunk_start: chunk_end], chunk_file, 'fasta')
        chunk_file.flush()
        mafft_cmds = ['mafft', '--auto', '--thread', str(threads),
                      '--keeplength', '--addfragments',
                      chunk_file.name, refprot_file.name]
        res = subprocess.run(mafft_cmds, capture_output=True)
        if res.returncode:
            raise RuntimeError(f"Alignment error:\n{res.stderr=}")
        else:
            with io.StringIO(res.stdout.decode('utf-8')) as f:
                chunk_alignment = list(Bio.SeqIO.parse(f, 'fasta'))
                # remove reference sequence, which should be first in file
                assert chunk_alignment[0].id == refprot.id, f"{chunk_alignment[0].id=}, {refprot.id=}"
                chunk_alignment = chunk_alignment[1:]
                assert len(chunk_alignment) == chunk_end - chunk_start
                alignment += chunk_alignment
print(f"Aligned {len(alignment)} sequences")
print(f"Writing alignment to {alignment_file}")
_ = Bio.SeqIO.write(alignment, alignment_file, 'fasta')

