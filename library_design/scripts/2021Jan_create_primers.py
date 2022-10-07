"""Script for creating gene assembly primers.
    
Primers containing an ambiguous codon at their center are tiled along a gene.

Jesse Bloom, 2013. 

Edited by Adam Dingens Nov 2015 to generate primers of differing lengths to all have a Tm of ~60C. Notes on edits below. 
This script first makes an ORIGINAL primer of specified length (default 37 bps). 
If the ORIGINAL primer has a Tm of greater than MaxPrimerTm, then nucleotides are trimmed off one by one (first 5', then 3', then 5' etc) until the Tm is less than MaxPrimerTm. Note that this could be over a degree less than the MaxPrimerTm. 
If the ORIGINAL primer has a Tm less than MinPrimerTm, then nucelotides are added one by one (first 3', then 5', then 3' etc) until the Tm is over MinPrimerTm. Note that this could be over a degree more than the MinPrimerTm
If the ORIGINAL primer has a Tm of less than MaxPrimerTm but greater than MinPrimerTm, it is not altered. 
The primers are constrained to be between MinPrimerlength and MaxPrimerLength bps long. The Tm of some MaxPrimerLength primers may not be > MinPrimerTemp, and the Tm of some MinPrimerLength primers may not be < MaxPrimerTm.

For command line arguments, run::

    python create_primers.py -h

The  Tm_NN command of the MeltingTemp Module of Biopython (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers. 
This calculation is based on nearest neighbor thermodynamics. nucelotides labeled N are given average values in the Tm calculation. 
It is possible to vary salt concentration and other addatives if needed.

Edited by Kate Crawford January 2021 to include options for `NNG` and `NNC`
ambiguous codons to simulate `NNS` mutagenesis as IDT oPools only allow for
`N` or `K` ambiguous nucelotides.
"""


import os
import sys
import math
import re
import argparse
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


def Parser():
    """Returns command line parser."""
    parser = argparse.ArgumentParser(
            description='Script by Adam Dingens and Jesse Bloom to design codon tiling primers with specific melting temperatures.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )

    parser.add_argument('sequencefile', help="the name of a file giving the sequence for which we are designing the primers. This file should only contain the sequence, and should not have any headers or other content. For the sequence, make the 5' and 3' ends that you do not want to mutate in lower case. Make the portion of the coding sequence that you want to tile with ambiguous codons in upper case. Typically, for example, the first site you mutate would be the initial start codon, so the first *ATG* would be the first upper case letter. You must have at least *(startprimerlength - 3) / 2* nucleotides in lower case at each end of the upper case sequence that you are mutating. This is because at least this much flanking sequence is needed to design primers of the indicated length; more sequence may be required if the primer at either end is extended beyond the startprimerlength.")
    parser.add_argument('primerprefix', help="prefix name to be added to each primer")
    parser.add_argument('firstcodon', type=int, help='number to assign to first codon in infile to mutagenize, used for primer naming')
    parser.add_argument('outfile', help='name of primer output file')
    parser.add_argument('--startprimerlength', type=int, help='starting primer length', default=37)
    parser.add_argument('--maxprimertm', type=float, help="Upper temperature limit for primers.", default=61)
    parser.add_argument('--minprimertm', type=float, help="Lower temperature limit for primers.", default=60)
    parser.add_argument('--minlength', type=int, help='Minimum primer length', default=25)
    parser.add_argument('--maxlength', type=int, help='Maximum primer length', default=51)
    parser.add_argument('--ambiguous_codon', choices={'NNN', 'NNS', 'NNK', 'NNC', 'NNG'},
                        default='NNN', help='What ambiguous codon to use'),
    parser.add_argument('--output', choices={'plates', 'opools'}, default='plates',
                        help='Format the final csv for ordering in plates or oligo pools.')

    return parser


def ReverseComplement(seq):
    """Returns reverse complement of sequence. Preserves upper/lower case."""
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n', 'S':'S', 's':'s', 'K':'M', 'k':'m'}
    rc = [d[nt] for nt in seq]
    rc.reverse()
    return ''.join(rc)


def CreateMutForOligosVarLength(seq, primerlength, prefix, firstcodon, maxprimertm, minprimertm, maxlength, minlength, ambiguous_codon):
    """Creates oligos to tile a gene and introduce ambiguous codon at each site.

    *seq* : sequence of the gene. The gene itself should be upper case. The
    flanking regions and start / stop codons should be lower case. 
    All upper case codons are randomized. The length of the lower
    case sequences at each end must be >= (primerlength - 3) / 2.0

    *primerlength* : length of primers. Must be an odd number, so that equal length
    flanking on each side.

    *prefix* : string prefix attached to primer names.

    *firstcodon* : number assigned to first codon in primer name.

    *ambiguous_codon* : ambiguous codon to use.

    Tiles primers across the gene in the forward direction. The primers
    are all of length primerlength with ambiguous codon at the middle codon.
    Note that only upper case letters are randomized.
    Primers are named as follows:

    "%s-for-mut%d" % (prefix, i) -> 5' tiling primers, where i = 2, 3, ...
    In other words, the indices cover codons 2 and up.

    Returns a list of all these primers as *(name, sequence)* 2-tuples.
    """
    n = len(seq)
    assert primerlength % 2 == 1, "primer length not odd"
    initial_flanklength = (primerlength - 3) // 2
    upperseq = ''.join([nt for nt in seq if nt.istitle()])
    assert upperseq in seq, "upper case nucleotides not substring"
    assert len(upperseq) % 3 == 0, "length of upper case not multiple of 3"
    startupper = seq.index(upperseq)
    if startupper < initial_flanklength:
        raise ValueError("not enough 5' lower case flanking nucleotides")
    if n - len(upperseq) - startupper < initial_flanklength:
        raise ValueError("not enough 3' lower case flanking nucleotides")
    ncodons = len(upperseq) // 3
    primers = []
    for icodon in range(ncodons):
        i = startupper + icodon * 3
        primer = "%s%s%s" % (seq[i - initial_flanklength : i], ambiguous_codon, seq[i + 3 : i + 3 + initial_flanklength]) 
        name = "%s-for-mut%d" % (prefix, firstcodon + icodon)
        primerseq = Seq(primer)
        
        TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
        add_3 = True
        minus_5 = True
        flank5 = flank3 = initial_flanklength
        if float(TmNN) > float(maxprimertm):
            while float(TmNN) > float(maxprimertm) and len(primer) > minlength:
                if minus_5:
                    flank5 -= 1
                    primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                    minus_5 = False
                else:
                    flank3 -= 1
                    primer =  "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])

                    minus_5 = True
                if startupper < flank5:
                    raise ValueError("not enough 5' lower case flanking nucleotides")
                if n - len(upperseq) - startupper < flank3:
                    raise ValueError("not enough 3' lower case flanking nucleotides") #not sure if this is correct!

                
                primerseq = Seq(primer)
                TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                primerlength= len(primer)
        else: 
            if float(TmNN) < float(minprimertm):
                while float(TmNN) < float(minprimertm) and len(primer) < maxlength:
                    if add_3:
                        flank3 += 1
                        primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                        add_3 = False
                    else:
                        flank5 +=1
                        primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                        add_3 = True    
                    primerseq = Seq(primer)
                    TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                    primerlength= len(primer)
                    if startupper < flank5:
                        raise ValueError("not enough 5' lower case flanking nucleotides")
                    if n - len(upperseq) - startupper < flank3:
                        raise ValueError("not enough 3' lower case flanking nucleotides") 

            else:
                pass
        primers.append((name, primer))
    #print primers 
    return primers



def main():
    parser = Parser()
    args = vars(parser.parse_args())

    print("Read the following command line arguments")
    for (argname, argvalue) in args.items():
        print("\t%s = %s" % (argname, argvalue))


    primerlength = args['startprimerlength']

    if (primerlength <=3 ) or (primerlength % 2 == 0):
        raise ValueError("Does not appear to be valid primer length: %d" % primerlength)
    
    sequencefile = args['sequencefile']
    if not os.path.isfile(sequencefile):
        raise IOError("Cannot find sequencefile %s" % sequencefile)
    sequence = open(sequencefile).read() 
    sequence = sequence.replace(' ', '')
    sequence = sequence.replace('\n', '')
    print("Read a sequence of length %d from %s:\n%s" % (len(sequence), sequencefile, sequence))
    outfile = args['outfile']
    primerprefix = args['primerprefix']
    firstcodon =  args['firstcodon']
    print("The primers will be named with the prefix %s, and the first codon numbered as %d." % (primerprefix, firstcodon))

    # Design forward mutation primers
    mutforprimers = CreateMutForOligosVarLength(sequence, primerlength, primerprefix, firstcodon, args['maxprimertm'], args['minprimertm'], args['maxlength'], args['minlength'], args['ambiguous_codon'])
    print("Designed %d mutation forward primers." % len(mutforprimers))


    
    # Design reverse mutation primers
    mutrevprimers = [(name.replace('for', 'rev'), ReverseComplement(seq)) for (name, seq) in mutforprimers]
    print("Designed %d mutation reverse primers." % len(mutrevprimers))
   
    # Print out all of the primers
    primers = mutforprimers + mutrevprimers
    print("This gives a total of %d primers." % len(primers))


    print("\nNow writing these primers to %s" % outfile)

    if args['output'] == 'plates':
        iplate = 1
        f = open(outfile, 'w')
        for primers in [mutforprimers, mutrevprimers]:
            f.write("\r\nPlate %d\r\n" % iplate)
            n_in_plate = 0
            for (name, primer) in primers:
                f.write("%s, %s\r\n" % (name, primer))
                n_in_plate += 1
                if n_in_plate == 96:
                    n_in_plate = 0
                    iplate += 1
                    f.write("\r\nPlate %d\r\n" % iplate)
            if n_in_plate:
                iplate += 1
    elif args['output'] == 'opools':
        f = open(outfile, 'w')
        f.write("Pool name,Primer name,Ambiguous codon,Sequence\r\n")
        for primers in [mutforprimers, mutrevprimers]:
            if primers == mutforprimers:
                pool = f"{primerprefix}_ForPool"
            elif primers == mutrevprimers:
                pool = f"{primerprefix}_RevPool"
            for (name, primer) in primers:
                f.write("%s,%s,%s,%s\r\n" % (pool, name, args['ambiguous_codon'], primer))
    else:
        raise ValueError(f"Invalid output of {args['output']}")



main() # run the main program
