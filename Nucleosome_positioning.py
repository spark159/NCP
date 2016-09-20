#!/usr/bin/env python

# Nucleosome_positioning By Sangwoo Park and Daehwan Kim, September 2016
# Calculate probability of nucleosome positioning on given DNA sequence
# Yeast genome data Based on Pubmed, Saccharomyces cerevisiae S288c (assembly R64)
# Positioning data Based on Kristin, et al, Nature 2012
# Algorithm Based on Segal, et al, Nature 2006

import sys, os
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# read out yeast genome data
def read_genome(fname):
    genome = {}
    for line in open(fname, "rU"):
        if '>' in line:
            key=line.strip('>').strip('\n')
            s=''
        else:
            s += line.rstrip()
            genome[key] = s
    return genome

# get complementary sequence 
def get_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for base in rev_seq:
        if base == 'A':
            new_seq += 'T'
        if base == 'T':
            new_seq += 'A'
        if base == 'C':
            new_seq += 'G'
        if base == 'G':
            new_seq += 'C'
    return new_seq
     
# dinucleotide step analysis (A vs. G)
def AG_freq (NCP_seq_list, nucleosome_dna_len):
    Afreq=np.zeros(nucleosome_dna_len - 1); Gfreq=np.zeros(nucleosome_dna_len - 1)
    for seq in NCP_seq_list:
        for i in range(len(seq)-1):
            dint = seq[i:i+2]
            if dint in ['AA','AT','TA','TT']:
                Afreq[i] += 1.0/len(NCP_seq_list)
            elif dint in ['CC','CG','GC','GG']:
                Gfreq[i] += 1.0/len(NCP_seq_list)
    return Afreq, Gfreq

# dinucleotide step analysis (full)
def din_freq (NCP_seq_list, nucleosome_dna_len):
    freq = {}
    for seq in NCP_seq_list:
        for i in range(len(seq)-1):
            dint = seq[i:i+2]
            if dint not in freq:
                freq[dint] = np.zeros(nucleosome_dna_len)
            freq[dint][i] += 1.0 / len(NCP_seq_list)
    return freq

# get transition matrix of NCP positining
def get_NCP_m (NCP_seq_list, nucleosome_dna_len):    
    freq=din_freq(NCP_seq_list, nucleosome_dna_len);
    result=[]; base=["A", "T", "C", "G"]
    for i in range(nucleosome_dna_len - 1):
        matrix=np.zeros([4,4])
        for j in range(4):
            norm=0.0; row=np.zeros(4)
            for k in range(4):
                key=base[j]+base[k]
                row[k]= freq[key][i]; norm += freq[key][i]
            matrix[j]=row/float(norm)
        result.append(matrix)
    return result

# calculate probabilty of NCP positioning on given 147 bp sequence
def get_NCP_prob(seq, nucleosome_dna_len):
    if len(seq) != nucleosome_dna_len: return None
    global tm
    prob=0.25; dic={'A':0, 'T':1, 'C':2, 'G':3}
    for i in range(146):
        n=dic[seq[i]]; m=dic[seq[i+1]]
        prob *= tm[i][n][m]
    return prob

# calculate the sum of weight-factor of given sub sequence (forward)
def get_forward (seq, n, nucleosome_dna_len):
    F=[1.0]
    for i in range(1,n+1):
        if i < nucleosome_dna_len:
            F.append(F[-1]*1)
        else:
            F.append(F[-1]*1 + F[-nucleosome_dna_len]*get_NCP_prob(seq[i-nucleosome_dna_len:(i-1)+1]))
    return F

# calculate the sum of weight-factor of given sub sequence (reverse)
def get_reverse (seq, n, nucleosome_dna_len):
    R=[1.0]; N=len(seq)
    for i in range(1,n+1):
        if i < nucleosome_dna_len:
            R.append(R[-1]*1)
        else:
            R.append(R[-1]*1 + R[-nucleosome_dna_len]*get_NCP_prob(seq[N-i:N-i+nucleosome_dna_len]))
    return R

# calculate the probabilty of NCP would be on i th position of given sequence
def NCPprob (seq, i, nucleosome_dna_len):
    if len(seq[i-1:]) < nucleosome_dna_len:
        return 0.0
    return (get_forward(seq, i-1)[-1] * get_NCP_prob (seq[i-1:i+nucleosome_dna_len-1]) * get_reverse(seq, len(seq)-(i+nucleosome_dna_len-1))[-1]) / get_reverse(seq, len(seq))[-1]

# get NCP positioning probabiliy profile for a long sequence
def NCPprob_profile (seq, nucleosome_dna_len):
    F=get_forward(seq, len(seq)); R=get_reverse(seq, len(seq))
    profile=[]
    for i in range(1, len(seq)+1):
        if len(seq[i-1:]) < nucleosome_dna_len:
            prob=0.0;
        else:
            prob=(F[i-1] * get_NCP_prob (seq[i-1:i+nucleosome_dna_len-1]) * R[len(seq)-(i+nucleosome_dna_len-1)]) / R[-1]
        profile.append(prob)
    return profile

# get NCP occupancy profile of given sequence
def NCPoccupancy (profile, nucleosome_dna_len):
    result=[]
    for i in range(len(profile)):
        prob=0.0
        for j in range(nucleosome_dna_len):
            pos=i-j
            if pos <0:
                break
            else:
                prob += profile[pos]
        result.append(prob)
    return result


def NCP(nucleosome_dna_len,
        verbose):
    # read out yeast genome data
    ygenome = read_genome("data/scerevisiae.fa")

    # read out nucleosome positioning & sequence data
    NCP_list, NCP_seq_list = [], []
    for line in open("data/nature11142-s2.txt"):
        line = line.strip()
        try:
            chr, pos, score, noise= line.split()
            pos, score, noise = int(pos) - 1, float(score), float(noise)
            if noise <= 0.0:
                continue
            NCP_list.append([chr, pos, score / noise])
        except ValueError:
            continue

    for NCP in NCP_list:
        chr, pos, score = NCP
        st = pos - nucleosome_dna_len / 2; en = pos + nucleosome_dna_len / 2
        if st < 0 or en >= len(ygenome[chr]):
            continue
        seq = ygenome[chr][st:en+1]; #rev_comp_seq=get_comp(seq)[::-1]
        NCP_seq_list.append(seq); #NCP_seq.append(rev_comp_seq)

    tm = get_NCP_m(NCP_seq_list, nucleosome_dna_len)

    #freq=din_freq(NCP_seq)
    #newAfreq=freq['AA']+freq['AT']+freq['TA']+freq['TT']
    Afreq, Gfreq = AG_freq(NCP_seq_list, nucleosome_dna_len)
    #print Afreq
    plt.figure(); plt.plot(Afreq); plt.show()
    #plt.figure(); plt.plot(Gfreq)
    #profile=NCPprob_profile (ygenome["chrII"])
    #occupancy=NCPoccupancy (profile)
    #plt.figure(); plt.plot(profile)
    #plt.figure(); plt.plot(occupancy)
    #print ygenome['chrI']


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Nucleosome center positioning')
    parser.add_argument('--nucleosome-dna-len',
                        dest='nucleosome_dna_len',
                        type=int,
                        default=147,
                        help='Flanking sequence length (both sides including the center, default: 147')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if args.nucleosome_dna_len % 2 == 0:
        print >> sys.stderr, "Error: please use an odd number for --nucleosome-dna-len, perhaps %d instead of %d" % (args.nucleosome_dna_len + 1, args.nucleosome_dna_len)
        sys.exit(1)

    NCP(args.nucleosome_dna_len,
        args.verbose)
