#!/usr/bin/env python

# Nucleosome_positioning by Sangwoo Park and Daehwan Kim, September 2016
# Calculate probability of nucleosome positioning on given DNA sequence
# Yeast genome data Based on Pubmed, Saccharomyces cerevisiae (UCSC -SAC2)
#    http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/
# Positioning data Based on Kristin, et al, Nature 2012
# Algorithm Based on Segal, et al, Nature 2006

import sys, os
import random
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Import some of our own modules
from models import HMM

# read out yeast genome data
# output is a dictionary, e.g., {"chrI": ACCATG, "chrII": TGCTT}
def read_genome(fname):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic

# get complementary sequence
# e.g., ACG => CGT
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
                Afreq[i] += 1.0
            elif dint in ['CC','CG','GC','GG']:
                Gfreq[i] += 1.0
                
    return Afreq / len(NCP_seq_list), Gfreq / len(NCP_seq_list)

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
        markov_order,
        verbose):
    # read out yeast genome data
    ygenome = read_genome("data/scerevisiae.fa")

    # read out nucleosome positioning & sequence data
    def read_NCP_list(NCP_fname):
        NCP_list = []
        for line in open(NCP_fname):
            line = line.strip()
            try:
                chr, pos, score, noise= line.split()
                pos, score, noise = int(pos) - 1, float(score), float(noise)
                if noise <= 0.0:
                    continue
                NCP_list.append([chr, pos, score / noise])
            except ValueError:
                continue
        return NCP_list
    
    NCP_list = read_NCP_list("data/nature11142-s2.txt")

    def read_NCP_sequence(genome, NCP):
        chr, pos, score = NCP
        st = pos - nucleosome_dna_len / 2; en = pos + nucleosome_dna_len / 2
        if st < 0 or en >= len(genome[chr]):
            return ""
        seq = genome[chr][st:en+1]; #rev_comp_seq=get_comp(seq)[::-1]
        return seq
            
    def read_NCP_sequences(genome, NCP_list):
        NCP_seq_list = []
        for NCP in NCP_list:
            seq = read_NCP_sequence(genome, NCP)
            if seq == "":
                continue
            NCP_seq_list.append(seq);
        return NCP_seq_list

    NCP_seq_list = read_NCP_sequences(ygenome, NCP_list)

    # Randomly select 99% of NCPs as a training set and the rest, 1%, as a validation set
    NCP_random_list = NCP_list[:] # deep copy
    random.shuffle(NCP_random_list)
    num99 = len(NCP_random_list) * 99 / 100
    NCP_tlist, NCP_vlist = NCP_random_list[:num99], NCP_random_list[num99:]
    NCP_seq_tlist, NCP_seq_vlist = read_NCP_sequences(ygenome, NCP_tlist), read_NCP_sequences(ygenome, NCP_vlist)

    mm = HMM.MarkovModel(markov_order)
    mm.train(ygenome, NCP_seq_tlist)
    num_test, num_correct = 0, 0
    for NCP in NCP_vlist:
        chr, pos = NCP[:2]
        max_i, max_score = -1, -sys.float_info.max
        for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
            seq = read_NCP_sequence(ygenome, [chr, i, 0.0])
            cur_score = mm.predict(seq)
            
            # comp_seq = get_comp(seq)
            # cur_score = max(mm.predict(seq), mm.predict(comp_seq))
            if max_score < cur_score:
                max_i = i
                max_score = cur_score

        num_test += 1
        if pos == max_i:
            num_correct += 1

        # DK - for debugging purposes
        if pos != max_i and False:
            print NCP
            print "predicted: %d, score: %f" % (max_i, max_score)
            print "true: %d, score: %f" % (pos, mm.predict(read_NCP_sequence(ygenome, [chr, pos, 0.0])))
            for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
                seq = read_NCP_sequence(ygenome, [chr, i, 0.0])
                cur_score = mm.predict(seq)
                print "\t%f at %d" % (cur_score, i)

            chr_seq = ygenome[chr]
            chr_len = len(chr_seq)
            left = max(0, pos - 1000)
            right = min(chr_len, pos + 1000)
            seq = chr_seq[left:right]
            F = mm.get_forward(seq)
            R = mm.get_reverse(seq)

            print F
            # print R

            max_i, max_score = -1, -sys.float_info.max
            for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
                cur_score = mm.prob_i(seq, i - left, F, R)
                if max_score < cur_score:
                    max_i = i
                    max_score = cur_score

            print "Global: %d, score: %f" % (max_i, max_score)

        
            sys.exit(1)
            
    print "%d-order Markov Model: %.2f%% (%d/%d)" % (markov_order, float(num_correct)/num_test*100, num_correct, num_test)
    # mm.help()
    sys.exit(1)
    

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
    parser.add_argument('--markov-order',
                        dest='markov_order',
                        type=int,
                        default=1,
                        help='Markov Model order (default: 1)')    
    parser.add_argument('--seed',
                        dest='seed',
                        type=int,
                        default=1,
                        help='Random seeding value (default: 1)')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if args.nucleosome_dna_len % 2 == 0:
        print >> sys.stderr, "Error: please use an odd number for --nucleosome-dna-len, perhaps %d instead of %d" % (args.nucleosome_dna_len + 1, args.nucleosome_dna_len)
        sys.exit(1)

    random.seed(args.seed)

    NCP(args.nucleosome_dna_len,
        args.markov_order,
        args.verbose)
