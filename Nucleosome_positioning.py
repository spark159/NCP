# Nucleosome_positioning By Sangwoo Park, September 2016
# Calculate probability of nucleosome positioning on given DNA sequence
# Yeast genome data Based on Pubmed, Saccharomyces cerevisiae S288c (assembly R64)
# Positioning data Based on Kristin, et al, Nature 2012
# Algorithm Based on Segal, et al, Nature 2006

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# read out yeast genome data
ygenome={}
f = open("data/scerevisiae.fa", "rU")
for line in f:
    if '>' in line:
        key=line.strip('>').strip('\n')
        s=''
    else:
        s += line.rstrip()
        ygenome[key] = s
f.close()

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

# read out nucleosome positioning & sequence data
NCP_pos=[]; NCP_seq=[]
f=open("data/nature11142-s2.txt")
for line in f:
    if len(line.split()) ==4:
        NCP_pos.append(line.split())
for i in range(len(NCP_pos)):
    key=NCP_pos[i][0]; pos=int(NCP_pos[i][1]); q=float(NCP_pos[i][3])
    st=(pos-1)-73; en=(pos-1)+73
    #st=(pos-1); en=(pos-1) + 146
    if q > 0 and st >=0 and en <= len(ygenome[key])-1:
        seq=ygenome[key][st:en+1]; #rev_comp_seq=get_comp(seq)[::-1]
        NCP_seq.append(seq); #NCP_seq.append(rev_comp_seq)
f.close()

# dinucleotide step analysis (A vs. G)
def AG_freq (NCP_seq):
    Afreq=np.zeros(146); Gfreq=np.zeros(146)
    for seq in NCP_seq:
        for i in range(len(seq)-1):
            if seq[i:i+2] in ['AA','AT','TA','TT']:
                Afreq[i] += 1.0/len(NCP_seq)
            if seq[i:i+2] in ['CC','CG','GC','GG']:
                Gfreq[i] += 1.0/len(NCP_seq)
    return Afreq, Gfreq

# dinucleotide step analysis (full)
def din_freq (NCP_seq):
    freq={'AA':np.zeros(146), 'TT':np.zeros(146), 'AC':np.zeros(146), 'GT':np.zeros(146),'AG':np.zeros(146), 'CT':np.zeros(146), 'AT':np.zeros(146), 'CA':np.zeros(146), 'TG':np.zeros(146), 'CC':np.zeros(146), 'GG':np.zeros(146), 'CG':np.zeros(146), 'GA':np.zeros(146), 'TC':np.zeros(146), 'GC':np.zeros(146), 'TA':np.zeros(146)}
    for seq in NCP_seq:
        for i in range(len(seq)-1):
            freq[seq[i:i+2]][i] += 1.0 / len(NCP_seq)
    return freq

# get transition matrix of NCP positining
def get_NCP_m (NCP_seq):    
    freq=din_freq(NCP_seq); 
    result=[]; base=["A", "T", "C", "G"]
    for i in range(146):
        matrix=np.zeros([4,4])
        for j in range(4):
            norm=0.0; row=np.zeros(4)
            for k in range(4):
                key=base[j]+base[k]
                row[k]= freq[key][i]; norm += freq[key][i]
            matrix[j]=row/float(norm)
        result.append(matrix)
    return result

tm=get_NCP_m(NCP_seq)

# calculate probabilty of NCP positioning on given 147 bp sequence
def get_NCP_prob(seq):
    if len(seq) != 147: return None
    global tm
    prob=0.25; dic={'A':0, 'T':1, 'C':2, 'G':3}
    for i in range(146):
        n=dic[seq[i]]; m=dic[seq[i+1]]
        prob *= tm[i][n][m]
    return prob

# calculate the sum of weight-factor of given sub sequence (forward)
def get_forward (seq, n):
    F=[1.0]
    for i in range(1,n+1):
        if i <= 146:
            F.append(F[-1]*1)
        else:
            F.append(F[-1]*1 + F[-147]*get_NCP_prob(seq[i-147:(i-1)+1]))
    return F

# calculate the sum of weight-factor of given sub sequence (reverse)
def get_reverse (seq, n):
    R=[1.0]; N=len(seq)
    for i in range(1,n+1):
        if i <= 146:
            R.append(R[-1]*1)
        else:
            R.append(R[-1]*1 + R[-147]*get_NCP_prob(seq[N-i:(N-i+146)+1]))
    return R

# calculate the probabilty of NCP would be on i th position of given sequence
def NCPprob (seq, i):
    if len(seq[i-1:]) < 147:
        return 0.0
    return (get_forward(seq, i-1)[-1] * get_NCP_prob (seq[i-1:(i+145)+1]) * get_reverse(seq, len(seq)-(i+146))[-1])/ get_reverse(seq, len(seq))[-1]

# get NCP positioning probabiliy profile for a long sequence
def NCPprob_profile (seq):
    F=get_forward(seq, len(seq)); R=get_reverse(seq, len(seq))
    profile=[]
    for i in range(1, len(seq)+1):
        if len(seq[i-1:]) < 147:
            prob=0.0;
        else:
            prob=(F[i-1] * get_NCP_prob (seq[i-1:(i+145)+1]) * R[len(seq)-(i+146)])/ R[-1]
        profile.append(prob)
    return profile

# get NCP occupancy profile of given sequence
def NCPoccupancy (profile):
    result=[]
    for i in range(len(profile)):
        prob=0.0
        for j in range(147):
            pos=i-j
            if pos <0:
                break
            else:
                prob += profile[pos]
        result.append(prob)
    return result

#freq=din_freq(NCP_seq)
#newAfreq=freq['AA']+freq['AT']+freq['TA']+freq['TT']
#[Afreq,Gfreq]=AG_freq (NCP_seq)
#print Afreq
#plt.figure(); plt.plot(Afreq); plt.plot(Gfreq); plt.show()
#plt.figure(); plt.plot(Gfreq); plt.show()
profile=NCPprob_profile (ygenome["chrII"])
occupancy=NCPoccupancy (profile)
#plt.figure(); plt.plot(profile); plt.show()
plt.figure(); plt.plot(occupancy); plt.show()
