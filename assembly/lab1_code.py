# read input file
def read_input (filename):
    seq_reads={}
    for line in open(filename):
        if line[0] == '>':
            key=int(line[1:-1])
            seq_reads[key]=''
        else:
            seq_reads[key] += line.strip()
    return seq_reads
seq_reads=read_input("lab01.fasta")

# make complementary sequence
def comp_seq (seq):
    comp_base={'a':'t','t':'a','c':'g','g':'c'}
    comp_seq=''
    for base in seq:
        comp_seq += comp_base[base]
    return comp_seq[::-1]

# give overlap information for two sequences
def overlap (seq1, seq2):
    N=len(seq1); comp_seq2=comp_seq(seq2)
    num=0; pos='null'; ori='null'
    for i in range(-(N-40),N-40+1):
        st1=max(0,i) ; en1=min(N,N+i)
        st2=max(0,-i); en2=min(N,N-i)
        if seq1[st1:en1] == seq2[st2:en2]:
            if en1-st1 > num:
                pos=i;ori='F';num=en1-st1
        if seq1[st1:en1] == comp_seq2[st2:en2]:
            if en1-st1 > num:
                pos=i;ori='R';num=en1-st1
    if num ==0:
        return False, pos, ori
    return True, pos,ori

# give overlap search output list
def pairwise(seq_reads):
    N=len(seq_reads); output=''
    result=[]
    for i in range(1,N):
        for j in range(i+1,N+1):
            I, pos, ori = overlap(seq_reads[i],seq_reads[j])
            if I:
                output += '%03d  %03d  %s  %d  \n' % (i,j,ori,pos)
                result.append([i,j,ori,pos])
    f=open('lab01.olaps','w')
    f.write(output)   
    return result

result=pairwise(seq_reads)
