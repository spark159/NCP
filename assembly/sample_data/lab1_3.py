
# open fasta file
def fasta_input (filename):
    seq_reads={}
    for line in open(filename):
        if line[0] == '>':
            key=int(line[1:-1])
            seq_reads[key]=''
        else:
            seq_reads[key] += line.strip()
    return seq_reads
seq_info=fasta_input("sample.fasta")

# open unitig file
def unis_input (file_name):
    unis_dic={}
    for line in open(file_name):
        linelist=line.strip().split()
        if linelist[0]=="UNI":
            uid=int(linelist[1])
            unis_dic[uid]=[]
        else:
            num,ori,pos=linelist
            unis_dic[uid].append([int(num),ori,int(pos)])
    return unis_dic
unis_info=unis_input("sample.unis")

# make complementary sequence
def comp_seq (seq):
    comp_base={'a':'t','t':'a','c':'g','g':'c'}
    comp_seq=''
    for base in seq:
        comp_seq += comp_base[base]
    return comp_seq[::-1]
                            
# generate unitig sequence
def unis_syn (seq_info, unis_info):
    unis_seq={}
    for uid in unis_info:
        seq=''
        i=0; n,o,p = unis_info[uid][i]
        while i < len(unis_info[uid]):
            target=seq_info[n]
            if o =='R':
                target=comp_seq(target)
            i+=1;
            if i==len(unis_info[uid]):
                np=len(target)
            else:
                nn, no, np = unis_info[uid][i]
            seq += target[0:np]
            n,o,p = nn,no,np
        unis_seq[uid]=seq
    return unis_seq
unis_seq=unis_syn (seq_info, unis_info)

