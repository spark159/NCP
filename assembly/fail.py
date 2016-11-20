# read overlap information
def read_input (fname):
    result=[]
    for line in open(fname):
        i,j,ori,pos=line.strip().split()
        result.append([int(i),int(j),ori,int(pos)])
    return result
overlap_data=read_input("lab01.olaps")

# merge two clusters including i j elements
def merge (cluster_set, i,j, offset):
    cpair=[]
    for target in [i,j]:
        for k in range(len(cluster_set)):
            if target  in cluster_set[k]:
                cpair.append(cluster_set.pop(k))
                break
    if len(cpair) ==2:
        if offset <0: cpair=cpair[::-1]
        cluster_set.append(cpair[0]+cpair[1])
    else:
        cluster_set.append(cpair[0])    
    return None

#find contigs information
def find_contigs (overlap_data):
    contigs=[[i] for i in range(1,301)]
    overlap=overlap_data[:]
    
    while len(overlap)>0:
        # find closest pair 
        min_k=0; min_off=overlap[min_k][3]
        for k in range(len(overlap)):
            i,j,ori,pos=overlap[k]
            if abs(pos) < abs(min_off):
                min_k=k; min_off=pos; min_i=i; min_j=j
        # delete used data and merge two clusters
        overlap.pop(min_k); merge(contigs, min_i, min_j, min_off)

    return contigs

contigs=find_contigs(overlap_data)
print contigs
