import copy

# read overlap information
def read_input (fname):
    result=[]
    for line in open(fname):
        i,j,ori,pos=line.strip().split()
        result.append([int(i),int(j),ori,int(pos)])
    return result
olaps_info=read_input("sample.olaps")

# get connectivity information
def get_web (olaps_info):
    web={}
    for olaps in olaps_info:
        i,j,ori,pos=olaps
        if not i in web:
            web[i]=[j]
        elif i in web:
            web[i].append(j)
        if not j in web:
            web[j]=[i]
        elif j in web:
            web[j].append(i)
    return web

# get subset of overlap information
def get_sublaps(web_info):
    global olaps_info
    result=[]
    for olaps in olaps_info:
        i,j,ori,pos=olaps
        try:
            if (j in web_info[i]) and (i in web_info[j]):
                result.append(olaps)
        except KeyError:
            continue
    return result
                                                
# find strand info for each reads
def get_strand (web_info,ref=None):
    info=get_sublaps(web_info)
    if ref:
        strand={ref[0]:ref[1]}
    else:
        strand={info[0][0]:1}
    while info:
        remain=[]
        for olaps in info:
            i,j,ori,pos=olaps
            if ori=='F':
                offset=1
            else:
                offset=-1
            if (i in strand) and (not (j in strand)):
                strand[j]=strand[i]*offset
            elif (not (i in strand)) and (j in strand):
                strand[i]=strand[j]*offset
            elif (i in strand) and (j in strand):
                continue
            else:
                remain.append(olaps)
        info=remain[:]
    return strand

# get directed graph
def get_flow (web_info, ref=None):
    olaps_info=get_sublaps(web_info)
    strand_info=get_strand(web_info, ref)
    result={}
    for olaps in olaps_info:
        i,j,ori,pos=olaps

        if strand_info[i]*pos >0:
            st=i; en=j
        else:
            st=j; en=i;
        pos=abs(pos)
        
        if not st in result:
            result[st]={'in':{}, 'out':{}}
        if not en in result:
            result[en]={'in':{}, 'out':{}}
        result[st]['out'][en]=pos
        result[en]['in'][st]=-pos
    return result
                                                                                                                                           
# check edge reduction of closed triple node
def check_redu (w3flow):
    a='NA'; b='NA'; c='NA'
    for node in w3flow:
        if len(w3flow[node]['out']) >1:
            a=node
        elif len(w3flow[node]['in']) >1:
            c=node
        else:
            b=node
    if a == 'NA' and b == 'NA':
        return False
    e1=w3flow[a]['out'][b]; e2=w3flow[b]['out'][c]; e3=w3flow[a]['out'][c]
    if e3 == e1+e2:
        return (a,c)
    return False

# pick one redundant edge 
def pick_redge (web_info):
    for node1 in web_info:
        if len(web_info[node1]) > 1:
            for i in range(len(web_info[node1])-1):
                for j in range(i+1,len(web_info[node1])):
                    node2=web_info[node1][i];node3=web_info[node1][j]
                    if node3 in web_info[node2]:
                        w3={node1:[node2,node3],node2:[node1,node3],node3:[node1,node2]}
                        w3flow=get_flow(w3)
                        edge=check_redu(w3flow)
                        if edge:
                            return edge
    return None
                            
# remove redundant edge connecting a pair of nodes
def remove_edge(web_info,pair):
    for i in range(len(pair)):
        for j in range(len(web_info[pair[i]])):
            if web_info[pair[i]][j] == pair[i-1]:
                dpoint=j; break
        del web_info[pair[i]][dpoint]
    return
                               
# get reduced web by removing redundant edges
def get_rweb (web_info):
    web=copy.deepcopy(web_info)
    while True:
        edge=pick_redge(web)
        if not edge:
           return web
        remove_edge(web,edge)
    return

# remove nodes in the web
def remove_nodes(web_info, node_list):
    for node in node_list:
        dpoint={}
        for n in  web_info[node]:
            for i in range(len(web_info[n])):
                if web_info[n][i]==node:
                    dpoint[n]=i; break
        for n in dpoint:
            del web_info[n][dpoint[n]]
        del web_info[node]
    return 

# get subweb information composed of input node list
def get_subweb (web, group, ex=None):
    extend=[]
    if ex:
        for other in web[ex]:
            if not other in group:
                extend.append(other)
    g=group+extend; subweb={}
    for node in g:
        subweb[node]=[]
        for other in web[node]:
            if other in g:
                subweb[node].append(other)
    return subweb

# find the nearest neighbor flow of one point.
def flow_pt (web, pt, strand):
    subweb=get_subweb(web,[pt], ex=pt)
    ptflow=get_flow(subweb, ref=[pt,strand])
    ptstrand=get_strand(subweb, ref=[pt,strand])
    return ptflow, ptstrand

# follow the flow in the graph to find one cluster.
def follow_flow (web, st):
    if len(web[st]) > 2:
        return None, None
    
    # forward direction
    pt=copy.deepcopy(st); strand=1; pos=0; group=[[pt,strand,pos]]
    while True:
        ptflow, ptstrand=flow_pt(web,pt,strand)
        if not ptflow:
            return group, [pt]
        if (len(ptflow[pt]['in']) + len(ptflow[pt]['out'])) > 2 and len(ptflow[pt]['in']) > 1:
            group=group[:-1]; break
        if len(ptflow[pt]['out'])==0 or len(ptflow[pt]['out']) > 1:
            break
        npt=ptflow[pt]['out'].keys()[0]; pos=ptflow[pt]['out'][npt]
        strand=ptstrand[npt]; pt=npt
        group.append([pt,strand,pos]);

    # reverse direction
    pt=copy.deepcopy(st); strand=1; rgroup=[]
    while True:
        ptflow, ptstrand=flow_pt(web,pt,strand)
        if (len(ptflow[pt]['in']) + len(ptflow[pt]['out'])) > 2 and len(ptflow[pt]['out']) > 1:
            rgroup=rgroup[:-1]; break
        if len(ptflow[pt]['in'])==0 or len(ptflow[pt]['in']) > 1:
            break
        npt=ptflow[pt]['in'].keys()[0]; pos=ptflow[pt]['in'][npt]
        strand=ptstrand[npt]; pt=npt
        rgroup.append([pt,strand,pos]);
    rgroup=rgroup[::-1]

    cluster=rgroup+group; node_list=[]
    for node_info in cluster:
        node_list.append(node_info[0])
    
    return cluster, node_list

# pick one node without branches in the graph
def pick_node (web):
    for node in web:
        if len(web[node]) < 3:
            return node
    return None

# find unitigs
def find_unis (olaps_info):
    web=get_web(olaps_info)
    rweb=get_rweb(web)
    unis_list=[]
    while rweb:
        st=pick_node(rweb)
        unis,node_list=follow_flow(rweb,st) 
        unis_list.append(unis); remove_nodes(rweb,node_list)
    return unis_list

# print out unitig list in the format
def printout (unis_list):
    output=''
    for i in range(len(unis_list)):
        output += 'UNI  %02d    %d    %d \n' % (i+1, len(unis_list[i]), 1)
        for j in range(len(unis_list[i])):
            if j ==0:
                num, strand, pos = unis_list[i][j]; mark='F'
                output += '  %03d   %s   %d \n' % (num, mark, 0)
                continue
            nnum,nstrand,npos=unis_list[i][j]
            if nstrand==1:
                nmark='F'
            else:
                nmark='R'
            if pos < 0:
                ref = pos
            else:
                ref = 0
            loc=npos-ref
            output += '  %03d   %s   %d \n' % (nnum, nmark, loc)
            num=nnum; strand=nstrand; pos=npos
    print output.strip()
    return

unis_list=find_unis(olaps_info)
printout(unis_list)
