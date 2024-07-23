import graphviz as gz
#https://www.graphviz.org/Packages/stable/windows/10/cmake/Release/x64/
import os
import numpy as np

def addEdge(a, b):
    global edgeLinks
    if a not in edgeLinks: edgeLinks[a] = set()
    if b not in edgeLinks: edgeLinks[b] = set()
    edgeLinks[a].add(b)
    edgeLinks[b].add(a)

def plot_causal_graph(A,max_lag,f_name='model'):
    global edgeLinks
    edgeLinks = dict()
    edge=[]
    dot = gz.Digraph()
    A=np.array(A)+0;
    N=A.shape;
    tmp=A[:,0:N[0]]+0;
    A[:,0:N[0]]=A[:,N[0]*(max_lag-1):]+0;
    A[:,N[0]*(max_lag-1):]=tmp+0;
    A=A.T; nodes=[];
    for ix in range(N[0]*max_lag):
        for iy in range(N[0]):
            if ix!=iy:
                if np.isnan(A[ix,iy]):
                    tmp=0
                else:
                        addEdge(str(ix), str(iy+N[0]*2))
                        dot.edge(str(ix), str(iy+N[0]*2),str(np.ceil(A[ix,iy]*10000)/100))
                        nodes.append([ix,iy+N[0]*2])



    for i in np.unique(nodes):
        k=int(np.floor(i/N[0]))+1;
        if k<max_lag:
            dot.node(str(i), str(i-(k-1)*N[0])+'(-'+str(k)+')')
        else:
            dot.node(str(i), str(i-(k-1)*N[0]))

    #dot.edges(edge)
    dot.render(f_name, view=False)
    os.system('dot -Tpng '+f_name+' -o '+f_name+'.png')# system need install graphiz 

def causal_graph(causal_matrix,significance=None,c_threshold=0.001,s_threshold=0.99,name='model'):
    a=abs(causal_matrix)
    if significance is None:
        1
    else:
        b=abs(significance)
        a[np.where(b>s_threshold)]=np.NaN
    
    a[np.where(a<c_threshold)]=np.NaN
    N=a.shape;
    c=np.reshape(a[:,:,::1],[N[0],N[0]*N[2]]).T
    for i in range(N[0]):
        c[i,i]=np.NaN

    data=c.transpose()
    plot_causal_graph(data,N[2],name)