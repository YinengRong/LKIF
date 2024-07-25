# case 6 (Large scale computing cost):
import time
from tqdm import tqdm
import numpy as np
from LK_Info_Flow import multi_causality_est

#generate data
NN=[200,1000,5000,10000]  #time steps
MM=[10,50,100]  #var numbers
nlag=0; #time delay
for iN in range(3):
    N=NN[iN];#time steps
    for iM in range(3):
        M=MM[iM];# varilabe number
        jj=0;
        tcost=[];IF=[];A_gt=[]
        print([iN,iM])
        for icase in tqdm(range(100)):
        
            b=(np.random.randint(4)+1)/20;

            max_x=10 #if x >10 then regard the model is unstable
            while max_x>=10:
                A=np.random.normal(size=(M,M));
                A[np.where(A<1)]=0;
                A[np.where(A>=1)]=0.2;
                A=A*(1-np.eye(M,M))+0.2*np.eye(M,M);
                tmp=np.ones((len(A[np.where(A==0.2)]),))
                tmp[::2]=-tmp[::2]
                A[np.where(A==0.2)]=A[np.where(A==0.2)]*tmp

                xx=np.zeros([N+1,M])
                xx[0,:]=np.random.normal(size=(M,))/2;
                for i in range(N):
                    xx[i+1]=np.dot(A,xx[i])+b*np.random.normal(size=(M,));
                #max_x=np.max(abs(xx));
                if np.isnan(np.max(abs(xx))):
                    print(jj)
                    max_x=10
                    jj=jj+1;
                else:
                    max_x=np.max(abs(xx))
#initialization
            xx1=np.zeros([N-nlag,M+M*nlag])
            xx1[:,:M]=xx[:N-nlag,:];
            for i in range(nlag):
                xx1[:,(1+i)*M:(2+i)*M]=xx[i+1:N-nlag+i+1,:]
#calculate information flow            
            ts=time.time()
            cau=multi_causality_est(X=xx1.T);
            #T21=cau;
            te=time.time()
            tcost.append(te-ts)
            IF.append(cau)
            A_gt.append(A)
#save results
        np.savez('result_%04d_%04d.npz'%(M,N), **{'A':A_gt, 'IF':IF, 'tcost':tcost})


#load results
NN=[200,1000,5000]
MM=[10,50,100]
nlag=0; #time delay
HD=[];Tcost=[];
for iN in range(3):
    N=NN[iN];#time steps
    for iM in range(3):
        M=MM[iM];# varilabe number
        data = np.load('result_%04d_%04d.npz'%(M,N),allow_pickle=True)
        A0=data['A']
        IF=data['IF']
        tcost=data['tcost']

        hd=[];
        for icase in range(100):
            a=IF[icase]['nIF'].squeeze();
            b=IF[icase]['p'].squeeze();
            A=A0[icase]+0;
            #A=A_gt[icase]+0;
            A[np.where(abs(A)>0.1)]=1
            A[np.where(abs(A)<0.1)]=0
            a[np.where(b<(1-1/(N-M)))]=np.NaN
            M=a.shape[0]

            a[np.where(abs(a)<(1/M**2)/100)]=np.NaN
            a[~np.isnan(a)]=1
            a[np.isnan(a)]=0
            hd.append(np.sum(abs(a-A)))
        HD.append(np.array(hd)/M**2)
        Tcost.append(tcost)


#plot costs
import matplotlib.pyplot as plt
fig = plt.figure()
xticklabels=[' dim= 10           \n sample= 200               ',\
             '50\n200','100\n200','10\n1000','50\n1000','100\n1000','10\n5000','50\n5000','100\n5000']
plt.figure(figsize=(12, 3))
bp = plt.boxplot(np.array(HD).T,labels=xticklabels)
axes = plt.gca()
left, right = axes.get_xlim()
axes.hlines(y=[0.01,0.1], xmin=left, xmax=right, linestyles='dashed')
plt.ylim(0,0.2)
plt.ylabel('a) normalized structural\nHamming distance')
plt.show()


fig = plt.figure()
xticklabels=[' dim= 10           \n sample= 200               ',\
             '50\n200','100\n200','10\n1000','50\n1000','100\n1000','10\n5000','50\n5000','100\n5000']
plt.figure(figsize=(12, 3))
bp = plt.boxplot(np.array(Tcost).T,labels=xticklabels)
# axes = plt.gca()
# left, right = axes.get_xlim()
# axes.hlines(y=[0.01,0.1], xmin=left, xmax=right, linestyles='dashed')
plt.ylim(0,0.15)
plt.ylabel('b) Runtime (s)')
plt.show()