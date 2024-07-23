# case 8(temperal varying):
import numpy as np
from tqdm import tqdm
from LK_Info_Flow import multi_causality_est
a11=0.3;a12=-0.4;a22=0.7;a21=0;b1=0.4;b2=0.5;
x=np.zeros([100000,])
y=np.zeros([100000,])
x[0]=0.4;y[0]=0.3;
for i in range(50000-1):
    x[i+1]=a11*x[i]+a21*y[i]+b1*np.random.uniform(-1,1);
    y[i+1]=a12*x[i]+a22*y[i]+b2*np.random.uniform(-1,1);

print ("x(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a11, a21,b1))
print ("y(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a12, a22,b2))

a11=0.3;a12=0.0;a22=0.7;a21=0.3;b1=0.4;b2=0.5;
for i in range(50000-1,100000-1):
    x[i+1]=a11*x[i]+a21*y[i]+b1*np.random.uniform(-1,1);
    y[i+1]=a12*x[i]+a22*y[i]+b2*np.random.uniform(-1,1);

print ("x(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a11, a21,b1))
print ("y(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a12, a22,b2))




print('caculating:')
window_size=1000;
T=np.zeros([100000,]);E99=np.zeros([100000,]);T1=np.zeros([100000,]);E991=np.zeros([100000,])
for i in tqdm(range(10000,90000)):
    tmp=np.zeros([2,window_size*2]);
    tmp[0]=x[i-window_size:i+window_size];
    tmp[1]=y[i-window_size:i+window_size];
    cau=multi_causality_est(tmp);
    T21=cau['IF'].squeeze();
    e99=cau['err_e90'].squeeze();
    T[i]=T21[0,1];E99[i]=e99[0,1];
    T1[i]=T21[1,0];E991[i]=e99[1,0];




from matplotlib import pyplot as plt
# import matplotlib.pyplot as plt
plt.plot(np.arange(20000,80001),  T[20000:80001],color='r')
plt.plot(np.arange(20000,80001),  np.zeros([60001,]),color='k')
plt.fill_between(np.arange(20000,80001),
                 T[20000:80001]-E99[20000:80001],
                 T[20000:80001]+E99[20000:80001],
                 color='b',
                 alpha=0.2)
plt.xlim([20000,80000])
plt.title('information flow from $X_2$ to $X_1$ with 90 precent significant test')
plt.show()