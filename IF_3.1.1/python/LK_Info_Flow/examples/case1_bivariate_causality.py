# case 1 (bivariate causality):
import numpy as np
np.random.seed(0)  
from LK_Info_Flow import causal
#import matplotlib.pyplot as plt

#generate dataset
a11=0.3;a12=-0.4;a22=0.7;a21=0;b1=0.5;b2=0.5;
x=np.zeros([10000,])
y=np.zeros([10000,])
x[0]=0.4;y[0]=0.3;
for i in range(10000-1):
    x[i+1]=a11*x[i]+a21*y[i]+b1*np.random.normal();
    y[i+1]=a12*x[i]+a22*y[i]+b2*np.random.normal();

#print the structure of system
print ("x(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a11, a21,b1))
print ("y(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W" % (a12, a22,b2))


#initialization
X=np.ones((2,10000));
X[0]=x;X[1]=y;

#calculate the causality
causal_graph=causal.multi_causality_est_OLS(X[:,200:]);


#information flow
IF=np.squeeze(causal_graph.get('IF'));

#normalized information flow
nIF=np.squeeze(causal_graph.get('nIF'));

#significant test 
e99=np.squeeze(causal_graph.get('err_e99'));

## plot heatmap
# print(nIF[0,1])
# plt.matshow(nIF, 0,cmap='RdYlGn',vmin=-0.01,vmax=0.01)


# show the results
if abs(IF[0,1])>e99[0,1]:
    print('x -> y percent:%5.2f' % (nIF[0,1]*100)+'%')
else:
    print('x not -> y')


if abs(IF[1,0])>e99[1,0]:
    print('y -> x precent:%5.2f' % (nIF[1,0]*100)+'%')
else:
    print('y not -> x')
