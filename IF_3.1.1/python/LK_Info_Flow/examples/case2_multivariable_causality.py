# case 2 (multivariable causality):
import numpy as np
np.random.seed(0)  
from LK_Info_Flow import causal

# generate dataset
a11=0.3;a21=0;  a31=0;  a41=0;  b1=0.4;
a12=0.5;a22=0.7;a32=0.1;a42=0;  b2=0.5;
a13=0;  a23=0.3;a33=0.5;a43=0;  b3=0.6;
a14=0.2;a24=0.4;a34=0.3;a44=0.1;  b4=0.3;
xx=np.zeros([1001,4])
xx[0]=np.array([0.4,0.5,0.6,0.7])
for i in range(1000):
    xx[i+1,0]=a11*xx[i,0]+a21*xx[i,1]+a31*xx[i,2]+a41*xx[i,3]+b1*np.random.normal();
    xx[i+1,1]=a12*xx[i,0]+a22*xx[i,1]+a32*xx[i,2]+a42*xx[i,3]+b2*np.random.normal();
    xx[i+1,2]=a13*xx[i,0]+a23*xx[i,1]+a33*xx[i,2]+a43*xx[i,3]+b3*np.random.normal();
    xx[i+1,3]=a14*xx[i,0]+a24*xx[i,1]+a34*xx[i,2]+a44*xx[i,3]+b4*np.random.normal();

#print the structure of system
print ("x1(i+1)=%.2f * x1(i) + %.2f * x2(i)+%.2f * x3(i) + %.2f * x4(i) + %.2f W" % (a11, a21, a31,a41,b1))
print ("x2(i+1)=%.2f * x1(i) + %.2f * x2(i)+%.2f * x3(i) + %.2f * x4(i) + %.2f W" % (a12, a22, a32,a42,b2))
print ("x3(i+1)=%.2f * x1(i) + %.2f * x2(i)+%.2f * x3(i) + %.2f * x4(i) + %.2f W" % (a13, a23, a33,a43,b3))
print ("x4(i+1)=%.2f * x1(i) + %.2f * x2(i)+%.2f * x3(i) + %.2f * x4(i) + %.2f W" % (a14, a24, a34,a44,b4))

#initialization
Nxx=np.array(np.shape(xx));
IF=np.zeros([Nxx[1],Nxx[1]]);
Nxx=np.array(np.shape(xx));
t=np.linspace(0,1000,1001) #time series

#calculate the causality
cau2=causal.multi_causality_est_OLS(X=xx.T,series_temporal_order=t)


#information flow
IF = np.squeeze(cau2.get('IF'))
#normalized information flow
nIF= np.squeeze(cau2.get('nIF'));
#significant test 
err= np.squeeze(cau2.get('err_e99'));



#qualitative results 
for i in range(4):
    for j in range(4):
        if abs(IF[i,j])>err[i,j]:
            IF[i,j]=1;

print('xj -> xi:')
f='    j';
Nxx=np.array(np.shape(xx));
for i in range(Nxx[1]):
    f='    '+f;
print(f)
f='    %6d'%1;
for i in range(Nxx[1]-1):
    f=f+'%5d'%(i+2);
print(f)
for j in range(Nxx[1]):
    if j== np.floor(Nxx[1]/2):
        f=' i';
    else:
        f='  ';
    f=f+'%5d'%(j+1)
    for i in range(Nxx[1]):
        f=f+'  %d  ' %(IF[j,i])
    print(f)


#quantitative causality
print('Tj->i:')
f='    j';
Nxx=np.array(np.shape(xx));
for i in range(Nxx[1]):
    f='     '+f;
print(f)
f='        %7d'%1;
for i in range(Nxx[1]-1):
    f=f+'%9d'%(i+2);
print(f)
for j in range(Nxx[1]):
    if j== np.floor(Nxx[1]/2):
        f=' i';
    else:
        f='  ';
    f=f+'%5d'%(j+1)
    for i in range(Nxx[1]):
        f=f+' %8.4f ' %(nIF[j,i])
    print(f)


#significant test (err99)
print('e99:')
f='    j';
Nxx=np.array(np.shape(xx));
for i in range(Nxx[1]):
    f='     '+f;
print(f)
f='        %7d'%1;
for i in range(Nxx[1]-1):
    f=f+'%9d'%(i+2);
print(f)
for j in range(Nxx[1]):
    if j== np.floor(Nxx[1]/2):
        f=' i';
    else:
        f='  ';
    f=f+'%5d'%(j+1)
    for i in range(Nxx[1]):
        f=f+' %8.4f ' %(err[j,i])
    print(f)


#significant test (p-value)
print('p-value:')
p=np.squeeze(cau2.get('p'));
f='    j';
Nxx=np.array(np.shape(xx));
for i in range(Nxx[1]):
    f='     '+f;
print(f)
f='        %7d'%1;
for i in range(Nxx[1]-1):
    f=f+'%9d'%(i+2);
print(f)
for j in range(Nxx[1]):
    if j== np.floor(Nxx[1]/2):
        f=' i';
    else:
        f='  ';
    f=f+'%5d'%(j+1)
    for i in range(Nxx[1]):
        f=f+' %8.4f ' %(p[j,i])
    print(f)

    
    
