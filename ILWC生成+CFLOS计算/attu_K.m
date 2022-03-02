function res=attu_K(f,m)
lambda=3*10^14/f;
qext=zeros(1,10000);
num=zeros(1,10000);
kext=zeros(1,10000);
for r=0.01:0.01:100
    a=2*pi*r/lambda;
    index=floor((r-0.01)*100+1);
    [~,~,qext_temp,~,~,~]=mie(a,1.5+0.02*1i,pi/5);
    qext(index)=qext_temp;
    num(index)=2.373*r^6*exp(-1.5*r);
    kext(index)=qext(index)*pi*r^2;
end
sum=0;
for i=1:10000
    sum=sum+kext(i)*num(i)*0.01;
end
res=sum*4343*1e-6;





