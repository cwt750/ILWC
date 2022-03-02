clear
clc
the_value=5;
ILWC_length=51;
AAA=zeros(ILWC_length^2,the_value);
fin1=zeros(ILWC_length^2,the_value);
fin2=zeros(ILWC_length^2,the_value);
G=zeros(ILWC_length^2,the_value);
%设定输入参数
b1=0.000717;b2=0.0000201;
r1=0.349;r2=0.830;
m=xlsread('data\m.xlsx');o=xlsread('data\o.xlsx');PCLW=xlsread('data\PCLW.xlsx');
t=600;
%画图+画图坐标轴
a=(0:1:ILWC_length);
b=(0:1:ILWC_length);
[x,y]=meshgrid(a,b);
name1={'0','10','20','30','40'};
name2=fliplr(name1);
sitex=[0:1:ILWC_length];
%导入动态特性参数B1,B2
B01=-b1*ones(1,ILWC_length^2);
B1=diag(B01);
B02=-b2*ones(1,ILWC_length^2);
B2=diag(B02);      
%求解S
A=zeros(ILWC_length^2,2);
d=zeros(ILWC_length^2,ILWC_length^2);
C=zeros(ILWC_length^2,ILWC_length^2);
F1=zeros(ILWC_length^2,ILWC_length^2);
F2=zeros(ILWC_length^2,ILWC_length^2);
S1=zeros(ILWC_length^2,ILWC_length^2);
S2=zeros(ILWC_length^2,ILWC_length^2);
for i=1:ILWC_length^2
    A(i,:)=[((i-mod(i,ILWC_length))/ILWC_length),(mod(i,ILWC_length)-1)];
    if mod(i,ILWC_length)==0
        A(i,:)=[(i/ILWC_length-1),ILWC_length-1];
    end
end  
for i=1:ILWC_length^2
    for j=1:ILWC_length^2
        d(i,j)=sqrt(sum((A(i,:)-A(j,:)).^2));
    end
end
C=0.35*exp(-d/7.8)+0.65*exp(-d/225.3);
F1=2*b1*C; 
F2=2*b2*C;
S1=chol(F1,'lower');
S2=chol(F2,'lower');
%导入X0
X01=zeros(ILWC_length^2,the_value+1);
X02=zeros(ILWC_length^2,the_value+1);
T1=S1*S1';
T2=S2*S2';
sigma1=T1/(2*b1);
sigma2=T2/(2*b2);
% mu1=(qfuncinv(PCLW/2)/(r1+r2))*ones(2500,1);
% mu2=(qfuncinv(PCLW/2)/(r1+r2))*ones(2500,1);
mu1=zeros(ILWC_length^2,1);
mu2=zeros(ILWC_length^2,1);
X01(:,1)=mvnrnd(mu1,sigma1);
X02(:,1)=mvnrnd(mu2,sigma2);
for KKK=1:the_value
%求解布朗运动
N=10;
w=zeros(ILWC_length^2,N+1);
r=t/N;
w(:,1) = 0;
for i=1:ILWC_length^2
w(i,2)=sqrt(r)*randn;
for j=3:N+1
    w(i,j)=w(i,j-1)+sqrt(r)*randn; % √r*N(0,1)
end
end
ww=zeros(ILWC_length^2,N);
for i=1:N
    ww(:,i)=(w(:,i+1)-w(:,i))/r;
end
%求解X
X1=zeros(ILWC_length^2,1);
X2=zeros(ILWC_length^2,1);
sum1=zeros(ILWC_length^2,N);
sum2=zeros(ILWC_length^2,N);
q1=zeros(ILWC_length^2,N);
q2=zeros(ILWC_length^2,N);
for i=1:N
    t1=(i-1)*r;
    t2=i*r;
    for e=1:ILWC_length^2
        for f=1:e
            sum1(e,i)=sum1(e,i)+S1(e,f)*ww(f,i);
            sum2(e,i)=sum2(e,i)+S2(e,f)*ww(f,i);
        end
    end
    q1(:,i)=sum1(:,i)/b1*(exp(b1*t2)-exp(b1*t1));
    q2(:,i)=sum2(:,i)/b2*(exp(b2*t2)-exp(b2*t1));
end    
for i=1:N
    fin1(:,KKK)=fin1(:,KKK)+q1(:,i);
    fin2(:,KKK)=fin2(:,KKK)+q2(:,i);
end
X1=exp(-b1*t)*X01(:,KKK)+exp(-b1*t)*fin1(:,KKK);
X2=exp(-b2*t)*X02(:,KKK)+exp(-b2*t)*fin2(:,KKK);
X01(:,KKK+1)=X1;
X02(:,KKK+1)=X2;
%求解G
G(:,KKK)=r1*X1+r2*X2;    
%求解ILWC
L=zeros(ILWC_length,ILWC_length);
ILWC=zeros(ILWC_length,ILWC_length);
for i=1:ILWC_length
    for j=1:ILWC_length
        L(i,j)=G(ILWC_length*(j-1)+i,KKK);
        v=qfunc(L(i,j));
        if(L(i,j)<=qfuncinv(PCLW(i,j)))
            ILWC(i,j)=0;
        else
            ILWC(i,j)=exp(m(i,j)+o(i,j)*qfuncinv(v/PCLW(i,j)));
        end
    end
end
aaa=ILWC(:);
for lll=1:length(aaa)
    AAA(lll,KKK)=aaa(lll);
end
if(mod(KKK,100)==0)
    figure(KKK)
    geoshow(x,y,ILWC,'Displaytype','texturemap');
    hcb=colorbar('eastoutside');
    [cmin,cmax]=caxis;
    caxis([cmin,cmax]);
    set(get(hcb,'Xlabel'),'String','ILWC/mm');
end
end
dlmwrite('AAA.txt',AAA,'delimiter',' ');
movefile('AAA.txt','data');