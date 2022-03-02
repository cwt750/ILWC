%%
%读取参数
clear
clc
AAA=load('data\AAA.txt');
res_x=xlsread('data\res_x.xlsx');
res_y=xlsread('data\res_y.xlsx');
res_length=xlsread('data\res_length.xlsx');
%%
%设定输入参数
ILWC_length=51;
h0=1;
hjizhan=0.6;
angle2=45;
the_value=size(AAA,2);
f=3*10^14/1.06;
m=1.33+2.89e-6i;
attu=attu_K(f,m);
%%
%初始化数组
ILWC=zeros(ILWC_length,ILWC_length);
htop=zeros(the_value,length(res_x));
c1=zeros(1,length(res_x));
c2=zeros(1,length(res_x));
length_sum=zeros(1,length(res_x)+1);
h0_grid=zeros(1,length(res_x));
htop_grid=zeros(1,length(res_x));
CFLOS=zeros(the_value,length(res_x));
CFLOS_sum=zeros(1,the_value);
%%
for KKK=1:the_value
    %将ILWC数据从1*2500转为50*50
    for i=1:ILWC_length
        for j=1:ILWC_length
            ILWC(j,i)=AAA((i-1)*50+j,KKK);
        end
    end
    %求解对应网格的云顶高度
    syms asd
    for i=1:length(res_x)
        if(ILWC(res_x(i),res_y(i))<0.0217558)
            htop(KKK,i)=0;
        else
            c1(i)=4.27*exp(-4.93*(ILWC(res_x(i),res_y(i))+0.06))+54.12*exp(-61.25*(ILWC(res_x(i),res_y(i))+0.06))+1.71;
            c2(i)=3.17*c1(i)^(-3.04)+0.074;
            xxx=igamma(c1(i),(asd)/c2(i))/gamma(c1(i));
            htop(KKK,i)=vpasolve(xxx==0.06,[0 5])+h0;
        end
    end
    %判断每个网格激光链路穿越了多少距离的LWC以及有多少衰减
    for i=2:length(res_x)+1
        length_sum(i)=length_sum(i-1)+res_length(i-1);
    end
    for i=1:length(res_x)
        h0_grid(i)=hjizhan+length_sum(i)*tand(angle2);
        htop_grid(i)=hjizhan+length_sum(i+1)*tand(angle2);
    end
    for i=1:length(res_x)
        if(htop(KKK,i)==0)
            CFLOS(KKK,i)=0;
        end
        if(htop(KKK,i)~=0&&(htop_grid(i)<=h0||h0_grid(i)>=htop(KKK,i)))
            CFLOS(KKK,i)=0;
        end
        if(htop(KKK,i)~=0&&(h0_grid(i)>=h0&&htop_grid(i)<=htop(KKK,i)))
            CFLOS(KKK,i)=res_length(i)*tand(angle2)*attu;
        end
        if(htop(KKK,i)~=0&&((h0_grid(i)>=h0&&h0_grid(i)<=htop(KKK,i))&&htop_grid(i)>=htop(KKK,i)))
            CFLOS(KKK,i)=(htop(KKK,i)-h0_grid(i))/(htop_grid(i)-h0_grid(i))*res_length(i)*tand(angle2)*attu;
        end
        if(htop(KKK,i)~=0&&(h0_grid(i)<h0&&htop_grid(i)>=htop(KKK,i)))
            CFLOS(KKK,i)=(htop(KKK,i)-h0)/(htop_grid(i)-h0_grid(i))*res_length(i)*tand(angle2)*attu;
        end
        if(htop(KKK,i)~=0&&(h0_grid(i)<h0&&(htop_grid(i)>=h0&&htop_grid(i)<=htop(KKK,i))))
            CFLOS(KKK,i)=(htop_grid(i)-h0)/(htop_grid(i)-h0_grid(i))*res_length(i)*tand(angle2)*attu;
        end
    end
    CFLOS_sum(KKK)=sum(CFLOS(KKK,:));
end