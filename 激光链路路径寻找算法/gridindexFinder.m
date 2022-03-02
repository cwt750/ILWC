clear;clc;clf;
ILWC_length=51;
angle1=115;
d=5;
x = linspace(0,ILWC_length,ILWC_length+1);
y = linspace(0,ILWC_length,ILWC_length+1);
[X,Y] = meshgrid(x,y);
line(X,Y,'color','b');
line(X',Y','color','b');
axis equal;
axis([0 ILWC_length 0 ILWC_length]);
set(gca,'xtick',0:ILWC_length);

gridindex = reshape(1:ILWC_length^2,ILWC_length,ILWC_length)';
numposx = 0.5*(X(1:end-1,2:end)+X(1:end-1,1:end-1))-0.1;
numposy = 0.5*(Y(2:end,1:end-1)+Y(1:end-1,1:end-1));
for i = 1 : ILWC_length
    for j = 1 : ILWC_length
        text(numposx(i,j),numposy(i,j),num2str(gridindex(i,j)));
    end
end

P1 = [(ILWC_length)/2 (ILWC_length)/2];
P2 = [(ILWC_length)/2+sind(angle1)*d (ILWC_length)/2+cosd(angle1)*d];
segs = calLength(P1,P2);
line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');

data = cell2mat(struct2cell(segs));
res_x=zeros(1,size(data,3));
res_y=zeros(1,size(data,3));
res_length=zeros(1,size(data,3));
for i=1:length(res_x)
    res_x(i)=ILWC_length+1-data(3,1,i);
    res_y(i)=data(2,1,i);
    res_length(i)=data(1,1,i);
end
res_x( find(res_length<1e-10))=[];
res_y( find(res_length<1e-10))=[];
res_length( find(res_length<1e-10))=[];
xlswrite('data\res_x.xlsx',res_x);
xlswrite('data\res_y.xlsx',res_y);
xlswrite('data\res_length.xlsx',res_length);
% res_x=zeros(1,height(struct2table(segs)));
% res_y=zeros(1,height(struct2table(segs)));
% for i=1:length(res_x)
%     res_x(i)=segs(i).index_x;
%     res_y(i)=11-segs(i).index_y;
% end

% display('所经过的网格序号\长度分别为:');
% for i = 1 : size(segs,2)
%     display(['序号: ' num2str(segs(i).index)]);
%     display(['长度: ' num2str(segs(i).length)]);
% end
