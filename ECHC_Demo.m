%%%%%COIL-100 dataset
%if 0
clear
disp('COIL-100 dataset')
data=importdata('Coil100.mat');
datalabels=data.datalabels;data=data.data;
K=numel(unique(datalabels));
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3 %%%%%repeat 3 times
BINS = cls(:,5)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

clear; %%%%%COIL-40 dataset
disp('COIL-40 dataset')
data=importdata('Coil40.mat');
datalabels=data.gt;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,4)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

clear %%%%%COIL-20 dataset
disp('COIL-20 dataset')
data=importdata('COIL20.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,3)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

clear %%%% FRGCv2 dataset
disp('FRGCv2 dataset')
data=importdata('FRGCv2.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,1)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

clear%%%%%%%%%UMIST dataset
disp('UMIST dataset')
data=importdata('UMIST.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,4)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

clear%%%%%%%%%EYaleB dataset
disp('EYaleB dataset')
data=importdata('EYaleB.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,1)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)

%end

clear %%%%%miceprotein dataset
disp('miceprotein dataset')
data=importdata('miceprotein.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,2)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)


clear %%%%%%%segment dataset
disp('segment dataset')
data=importdata('segment.mat');
datalabels=data.datalabels;data=data.data;
DM=pdist2(data,data,'cosine');
cls=CHC_fg(DM,0);
stddata=std(data);loc=find(stddata==0);data(:,loc)=[];
for j=1:1:3
BINS = cls(:,3)'; K=numel(unique(datalabels));
t=full(ind2vec(BINS));
commandhistory
x=data';
net=patternnet(100);
net=train(net,x,t);
commandhistory
mapdata=zeros(size(data,1),100);
for i=1:1:size(data,1)
mapdata(i,:)=(net.IW{1}*x(:,i)+net.b{1})';
end
mapDM=pdist2(mapdata,mapdata,'cosine');NM=zeros(1,50);AC=zeros(1,50);
parfor i=1:50
W=gac_createdigraph(mapDM, i);
W=full(W);WD = 1./W;
cl1=CHC(WD,K-1);NM(i)=nmi(cl1,datalabels);AC(i)=accuracy(cl1,datalabels+1)/100;
end
max_NM(j)=max(NM),max_AC(j)=max(AC)
end
average_NM = mean(max_NM),average_AC = mean(max_AC)