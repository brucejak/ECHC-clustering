function [Idx,Idx_with_noise,cutnum,cutlink_ori,p,firstlayer_loc_onsortp,mass,R,cutlinkpower_all,layer_border_ALL,order,maxind,C2C_all,eachlayerlabel,m_and_d_all,labels2] = TorqueClustering(ALL_DM,NAB)
% input: ALL_DM: distance matrix of n x n, where n means the number of total data samples;
%       NAB(number of abnormal connections): 
%                 when NAB=0, exploit the automatic mode to determine the number of clusters;
%                 when NAB��0��remove the top NAB connections with the largest torque to obtain the clusters.                       
% output: Idx: final clusters labels;
%         cutnum: number of abnormal connections that are removed when exploiting the automatic mode.
datanum=size(ALL_DM,1);cutlinkpower_all=[];C2C_all=[];m_and_d_all=[];eachlayerlabel=[];
ljmat=sparse(datanum,datanum);
%community=cell(1,datanum);
dataloc=1:1:datanum;
%for i=1:1:datanum
%community{1,i}=i;
%end
community=arrayfun(@(x) {dataloc(x)},1:length(dataloc));
commu_DM=ALL_DM;
%commu_DM1=commu_DM;
community_num=datanum;
% i=1:1:community_num
%commu_DM1(i,i)=NaN;
%end

G=sparse(community_num,community_num);
%for i=1:1:community_num
%gd=find(commu_DM1(i,:)==min(commu_DM1(i,:)));
   
        %G(i,gd(1))=1;
        %G(gd(1),i)=1;
%end
[~,Order]=sort(commu_DM,2);
neiborloc=cell(1,community_num);
for i=1:1:community_num
    G(i,Order(i,2))=1;
    G(Order(i,2),i)=1;
    neiborloc{1,i}=Order(i,2);
end
SG=graph(G);
BINS = conncomp(SG);

%linklocCell=[];
%[a,b,~]=find(G);
%for i=1:1:community_num
%neiborloc{1,i}=b(a==i)';
%end
cur_NC=numel(unique(BINS));
disp(['The number of clusters in this layer is: ',num2str(cur_NC)])
[ljmat,cutlinkpower,C2C,m_and_d]=Updateljmat(ljmat,neiborloc,community,commu_DM,G,ALL_DM);
[cutlinkpower,ljmat,C2C]=uniqueZ(cutlinkpower,ljmat,C2C);
firstlayer_conn_num=size(cutlinkpower,1);
cutlinkpower_all=[cutlinkpower_all;cutlinkpower];
C2C_all=[C2C_all;C2C];
m_and_d_all=[m_and_d_all;m_and_d];
layer_border_ALL=[firstlayer_conn_num];
%th=1;

ljmat_layerG = graph(ljmat);
BINS_layer = conncomp(ljmat_layerG);
eachlayerlabel=[eachlayerlabel BINS_layer'];

old_G=G;
while 1
    
if numel(unique(BINS))==1 %% 2021.9.23����
break
end
   
Idx=BINS;uni_Idx=unique(Idx);
num_uni_Idx=length(uni_Idx);

community_new=cell(1,num_uni_Idx);
for i=1:1:num_uni_Idx
uniloc=(uni_Idx(i)==Idx);
community_new{1,i}=[community{uniloc}];
end
community=community_new;
community_num=size(community,2);
%commu_DM=[];
commu_DM=zeros(community_num,community_num);
%linklocCell=cell(community_num,community_num);
for i=1:1:community_num
for j=1:1:community_num
commu_DM(i,j)=ps2psdist(community{1,i},community{1,j},ALL_DM);
end
end
%commu_DM1=commu_DM;
for i=1:1:community_num
commu_DM(i,i)=-inf;
end

G=sparse(community_num,community_num);
%Order=[];
[~,Order]=sort(commu_DM,2);
neiborloc=cell(1,community_num);
for i=1:1:community_num
%gd=find(commu_DM1(i,:)==min(commu_DM1(i,:)));
    if numel(community{i})<=numel(community{Order(i,2)})
        G(i,Order(i,2))=1;
        G(Order(i,2),i)=1;
        neiborloc{1,i}=Order(i,2);
    end  
end
if isequal(G,old_G)
    break;
end
    
SG=graph(G);
BINS = conncomp(SG);
%eachlayerlabel=[eachlayerlabel BINS'];
old_G=G;
%[a,b,~]=find(G);
%for i=1:1:community_num
%neiborloc{1,i}=b(a==i)';
%end
%numel(unique(BINS))
cur_NC=numel(unique(BINS));
disp(['The number of clusters in this layer is: ',num2str(cur_NC)])

[ljmat,cutlinkpower,C2C,m_and_d]=Updateljmat(ljmat,neiborloc,community,commu_DM,G,ALL_DM);
[cutlinkpower,ljmat,C2C]=uniqueZ(cutlinkpower,ljmat,C2C);
cutlinkpower_all=[cutlinkpower_all;cutlinkpower];
C2C_all=[C2C_all;C2C];
m_and_d_all=[m_and_d_all;m_and_d];
layer_border=size(cutlinkpower_all,1);
layer_border_ALL=[layer_border_ALL,layer_border];

ljmat_layerG = graph(ljmat);
BINS_layer = conncomp(ljmat_layerG);
eachlayerlabel=[eachlayerlabel BINS_layer'];

if numel(unique(BINS))==1
break
end
end
mass=cutlinkpower_all(:,5).*cutlinkpower_all(:,6);
R=cutlinkpower_all(:,7).^2;

%norm_mass=mapminmax(mass',0.001,0.999);
%norm_R=mapminmax(R',0.001,0.999);

p=mass.*R;

R_mass=R./mass;



[~,order]=sort(p,'descend');%loc=1:numel(p);
[~,order_2]=sort(order);%%
firstlayer_loc_onsortp=order_2(1:firstlayer_conn_num);

if NAB==0
[cutnum, maxind] =Nab_dec(p,mass,R,firstlayer_loc_onsortp);
   % cutlink=cutlinkpower_all(order(1:cutnum),:);
   %  cutlink_ori=cutlink;
   % cutlink(:,[1 2 5 6 7])=[];
   noise_loc=intersect(intersect((find(R>mean(R))), (find(mass<mean(mass)))), (find(R_mass>mean(R_mass))));
   %noise_loc=intersect((find(R>mean(R))), (find(mass<mean(mass))));
   cutlink1=cutlinkpower_all(order(1:cutnum),:);
   cutlink2=cutlinkpower_all((union(order(1:cutnum),noise_loc)),:);
   cutlink_ori=cutlink1;
   cutlink1(:,[1 2 5 6 7])=[];
   cutlink2(:,[1 2 5 6 7])=[];
%end
else
    cutnum = NAB;
  %  cutlink=cutlinkpower_all(order(1:cutnum),:);
  %  cutlink_ori=cutlink;
  %  cutlink(:,[1 2 5 6 7])=[];
  maxind=0;
  noise_loc=intersect(intersect((find(R>mean(R))), (find(mass<mean(mass)))), (find(R_mass>mean(R_mass))));
   cutlink1=cutlinkpower_all(order(1:cutnum),:);
   cutlink2=cutlinkpower_all((union(order(1:cutnum),noise_loc)),:);
   cutlink_ori=cutlink1;
   cutlink1(:,[1 2 5 6 7])=[];
   cutlink2(:,[1 2 5 6 7])=[];
end
ljmat1=ljmat;
%cutlinknum=size(cutlink,1);
%for i=1:1:cutlinknum
%ljmat(cutlink(i,1),cutlink(i,2))=0;
%ljmat(cutlink(i,2),cutlink(i,1))=0;
%end
%ljmat_G=graph(ljmat);
%BINS = conncomp(ljmat_G);
%Idx=BINS';
%NC=numel(unique(BINS));
cutlinknum1=size(cutlink1,1);
for i=1:1:cutlinknum1
ljmat(cutlink1(i,1),cutlink1(i,2))=0;
ljmat(cutlink1(i,2),cutlink1(i,1))=0;
end
ljmat_G=graph(ljmat);
BINS = conncomp(ljmat_G);
labels1=BINS';

cutlinknum2=size(cutlink2,1);
for i=1:1:cutlinknum2
ljmat1(cutlink2(i,1),cutlink2(i,2))=0;
ljmat1(cutlink2(i,2),cutlink2(i,1))=0;
end
ljmat1_G=graph(ljmat1);
BINS = conncomp(ljmat1_G);
labels2=BINS';

Idx=labels1;
%Idx_with_noise=Final_label(labels1,labels2);
 Idx_with_noise=[];
end
