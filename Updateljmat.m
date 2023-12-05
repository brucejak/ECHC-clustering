function [new_ljmat,cutlinkpower,C2C,m_and_d] = Updateljmat(old_ljmat,neiborloc,community,commu_DM,G,ALL_DM)
%NNC_UPDATELJMAT Summary of this function goes here
%   Detailed explanation goes here
% On the basis of GC_updateljmat_2, consider the cluster whose mass is one.
% ���£���ԭ�л����ϣ��޳���neiborloc�ж���Ŀ����ԣ�ֻ�ж��У����ߣ��ޡ�
%�汾6���ú���ȥ�����������������������
pd=numel(community{1});
if pd>1
community_num=numel(community);

    %cutlinknum=0;
%for i=1:1:community_num
    %cutlinknum=cutlinknum+numel(neiborloc{i});
%end
cutlinknum = sum(~cellfun('isempty', neiborloc));
cutlinkpower=zeros(cutlinknum,7);
C2C=cell(cutlinknum,2);
m_and_d=cell(community_num,3);
th=1;
    for i=1:1:community_num
        if isempty(neiborloc{i})
            m_and_d{i,1}=community{i};
            m_and_d{i,2}=numel(community{i});
        end
        if ~isempty(neiborloc{i})
            [linkloc1,linkloc2]=mindisttwinsloc(community{i},community{neiborloc{i}},ALL_DM);
            xx=min(community{i});
            yy=min(community{neiborloc{i}});
            C2C{th,1}=community{i};
            C2C{th,2}=community{neiborloc{i}};
            m_and_d{i,1}=community{i};
            m_and_d{i,2}=numel(community{i});
            m_and_d{i,3}=commu_DM(i,neiborloc{i});
            %linkloc1=xx(1);
           % linkloc2=yy(1);
             %linkloc1 = community{i}(1);
            %linkloc2 = community{neiborloc{i}}(1);
            %linkloc1=linklocCell{i,neiborloc{i}}(1);
            %linkloc2=linklocCell{i,neiborloc{i}}(2);
            old_ljmat(linkloc1,linkloc2)=1;
            old_ljmat(linkloc2,linkloc1)=1;
            %[th] = row0loc(cutlinkpower);
            cutlinkpower(th,1)=xx;
            cutlinkpower(th,2)=yy;
            cutlinkpower(th,3)=linkloc1;
            cutlinkpower(th,4)=linkloc2;
            cutlinkpower(th,5)=length(community{i});
            cutlinkpower(th,6)=length(community{neiborloc{i}});
            cutlinkpower(th,7)=commu_DM(i,neiborloc{i});
            th=th+1;
        end
    end
end

if pd==1
    community_num=numel(community);
cutlinkpower=zeros(community_num,7);%%�޸�Ϊ7
C2C=cell(community_num,2);
m_and_d=cell(community_num,3);
th=1;
    for i=1:1:community_num
            %[linkloc1,linkloc2]=maxdisttwinsloc(community{i},community{neiborloc{i}},ALL_DM);
             linkloc1 = community{i}(1);
            linkloc2 = community{neiborloc{i}}(1);
            C2C{th,1}=community{i};
            C2C{th,2}=community{neiborloc{i}};
            m_and_d{i,1}=community{i};
            m_and_d{i,2}=numel(community{i});
            m_and_d{i,3}=commu_DM(i,neiborloc{i});
            %linkloc1=linklocCell{i,neiborloc{i}}(1);
            %linkloc2=linklocCell{i,neiborloc{i}}(2);
            %old_ljmat(linkloc1,linkloc2)=1;
            %old_ljmat(linkloc2,linkloc1)=1;
            %[th] = row0loc(cutlinkpower);
            cutlinkpower(th,1)=linkloc1;
            cutlinkpower(th,2)=linkloc2;
            cutlinkpower(th,3)=linkloc1;
            cutlinkpower(th,4)=linkloc2;
            cutlinkpower(th,5)=1;
            cutlinkpower(th,6)=1;
            cutlinkpower(th,7)=commu_DM(i,neiborloc{i});
            th=th+1;
    end
        old_ljmat=G;
end

new_ljmat=old_ljmat;

    
    
end

