%======================����agent��������============================
clc
clear
close all
%===============С�����������=============
celllength=20;
celllength1=celllength;
celllength2=celllength;
pneighborlink=1;%���ھӵ����Ӹ���
pcut=0.5;%�ϵ���������
link=zeros(celllength1,celllength2,celllength1,celllength2);
%===============����������С��������================
%==============����������֮�����ϵ====================
for i=1:celllength1
    for j=1:celllength2
        for m=1:celllength1
            for n=1:celllength2
                link(i,j,i,j)=0;
                if link(m,n,i,j)==0
                    if rand<=pneighborlink%��һ���ĸ��ʽ������ھӵ�����
                        if (abs(i-m)==1&&j==n)||(abs(j-n)==1&&i==m)
                            link(m,n,i,j)=1;
                            link(i,j,m,n)=1;
                        end
                    end
                end
            end
        end
    end
end
%===========�ϵ�����===========================
for i=1:celllength1
    for j=1:celllength2
            if rand<pcut
                %==========�Ͽ�����=============
                if i>1&&i<celllength1&&j>1&&j<celllength2
                    randne=randsrc(1,1,[1:4]);%�ĸ��ھ�������һ��
                    if randne==1
                        link(i-1,j,i,j)=0;
                        link(i,j,i-1,j)=0;
                    elseif randne==2
                        link(i+1,j,i,j)=0;
                        link(i,j,i+1,j)=0;
                    elseif randne==3
                        link(i,j-1,i,j)=0;
                        link(i,j,i,j-1)=0;
                    elseif randne==4
                        link(i,j+1,i,j)=0;
                        link(i,j,i,j+1)=0;
                    end
                elseif i==1&&j>1&&j<celllength2
                    randne=randsrc(1,1,[1:3]);%�����ھ�������һ��
                    if randne==1
                        link(i+1,j,i,j)=0;
                        link(i,j,i+1,j)=0;
                    elseif randne==2
                        link(i,j-1,i,j)=0;
                        link(i,j,i,j-1)=0;
                    elseif randne==3
                        link(i,j+1,i,j)=0;
                        link(i,j,i,j+1)=0;
                    end
                elseif i>1&&i<celllength1&&j==1
                    randne=randsrc(1,1,[1:3]);%�����ھ�������һ��
                    if randne==1
                        link(i-1,j,i,j)=0;
                        link(i,j,i-1,j)=0;
                    elseif randne==2
                        link(i+1,j,i,j)=0;
                        link(i,j,i+1,j)=0;
                    elseif randne==3
                        link(i,j+1,i,j)=0;
                        link(i,j,i,j+1)=0;
                    end
                elseif i==celllength1&&j>1&&j<celllength2
                    randne=randsrc(1,1,[1:3]);%�����ھ�������һ��
                    if randne==1
                        link(i-1,j,i,j)=0;
                        link(i,j,i-1,j)=0;
                    elseif randne==2
                        link(i,j-1,i,j)=0;
                        link(i,j,i,j-1)=0;
                    elseif randne==3
                        link(i,j+1,i,j)=0;
                        link(i,j,i,j+1)=0;
                    end
                elseif i>1&&i<celllength1&&j==celllength2
                    randne=randsrc(1,1,[1:3]);%�����ھ�������һ��
                    if randne==1
                        link(i-1,j,i,j)=0;
                        link(i,j,i-1,j)=0;
                    elseif randne==2
                        link(i+1,j,i,j)=0;
                        link(i,j,i+1,j)=0;
                    elseif randne==3
                        link(i,j-1,i,j)=0;
                        link(i,j,i,j-1)=0;
                    end
                end
                %=========����Զ������===========
                ix=randsrc(1,1,[1:celllength1]);
                if i>1&&i<celllength1&&j>1&&j<celllength2
                    if ix==i-1||ix==i+1
                        jy=randsrc(1,1,[1:j-1,j+1,celllength2]);
                    elseif ix==i
                        jy=randsrc(1,1,[1:j-2,j+2,celllength2]); 
                    else
                        jy=randsrc(1,1,[1:celllength2]);
                    end
                    link(ix,jy,i,j)=1;
                    link(i,j,ix,jy)=1;
                elseif i==1&&j>1&&j<celllength2
                    if ix==i+1
                        jy=randsrc(1,1,[1:j-1,j+1,celllength2]);
                    elseif ix==i
                        jy=randsrc(1,1,[1:j-2,j+2,celllength2]); 
                    else
                        jy=randsrc(1,1,[1:celllength2]);
                    end
                    link(ix,jy,i,j)=1;
                    link(i,j,ix,jy)=1;
                elseif i==celllength1&&j>1&&j<celllength2
                    if ix==i-1
                        jy=randsrc(1,1,[1:j-1,j+1,celllength2]);
                    elseif ix==i
                        jy=randsrc(1,1,[1:j-2,j+2,celllength2]); 
                    else
                        jy=randsrc(1,1,[1:celllength2]);
                    end
                    link(ix,jy,i,j)=1;
                    link(i,j,ix,jy)=1;
                end
            end
    end
end
%======================��������=================
for i=1:celllength1
    for j=1:celllength2
        for m=1:celllength1
            for n=1:celllength2
                if link(m,n,i,j)==1
                   link(i,j,m,n)=1;
                   a=[i m];
                   b=[j n];
                   line(a,b);
                end
            end
        end
    end
end

for i=1:celllength1
    for j=1:celllength2
        axis([0,celllength1+1,0,celllength2+1]);
        text(i,j,' ','color','k');
    end
end

%save network celllength1 celllength2 link
