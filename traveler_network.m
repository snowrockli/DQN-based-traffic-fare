%======================生成agent复杂网络============================
clc
clear
close all
%===============小世界网络参数=============
celllength=20;
celllength1=celllength;
celllength2=celllength;
pneighborlink=1;%与邻居的连接概率
pcut=0.5;%断点重连概率
link=zeros(celllength1,celllength2,celllength1,celllength2);
%===============建立出行者小世界网络================
%==============建立出行者之间的联系====================
for i=1:celllength1
    for j=1:celllength2
        for m=1:celllength1
            for n=1:celllength2
                link(i,j,i,j)=0;
                if link(m,n,i,j)==0
                    if rand<=pneighborlink%以一定的概率建立与邻居的链接
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
%===========断点重连===========================
for i=1:celllength1
    for j=1:celllength2
            if rand<pcut
                %==========断开连接=============
                if i>1&&i<celllength1&&j>1&&j<celllength2
                    randne=randsrc(1,1,[1:4]);%四个邻居里面挑一个
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
                    randne=randsrc(1,1,[1:3]);%三个邻居里面挑一个
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
                    randne=randsrc(1,1,[1:3]);%三个邻居里面挑一个
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
                    randne=randsrc(1,1,[1:3]);%三个邻居里面挑一个
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
                    randne=randsrc(1,1,[1:3]);%三个邻居里面挑一个
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
                %=========建立远程连接===========
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
%======================绘制连接=================
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
