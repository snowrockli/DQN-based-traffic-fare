function [u1,u2]=f_utility(risk,rp,xx1,xx2,px1,px2,dn,alfa,beita,lamada,i,j)
%=====================设定参照点===========
u0(i,j)=min(rp)+risk(i,j)*(max(rp)-min(rp));
%==========方式1前景效用=====================
for k=1:dn
    if xx1(k)<=u0(i,j)
        vg1(k)=(u0(i,j)-xx1(k))^alfa;%收益
    else
        vg1(k)=-lamada*(xx1(k)-u0(i,j))^beita;%损失
    end
end
if xx1(1)>=u0(i,j)
    kstar=1;
elseif xx1(dn)<=u0(i,j)
    kstar=dn;
end
for k=2:dn
    if vg1(k-1)>=0&&vg1(k)<0
        kstar=k;%记录参考点之后的第一个位置
    end
end
Futility11=0;
for k=kstar:dn
    wp11=wwww(sum(px1(k:dn)))-wwww(sum(px1(k+1:dn)));
    Futility11=Futility11+vg1(k)*wp11;%损失
end
Futility12=0;
for k=1:kstar-1
    wp12=wwww(sum(px1(1:k)))-wwww(sum(px1(1:k-1)));
    Futility12=Futility12+vg1(k)*wp12;%收益
end
u1=Futility11+Futility12;%前景效用
%===================方式2前景效用==============
for k=1:dn
    if xx2(k)<=u0(i,j)
        vg2(k)=(u0(i,j)-xx2(k))^alfa;%收益
    else
        vg2(k)=-lamada*(xx2(k)-u0(i,j))^beita;%损失
    end
end
if xx2(1)>=u0(i,j)
    kstar=1;
elseif xx2(dn)<=u0(i,j)
    kstar=dn;
end
for k=2:dn
    if vg2(k-1)>=0&&vg2(k)<0
        kstar=k;%记录参考点之后的第一个位置
    end
end
Futility21=0;
for k=kstar:dn
    wp21=wwww(sum(px2(k:dn)))-wwww(sum(px2(k+1:dn)));
    Futility21=Futility21+vg2(k)*wp21;%损失
end
Futility22=0;
for k=1:kstar-1
    wp22=wwww(sum(px2(1:k)))-wwww(sum(px2(1:k-1)));
    Futility22=Futility22+vg2(k)*wp22;%收益
end
u2=Futility21+Futility22;%前景效用

end
