%=======================agent-logitģ��================
clc
clear
close all
load network
T_DQN=400;%ǿ��ѧϰ��������
keseip=0.6;%�۸�����ϵ��
theita=0.15;%Ч�ø�֪ϵ��
gama=0.8;%ǿ��ѧϰ��
epsilon(1)=0.3;%��������
final_epsilon=0.01;%�� epsilon С�ڸ�ֵʱ�����������ѡ����Ϊ
%basic_num=T_DQN/5;
basic_num=100;
%========ǰ��Ч�ò���=======
pc=0;%����ǿ��
rou=0.5;%���пɿ��Ե�Ҫ��,�������̬��
cl=0.98;%����ˮƽ
dn=5;%��������ȷ���
alfa=0.88;%���չ�̶ܳ�
beita=0.88;%����ƫ�ó̶�
lamada=2.25;%���չ��ϵ��
%==========��ʼ����=========
min_p1=0;
max_p1=15;
p1(1)=5;%������ʼ�۸�
min_p2=0;
max_p2=10;
p2(1)=1;%������ʼ�۸�
t1=1;%�����г�ʱ��
t2=3;%�����г�ʱ��
c1=2;%�������ʶȳɱ�
c2=3;%�������ʶȳɱ�
b1=4;%������λ�ɱ�
b2=1;%������λ�ɱ�

fai1=0.2;%��������ϵ��
fai2=0.3;%��������ϵ��
celllength=celllength1;%������Ⱥ���ģ
risk=rand(celllength,celllength);%����̬�Ⱦ���
risk(risk>0.5)=0.75;%��ʼ���հ���
risk(risk<=0.5)=0.25;%��ʼ���չ��
% risk(risk>0.5)=0.5;%��ʼ���հ���
% risk(risk<=0.5)=0.5;%��ʼ���չ��
risk_num1(1)=numel(risk(risk==0.75))/10;%���հ��������� 
risk_num2(1)=numel(risk(risk==0.25))/10;%���չ��������
for i=1:celllength
    for j=1:celllength
        ww(i,j,1,1)=rand;%��ʼ�����з�ʽѡ�����
        ww(i,j,2,1)=1-ww(i,j,1,1);
        route(i,j,1)=randsrc(1,1,[1:2]);%���䷽ʽ
    end
end
q1(1)=numel(find(route(:,:,1)==1))/10;
q2(1)=numel(find(route(:,:,1)==2))/10;
A=[-5,-1,0,1,5];%������:�۸���ڷ���
D=[];%����طż���
D_train=[];%Qֵѵ������
%====================��ʼ��������============
num_sample=T_DQN;
S_A=rand(num_sample,5);
Q=rand(num_sample,1);
net=newff(S_A',Q',11,{'logsig','purelin','traingd'});
net.trainParam.showWindow = 0;%�Ƿ�չʾ����
for t=1:T_DQN
    %===========����״̬�����Qֵ================
    %STATE=[p1(t),p2(t),q1(t),q2(t),risk_num1(t),risk_num2(t)];%״̬
    STATE=[p1(t),p2(t),q1(t),q2(t)];
    for a_n=1:length(A)
        x_S_A(:,a_n)=[STATE,A(a_n)]';
        y_Q(a_n)=sim(net,x_S_A(:,a_n)); %��ͬ������Qֵ
    end
    %===========ѡ����============
    if rand<epsilon(t)
        i_star=randi([1,length(A)],1,1);%���ѡ��
        A_select=A(i_star);
    else
        i_star=find(y_Q==max(y_Q));%ѡQֵ����
        A_select=A(i_star(randi([1,length(i_star)],1,1)));
    end
    epsilon(t+1)=epsilon(t)-0.001;
    %===============ִ��ѡ��Ķ������õ��µ�״̬===================
    %===========���¼۸�============================
    p1(t+1)=p1(t)+0.01*A_select*p1(t);%���������۸�
    
    
    pp1=p1(t);qq1=q1(t);uu1=p1(t+1);
    pp2=p2(t);qq2=q2(t);
    save data qq1 qq2 keseip pp1 pp2 b1 b2 uu1
    [x,fval2]=fminbnd('f2',min_p2,max_p2);
    p2(t+1)=x;%����ͨ�����ĸ��¼۸�
    %============���¿���===========================
    %===========����С����============
    g1(t)=keseip*p1(t+1)+t1+c1;
    g2(t)=keseip*p2(t+1)+t2+c2;
    XGMG1=abs(fai1*g1(t));XGMG2=abs(fai2*g2(t));XXGMG1(t)=XGMG1;XXGMG2(t)=XGMG2;%���㷽��
    AG1=g1(t)-sqrt(XGMG1)*norminv(0.5+0.5*cl,0,1);BG1=g1(t)+sqrt(XGMG1)*norminv(0.5+0.5*cl,0,1);%������������
    AG2=g2(t)-sqrt(XGMG2)*norminv(0.5+0.5*cl,0,1);BG2=g2(t)+sqrt(XGMG2)*norminv(0.5+0.5*cl,0,1);%������������
    for k=0:dn
        x1(k+1)=AG1+k*(BG1-AG1)/dn;%����С����߽�
        x2(k+1)=AG2+k*(BG2-AG2)/dn;%����С����߽�
    end
    for k=0:dn-1
        xx1(k+1)=AG1+(2*k+1)*(BG1-AG1)/(2*dn);%����С������ֵ
        xx2(k+1)=AG2+(2*k+1)*(BG2-AG2)/(2*dn);%����С������ֵ
        px1(k+1)=normcdf(x1(k+2),g1(t),sqrt(XGMG1))-normcdf(x1(k+1),g1(t),sqrt(XGMG1));%�������ʷֲ�
        px2(k+1)=normcdf(x2(k+2),g2(t),sqrt(XGMG2))-normcdf(x2(k+1),g2(t),sqrt(XGMG2));%�������ʷֲ�
    end
    rp(1)=g1(t)+sqrt(XGMG1)*norminv(rou,0,1);%����Ԥ��
    rp(2)=g2(t)+sqrt(XGMG2)*norminv(rou,0,1);%����Ԥ��
    %=============����Ч��======================
    for i=1:celllength
        for j=1:celllength
            [utility1,utility2]=f_utility(risk,rp,xx1,xx2,px1,px2,dn,alfa,beita,lamada,i,j);
            Futility(i,j,1,t)=utility1;%�����۸�Ч��
            Futility(i,j,2,t)=utility2;%�����۸�Ч��
            ET(i,j,t)=Futility(i,j,route(i,j,t),t);%ʵ��Ч��
            ww(i,j,1,t+1)=exp(theita*Futility(i,j,1,t))/(exp(theita*Futility(i,j,1,t))+exp(theita*Futility(i,j,2,t)));%�������䷽ʽѡ�����
            ww(i,j,2,t+1)=1-ww(i,j,1,t+1);%�������䷽ʽѡ�����
            if rand<=ww(i,j,1,t+1)
                route(i,j,t+1)=1;
            else
                route(i,j,t+1)=2;
            end
        end
    end
    %===============���²��յ�===================================
    for i=1:celllength
        for j=1:celllength
            sx=[];sy=[];sF=[];kk=1;
            for m=1:celllength
                for n=1:celllength
                    if link(m,n,i,j)==1
                        sx(kk)=m;%��¼�ھӺ�����
                        sy(kk)=n;%��¼�ھ�������
                        sF(kk)=ET(m,n,t);%��¼�ھ�ǰ��ֵ
                        kk=kk+1;
                    end
                end
            end
            sx(kk)=i;%��¼���������
            sy(kk)=j;%��¼����������
            sF(kk)=ET(i,j,t);
            %=================Ѱ��ǰ��Ч������============
            if numel(sF)==0
                risk(i,j)=risk(i,j);
            else
                kstar=find(sF==max(sF));
                %=================���²��յ�=============
                if rand<pc
                    risk(i,j)=risk(sx(kstar(1)),sy(kstar(1)));
                end
                %risk(i,j)=(1-pc)*risk(i,j)+pc*risk(sx(kstar(1)),sy(kstar(1)));
            end
        end
    end
    %=================��������============
    q1(t+1)=numel(find(route(:,:,t+1)==1))/10;
    q2(t+1)=numel(find(route(:,:,t+1)==2))/10;
    risk_num1(t+1)=numel(risk(risk==0.75))/10;%���հ��������� 
    risk_num2(t+1)=numel(risk(risk==0.25))/10;%���չ��������
    %===============���¾���طż���====================
    if p1(t+1)>min_p1&&p1(t+1)<max_p1
        REWARD=(p1(t+1)-b1)*q1(t+1)*0.01;
    elseif p1(t+1)<=min_p1||p1(t+1)>=max_p1
        REWARD=-1;
    end
    %NEW_STATE=[p1(t+1),p2(t+1),q1(t+1),q2(t+1),risk_num1(t),risk_num2(t)];
    NEW_STATE=[p1(t+1),p2(t+1),q1(t+1),q2(t+1)];
    D(t,:)=[STATE,A_select,REWARD,NEW_STATE];
    %===============�Ӿ���طż�����ѡȡһЩ����====================
    D_train_temp=[];
    if t>=basic_num%�жϼ������������Ƿ��㹻
        c=randperm(numel(1:t));%���´���˳��
        m=basic_num;%ѡ��m������
        for i=1:m
            D_train_temp(i,1:size(STATE,2)+1)=D(c(i),1:size(STATE,2)+1);
            STATE_next=D(c(i),size(STATE,2)+3:end);%��һ��״̬
            for a_n=1:length(A)
                x_next_S_A(:,a_n)=[STATE_next,A(a_n)]';
                y_next_Q(a_n)=sim(net,x_next_S_A(:,a_n)); %��ͬ������Qֵ
            end
            max_Q=max(y_next_Q);
            D_train_temp(i,size(STATE,2)+2)=D(c(i),size(STATE,2)+2)+gama*max_Q;%����Qֵ
        end
        %===============����Qֵѵ����=================
        [D_train_num,~]=size(D_train);
        if isempty(D_train)==1
            D_train=D_train_temp;
        else
            for i=1:m
                for j=1:D_train_num
                    if isequal(D_train_temp(i,1:size(STATE,2)+1),D_train(j,1:size(STATE,2)+1))==1
                        D_train(j,size(STATE,2)+2)=D_train_temp(i,size(STATE,2)+2);
                    end
                end
            end
            for i=1:m
                if ismember(D_train_temp(i,1:size(STATE,2)+1),D_train(:,1:size(STATE,2)+1),'rows')==0
                    [D_train_num,~]=size(D_train);
                    D_train(D_train_num+1,:)=D_train_temp(i,:);
                end
            end
        end
        %D_train=unique(D_train,'rows');
        net = train(net,D_train(:,1:size(STATE,2)+1)',D_train(:,size(STATE,2)+2)');
        net.trainParam.goal = 1e-5;
        net.trainParam.epochs = 300;
        net.trainParam.lr = 0.05;
        net.trainParam.showWindow = 0;%�Ƿ�չʾ����
    end   
end


























