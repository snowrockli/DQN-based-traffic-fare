%=======================logitģ��================
clc
clear
close all
T_DQN=400;%ǿ��ѧϰ��������
keseip=0.6;%�۸�����ϵ��
gama=0.8;%ǿ��ѧϰ��
epsilon(1)=0.3;%��������
final_epsilon=0.01;%�� epsilon С�ڸ�ֵʱ�����������ѡ����Ϊ
%basic_num=T_DQN/5;
basic_num=100;
beita=0.15;%logitϵ��
min_p1=0;
max_p1=15;
p1(1)=5;%������ʼ�۸�
min_p2=0;
max_p2=7;
p2(1)=1;%������ʼ�۸�
t1=1;%�����г�ʱ��
t2=3;%�����г�ʱ��
c1=2;%�������ʶȳɱ�
c2=3;%�������ʶȳɱ�
b1=4;%������λ�ɱ�
b2=1;%������λ�ɱ�

q1(1)=20;%������ʼ����
q2(1)=20;%������ʼ����
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
    STATE=[p1(t),p2(t),q1(t),q2(t)];%״̬
    for a_n=1:length(A)
        x_S_A(:,a_n)=[STATE,A(a_n)]';
        y_Q(a_n)=sim(net,x_S_A(:,a_n)); %��ͬ������Qֵ
    end
    %===========ѡ����============
    if rand<epsilon(t)
        i_star=randi([1,length(A)],1,1);%���ѡ��
        A_select=A(i_star(1));
    else
        i_star=find(y_Q==max(y_Q));%ѡQֵ����
        %A_select=A(i_star(1));
        A_select=A(i_star(randi([1,length(i_star)],1,1)));
        
%         y_select=y_Q/sum(y_Q);
%         for i=1:length(A)
%             if i==1
%                 y_select_cp(i)=y_select(i);
%             else
%                 y_select_cp(i)=y_select(i)+sum(y_select(1:i-1));
%             end
%         end
%         rand_a=rand;
%         for i=1:length(A)
%             if rand_a<y_select_cp(1)
%                 A_select=A(1);
%                 i_star=1;
%             elseif rand_a<=y_select_cp(i)&&rand_a>y_select_cp(i-1)
%                 A_select=A(i);
%                 i_star=i;
%             end
%         end
    end
    epsilon(t+1)=epsilon(t)-0.001;
    %===============ִ��ѡ��Ķ������õ��µ�״̬==============
    p1(t+1)=p1(t)+0.01*A_select*p1(t);%���������۸�
    
    %p1(t+1)=p1(t)+0.1*A_select;
    
    pp1=p1(t);qq1=q1(t);uu1=p1(t+1);
    pp2=p2(t);qq2=q2(t);
    save data qq1 qq2 keseip pp1 pp2 b1 b2 uu1
    [x,fval2]=fminbnd('f2',min_p2,max_p2);
    p2(t+1)=x;%����ͨ�����ĸ��¼۸�
    q1(t+1)=(q1(1)+q2(1))*exp(-beita*(keseip*p1(t+1)+t1+c1))/(exp(-beita*((keseip*p1(t+1)+t1+c1)))+exp(-beita*((keseip*p2(t+1)+t2+c2))));%������������
   
    q2(t+1)=(q1(1)+q2(1)-q1(t+1));%������������
    %===============���¾���طż���====================
    if p1(t+1)>min_p1&&p1(t+1)<max_p1
        REWARD=(p1(t+1)-b1)*q1(t+1)*0.01;
    elseif p1(t+1)<=min_p1||p1(t+1)>=max_p1
        REWARD=-1;
    end
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
