%=======================logit模型================
clc
clear
close all
T_DQN=400;%强化学习迭代次数
keseip=0.6;%价格敏感系数
gama=0.8;%强化学习率
epsilon(1)=0.3;%搜索策略
final_epsilon=0.01;%当 epsilon 小于该值时，将不在随机选择行为
%basic_num=T_DQN/5;
basic_num=100;
beita=0.15;%logit系数
min_p1=0;
max_p1=15;
p1(1)=5;%地铁初始价格
min_p2=0;
max_p2=7;
p2(1)=1;%公交初始价格
t1=1;%地铁行程时间
t2=3;%公交行程时间
c1=2;%地铁舒适度成本
c2=3;%公交舒适度成本
b1=4;%地铁单位成本
b2=1;%公交单位成本

q1(1)=20;%地铁初始人数
q2(1)=20;%公交初始人数
A=[-5,-1,0,1,5];%动作集:价格调节幅度
D=[];%经验回放集合
D_train=[];%Q值训练集合
%====================初始化神经网络============
num_sample=T_DQN;
S_A=rand(num_sample,5);
Q=rand(num_sample,1);
net=newff(S_A',Q',11,{'logsig','purelin','traingd'});
net.trainParam.showWindow = 0;%是否展示窗口

for t=1:T_DQN
    %===========输入状态，输出Q值================
    STATE=[p1(t),p2(t),q1(t),q2(t)];%状态
    for a_n=1:length(A)
        x_S_A(:,a_n)=[STATE,A(a_n)]';
        y_Q(a_n)=sim(net,x_S_A(:,a_n)); %不同动作的Q值
    end
    %===========选择动作============
    if rand<epsilon(t)
        i_star=randi([1,length(A)],1,1);%随机选择
        A_select=A(i_star(1));
    else
        i_star=find(y_Q==max(y_Q));%选Q值最大的
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
    %===============执行选择的动作，得到新的状态==============
    p1(t+1)=p1(t)+0.01*A_select*p1(t);%地铁调整价格
    
    %p1(t+1)=p1(t)+0.1*A_select;
    
    pp1=p1(t);qq1=q1(t);uu1=p1(t+1);
    pp2=p2(t);qq2=q2(t);
    save data qq1 qq2 keseip pp1 pp2 b1 b2 uu1
    [x,fval2]=fminbnd('f2',min_p2,max_p2);
    p2(t+1)=x;%公交通过博弈更新价格
    q1(t+1)=(q1(1)+q2(1))*exp(-beita*(keseip*p1(t+1)+t1+c1))/(exp(-beita*((keseip*p1(t+1)+t1+c1)))+exp(-beita*((keseip*p2(t+1)+t2+c2))));%地铁更新流量
   
    q2(t+1)=(q1(1)+q2(1)-q1(t+1));%公交更新流量
    %===============更新经验回放集合====================
    if p1(t+1)>min_p1&&p1(t+1)<max_p1
        REWARD=(p1(t+1)-b1)*q1(t+1)*0.01;
    elseif p1(t+1)<=min_p1||p1(t+1)>=max_p1
        REWARD=-1;
    end
    NEW_STATE=[p1(t+1),p2(t+1),q1(t+1),q2(t+1)];
    D(t,:)=[STATE,A_select,REWARD,NEW_STATE];
    %===============从经验回放集合中选取一些样本====================
    D_train_temp=[];
    if t>=basic_num%判断记忆池里的数据是否足够
        c=randperm(numel(1:t));%重新打乱顺序
        m=basic_num;%选出m个经验
        for i=1:m
            D_train_temp(i,1:size(STATE,2)+1)=D(c(i),1:size(STATE,2)+1);
            STATE_next=D(c(i),size(STATE,2)+3:end);%下一个状态
            for a_n=1:length(A)
                x_next_S_A(:,a_n)=[STATE_next,A(a_n)]';
                y_next_Q(a_n)=sim(net,x_next_S_A(:,a_n)); %不同动作的Q值
            end
            max_Q=max(y_next_Q);
            D_train_temp(i,size(STATE,2)+2)=D(c(i),size(STATE,2)+2)+gama*max_Q;%更新Q值
        end
        %===============更新Q值训练集=================
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
        net.trainParam.showWindow = 0;%是否展示窗口
    end
end
