%=======================agent-logit模型================
clc
clear
close all
load network
T_DQN=400;%强化学习迭代次数
keseip=0.6;%价格敏感系数
theita=0.15;%效用感知系数
gama=0.8;%强化学习率
epsilon(1)=0.3;%搜索策略
final_epsilon=0.01;%当 epsilon 小于该值时，将不在随机选择行为
%basic_num=T_DQN/5;
basic_num=100;
%========前景效用参数=======
pc=0;%交互强度
rou=0.5;%出行可靠性的要求,整体风险态度
cl=0.98;%置信水平
dn=5;%置信区间等分数
alfa=0.88;%风险规避程度
beita=0.88;%风险偏好程度
lamada=2.25;%风险规避系数
%==========初始参数=========
min_p1=0;
max_p1=15;
p1(1)=5;%地铁初始价格
min_p2=0;
max_p2=10;
p2(1)=1;%公交初始价格
t1=1;%地铁行程时间
t2=3;%公交行程时间
c1=2;%地铁舒适度成本
c2=3;%公交舒适度成本
b1=4;%地铁单位成本
b2=1;%公交单位成本

fai1=0.2;%地铁变异系数
fai2=0.3;%公交变异系数
celllength=celllength1;%出行者群体规模
risk=rand(celllength,celllength);%风险态度矩阵
risk(risk>0.5)=0.75;%初始风险爱好
risk(risk<=0.5)=0.25;%初始风险规避
% risk(risk>0.5)=0.5;%初始风险爱好
% risk(risk<=0.5)=0.5;%初始风险规避
risk_num1(1)=numel(risk(risk==0.75))/10;%风险爱好者人数 
risk_num2(1)=numel(risk(risk==0.25))/10;%风险规避者人数
for i=1:celllength
    for j=1:celllength
        ww(i,j,1,1)=rand;%初始化出行方式选择概率
        ww(i,j,2,1)=1-ww(i,j,1,1);
        route(i,j,1)=randsrc(1,1,[1:2]);%运输方式
    end
end
q1(1)=numel(find(route(:,:,1)==1))/10;
q2(1)=numel(find(route(:,:,1)==2))/10;
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
    %STATE=[p1(t),p2(t),q1(t),q2(t),risk_num1(t),risk_num2(t)];%状态
    STATE=[p1(t),p2(t),q1(t),q2(t)];
    for a_n=1:length(A)
        x_S_A(:,a_n)=[STATE,A(a_n)]';
        y_Q(a_n)=sim(net,x_S_A(:,a_n)); %不同动作的Q值
    end
    %===========选择动作============
    if rand<epsilon(t)
        i_star=randi([1,length(A)],1,1);%随机选择
        A_select=A(i_star);
    else
        i_star=find(y_Q==max(y_Q));%选Q值最大的
        A_select=A(i_star(randi([1,length(i_star)],1,1)));
    end
    epsilon(t+1)=epsilon(t)-0.001;
    %===============执行选择的动作，得到新的状态===================
    %===========更新价格============================
    p1(t+1)=p1(t)+0.01*A_select*p1(t);%地铁调整价格
    
    
    pp1=p1(t);qq1=q1(t);uu1=p1(t+1);
    pp2=p2(t);qq2=q2(t);
    save data qq1 qq2 keseip pp1 pp2 b1 b2 uu1
    [x,fval2]=fminbnd('f2',min_p2,max_p2);
    p2(t+1)=x;%公交通过博弈更新价格
    %============更新客流===========================
    %===========划分小区间============
    g1(t)=keseip*p1(t+1)+t1+c1;
    g2(t)=keseip*p2(t+1)+t2+c2;
    XGMG1=abs(fai1*g1(t));XGMG2=abs(fai2*g2(t));XXGMG1(t)=XGMG1;XXGMG2(t)=XGMG2;%计算方差
    AG1=g1(t)-sqrt(XGMG1)*norminv(0.5+0.5*cl,0,1);BG1=g1(t)+sqrt(XGMG1)*norminv(0.5+0.5*cl,0,1);%计算置信区间
    AG2=g2(t)-sqrt(XGMG2)*norminv(0.5+0.5*cl,0,1);BG2=g2(t)+sqrt(XGMG2)*norminv(0.5+0.5*cl,0,1);%计算置信区间
    for k=0:dn
        x1(k+1)=AG1+k*(BG1-AG1)/dn;%地铁小区间边界
        x2(k+1)=AG2+k*(BG2-AG2)/dn;%公交小区间边界
    end
    for k=0:dn-1
        xx1(k+1)=AG1+(2*k+1)*(BG1-AG1)/(2*dn);%地铁小区间中值
        xx2(k+1)=AG2+(2*k+1)*(BG2-AG2)/(2*dn);%公交小区间中值
        px1(k+1)=normcdf(x1(k+2),g1(t),sqrt(XGMG1))-normcdf(x1(k+1),g1(t),sqrt(XGMG1));%地铁概率分布
        px2(k+1)=normcdf(x2(k+2),g2(t),sqrt(XGMG2))-normcdf(x2(k+1),g2(t),sqrt(XGMG2));%公交概率分布
    end
    rp(1)=g1(t)+sqrt(XGMG1)*norminv(rou,0,1);%地铁预算
    rp(2)=g2(t)+sqrt(XGMG2)*norminv(rou,0,1);%公交预算
    %=============计算效用======================
    for i=1:celllength
        for j=1:celllength
            [utility1,utility2]=f_utility(risk,rp,xx1,xx2,px1,px2,dn,alfa,beita,lamada,i,j);
            Futility(i,j,1,t)=utility1;%地铁价格效用
            Futility(i,j,2,t)=utility2;%公交价格效用
            ET(i,j,t)=Futility(i,j,route(i,j,t),t);%实际效用
            ww(i,j,1,t+1)=exp(theita*Futility(i,j,1,t))/(exp(theita*Futility(i,j,1,t))+exp(theita*Futility(i,j,2,t)));%更新运输方式选择概率
            ww(i,j,2,t+1)=1-ww(i,j,1,t+1);%更新运输方式选择概率
            if rand<=ww(i,j,1,t+1)
                route(i,j,t+1)=1;
            else
                route(i,j,t+1)=2;
            end
        end
    end
    %===============更新参照点===================================
    for i=1:celllength
        for j=1:celllength
            sx=[];sy=[];sF=[];kk=1;
            for m=1:celllength
                for n=1:celllength
                    if link(m,n,i,j)==1
                        sx(kk)=m;%记录邻居横坐标
                        sy(kk)=n;%记录邻居纵坐标
                        sF(kk)=ET(m,n,t);%记录邻居前景值
                        kk=kk+1;
                    end
                end
            end
            sx(kk)=i;%记录自身横坐标
            sy(kk)=j;%记录自身纵坐标
            sF(kk)=ET(i,j,t);
            %=================寻找前景效用最大的============
            if numel(sF)==0
                risk(i,j)=risk(i,j);
            else
                kstar=find(sF==max(sF));
                %=================更新参照点=============
                if rand<pc
                    risk(i,j)=risk(sx(kstar(1)),sy(kstar(1)));
                end
                %risk(i,j)=(1-pc)*risk(i,j)+pc*risk(sx(kstar(1)),sy(kstar(1)));
            end
        end
    end
    %=================更新运量============
    q1(t+1)=numel(find(route(:,:,t+1)==1))/10;
    q2(t+1)=numel(find(route(:,:,t+1)==2))/10;
    risk_num1(t+1)=numel(risk(risk==0.75))/10;%风险爱好者人数 
    risk_num2(t+1)=numel(risk(risk==0.25))/10;%风险规避者人数
    %===============更新经验回放集合====================
    if p1(t+1)>min_p1&&p1(t+1)<max_p1
        REWARD=(p1(t+1)-b1)*q1(t+1)*0.01;
    elseif p1(t+1)<=min_p1||p1(t+1)>=max_p1
        REWARD=-1;
    end
    %NEW_STATE=[p1(t+1),p2(t+1),q1(t+1),q2(t+1),risk_num1(t),risk_num2(t)];
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
