clc
clear all

%01 算法输入+预处理模块
values = [10 3 2 8 9 7 4 5 9 4 5 6 5 4 8 3 4 9 2 7 5 3 6 8 2 4 9 5];
sizes  = [2 9 8 3 2 10 2 3 4 5 3 4 2 8 6 4 2 3 8 6 3 5 2 4 1 8 6 7];
max_size=30;%跟科研相关 约束条件
n=size(values,2);%只拿列数

%02 算法控制模块
PopSize=100;
Pro_Initial=0.55; %初代生成概率
Pro_Cross=0.7;
Pro_Muta=0.05;
elite_num=floor(PopSize*0.05);
st=200;%迭代次数
%03 初代生成
%Fitness=zeros(PopSize,1);
InitialPop=zeros(PopSize,n);%预设大小
k=0;%成功次数
KK=0;%尝试次数
size_all=[];%所有随机出来的sum_size历史记录
while k<PopSize
    Indi=rand(1,n)<Pro_Initial; %加速版，随机生成一列0~1的数然后做逻辑判断
    sum_size=sum(sizes.*Indi);
    KK=KK+1;%多做了一次尝试
    size_all=[size_all,sum_size];%本身在循环里面，第一位数在不断变长
    if sum_size<=max_size
    % max_size_sat=0;
    % while max_size_sat==0
    %     Indi=zeros(1,n);
    %     for i=1:1:n
    %         if rand()<0.3  %数学+科研
    %             Indi(i)=1;
    %         end
    %     end   
    %     sum_size=sum(sizes.*Indi);
    % 
    %     %数据量很大时用下面的方法求sum_size
    %     %size_all = zeros(1, 1e6); % 预先开100万长度，足够用
    %     %idx = 1; % 记录位置
    %     %size_all(idx)=sum_size;
    %     %idx =idx+1;
    %     %size_all = size_all(1:idx-1);
    % 
    %     KK=KK+1;%多做了一次尝试
    %     size_all=[size_all,sum_size];%本身在循环里面，第一位数在不断变长
    %     if sum_size<=max_size
    %         max_size_sat=1;
    %     end
    % end
      k=k+1;%合规次数加
      InitialPop(k,:)=Indi;%合规的记录下来
     %Fitness(k,:)=sum(Indi.*values);
    end
end
[k,KK] %循环之外，只输出最终结果

%04 选择：计算fitness（总Value值）、通过概率选择亲代
for cs=1:1:st % st次大循环


    % Fitness=zeros(PopSize,1);
    % for i=1:1:PopSize
    %   Fitness(i,:)=sum(InitialPop(i,:).*values);
    % end
    fitness=InitialPop*values';%点积能实现优化
    [sorted_fit,idx]=sort(fitness,'descend');%idx=排序后的原来位置编号
    sorted_Population=InitialPop(idx,:);% 用idx把种群也按适应度从高到低排好
    normal_num=PopSize-elite_num;   %普通种群数量
    elite_pop=sorted_Population(1:elite_num,:);%精英种群
    normal_pop=sorted_Population(elite_num+1:PopSize,:);%普通种群
    normal_fit=sorted_fit(elite_num+1:end);  %普通个体的分数
    sum_normal=sum(normal_fit);
    prob_normal=normal_fit/sum_normal;    %只算普通个体的概率
    selected_idx=randsample(normal_num,normal_num,true,prob_normal);%normal_num里，抽normal_num个，按prob_normal概率
    parents=normal_pop(selected_idx,:); %适应度高的会被抽中多次
    parents_fitness=parents*values';
    
 %05 交叉：怎么挑、怎么交叉、怎么保证合法性
    p_tmp=parents;
    [normal_num,gene_length]=size(p_tmp);
    child_pop=[];%子代种群
    
    while size(child_pop,1)<normal_num
        cross_points=randsample(1:normal_num,2);%randsample(范围, 抽取个数, 是否有放回)
        cross1=min(cross_points);  %小的在前
        cross2=max(cross_points);  %大的在后
        x01=p_tmp(cross1,:);%父1
        x02=p_tmp(cross2,:);%父2
        retry = 0;
        while all(x01==x02)&&retry<3%如果父母完全相同，最多重试3次 中间是and逻辑
            cross_points=randsample(1:normal_num,2);
            cross1=min(cross_points);
            cross2=max(cross_points);
            x01=p_tmp(cross1,:);
            x02=p_tmp(cross2,:);
            retry=retry+1;
        end
        if rand()<Pro_Cross
          dd=randi([1,gene_length-1]);%截断点
          x03=[x01(1:dd),x02(dd+1:end)];%子1
          x04=[x02(1:dd),x01(dd+1:end)];%子2
          sumx03=sum(sizes.*x03);
          sumx04=sum(sizes.*x04);
          while sumx03>max_size%是否超约束
              values03=x03.*values;
              % values03=zeros(1,gene_length); 
              % for i=1:1:gene_length
              %     if x03(i)==1
              %         values03(i)=values(i);%对确定背包的附上价值
              %     end
              % end
              minpos03=find(values03==min(values03(values03~=0)),1);%找非零价值的最小值所在位置
              x03(minpos03)=0;%改变价值最低的
              sumx03=sum(sizes.*x03);%更新求和项
          end
          while sumx04>max_size
              values04=x04.*values;
              % values04=zeros(1,gene_length); %尽量不用for循环，太慢了
              % for i=1:1:gene_length
              %     if x04(i)==1
              %         values04(i)=values(i);
              %     end
              % end
              minpos04=find(values04==min(values04(values04~=0)),1);
              x04(minpos04)=0;
              sumx04=sum(sizes.*x04);
          end
        else 
          x03=x01;
          x04=x02;
        end
        if isempty(child_pop)||~any(all(child_pop==x03,2))%判断是否是空集，中间是或关系，判断是否有行相等，~为取反操作
         child_pop=[child_pop;x03];
        end
        if isempty(child_pop)||~any(all(child_pop==x04,2))
         child_pop=[child_pop;x04];
        end
    end
    child_pop=child_pop(1:normal_num,:);
    
 %06 变异
    for i=1:size(child_pop,1)
        x=child_pop(i,:);  %取出第i个个体
        for j=1:gene_length
            if rand()<Pro_Muta
                x(j)=1-x(j);  %翻转：0变1，1变0
            end
        end
        sumx=sum(sizes.*x);
        while sumx>max_size %总量约束，逻辑同上
            values_x=x.*values;
            minpos=find(values_x==min(values_x(values_x~=0)),1);
            x(minpos)=0;
            sumx=sum(sizes.*x);
        end
        child_pop(i,:) = x;
    end
    InitialPop=[elite_pop;child_pop];

%07 循环、记录、出图
   big_value=max(InitialPop*values');
   fprintf('第%d代，最优适应度：%.2f\n',cs,big_value);
end