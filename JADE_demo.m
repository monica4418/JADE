function [optValue,bestP,bestfitnessG,PJILU]=JADE_demo(NP,Gmax,Pmin,Pmax,Dim,fun,fun_flag)
 
    %addpath('../CauchyFunction');
 
    %%初始化
%     Pmax=ma;%搜索上界
%     Pmin=mi;%搜索下界
%     NP=100;%种群规模
%     Dim=30;%个体维数
    G=1;%设置迭代器（当前迭代代数）
    uCR=0.5;%初始化交叉概率
    uF=0.5;%初始化缩放因子
    p0=0.06;
    top=p0*NP;%每代中最优的top个
    A=[];%初始归档种群为空集
    t=1;%记录归档种群A中的个体个数
   % Gmax=gen;%最大迭代次数
    %index=type;%测试函数类型
    
%     index=input('Input the index of TestFunction:');
 popsize=NP;
 P = repmat(Pmin, popsize, 1) + rand(popsize, Dim) .* (repmat(Pmax-Pmin, popsize, 1)); %种群初始化
    %%总体大循环
    for i=1:NP
          fitnessP(i)=fun(P(i,:)',fun_flag);%计算种群个体适应值
    end
    while G<Gmax
        Scr=[];%初始成功参加变异的个体的交叉概率为空集
        Sf=[];%初始成功参加变异的个体的缩放因子为空集
        n1=1;%记录Scr中的元素个数
        n2=1;%记录Sf中的元素个数
 
        for i=1:NP
            %对第G代的每个个体计算交叉概率和缩放因子
            CR(i)=normrnd(uCR,0.2);%服从正态分布
            F(i)=cauchyrnd(uF,0.2);%服从柯西分布0.5和0.3
            while (CR(i)>1||CR(i)<0)
                CR(i)=normrnd(uCR,0.2);
            end
            if (F(i)>1)
                F(i)=1;
            end
            while (F(i)<=0)
                F(i)=cauchyrnd(uF,0.2);
            end
        end
        [fitnessBestP,indexBestP]=min(fitnessP);
        bestP=P(indexBestP,:);
 
        %根据个体适应值进行排序，得到最好的前top个个体
        [~,indexSortP]=sort(fitnessP);
        for i=1:top
           bestTopP(i,:)=P(indexSortP(i),:); 
        end
 
        %%变异操作
        for i=1:NP   
            %从top个个体中随机选出一个作为Xpbest
            k0=randperm(top,1);
            Xpbest=bestTopP(k0,:);
            %从当前种群P中随机选出P1
            k1=randi(NP);
            P1=P(k1,:);
            while (k1==i||k1==k0)
                k1=randi(NP);
                P1=P(k1,:); 
            end
 
            %从P∪A中随机选出P2
            PandA=[P;A];
            [num,~]=size(PandA);
            k2=randi(num);
            P2=PandA(k2,:);
            while (k2==i||k2==k0||k2==k1)
                k2=randi(num);
                P2=PandA(k2,:); 
            end
 
            V(i,:)=P(i,:)+F(i).*(Xpbest-P(i,:))+F(i).*(P1-P2);   
       %      end
  %      end
 
        %%交叉操作
  %      for i=1:NP
            jrand=randi([1,Dim]); 
            for j=1:Dim
                k3=rand;
                if(k3<=CR(i)||j==jrand)
                    U(i,j)=V(i,j);
                else
                    U(i,j)=P(i,j);
                end
            end
  %      end
 
        %%边界处理
  %      for i=1:NP
           for j=1:Dim
              while (U(i,j)>Pmax||U(i,j)<Pmin)
                 U(i,j)=(Pmax-Pmin)*rand+Pmin; 
              end
           end
 %       end
 
        %%选择操作
  %      for i=1:NP
            fitnessU(i)=fun(U(i,:)',fun_flag);
            if(fitnessU(i)<fitnessP(i))
                A(t,:)=P(i,:);
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
                Scr(n1)=CR(i);
                Sf(n2)=F(i);
                t=t+1;
                n1=n1+1;
                n2=n2+1;
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
               bestfitnessG(G)=fitnessBestP;
               G=G+1;
     
        end
 
        %判断归档种群A的规模是否在NP之内，若大于，则随机移除个体使其规模保持NP
        [tA,~]=size(A);
        if tA>NP
            nRem=tA-NP;
            k4=randperm(tA,nRem);
            A(k4,:)=[]; 
            [tA,~]=size(A);
            t=tA+1;
        end
 
        %自适应参数，更新uCR和uF
        c=0.1;
        [~,ab]=size(Scr);
        if ab~=0
            newSf=(sum(Sf.^2))/(sum(Sf));
            uCR=(1-c)*uCR+c.*mean(Scr);
            uF=(1-c)*uF+c.*newSf;
        end
 

    end
 
    optValue=fitnessBestP;
        PJILU=P;
bestfitnessG=bestfitnessG(500:500:Gmax);
end
