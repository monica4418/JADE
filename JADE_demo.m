function [optValue,bestP,bestfitnessG,PJILU]=JADE_demo(NP,Gmax,Pmin,Pmax,Dim,fun,fun_flag)
 
    %addpath('../CauchyFunction');
 
    %%��ʼ��
%     Pmax=ma;%�����Ͻ�
%     Pmin=mi;%�����½�
%     NP=100;%��Ⱥ��ģ
%     Dim=30;%����ά��
    G=1;%���õ���������ǰ����������
    uCR=0.5;%��ʼ���������
    uF=0.5;%��ʼ����������
    p0=0.06;
    top=p0*NP;%ÿ�������ŵ�top��
    A=[];%��ʼ�鵵��ȺΪ�ռ�
    t=1;%��¼�鵵��ȺA�еĸ������
   % Gmax=gen;%����������
    %index=type;%���Ժ�������
    
%     index=input('Input the index of TestFunction:');
 popsize=NP;
 P = repmat(Pmin, popsize, 1) + rand(popsize, Dim) .* (repmat(Pmax-Pmin, popsize, 1)); %��Ⱥ��ʼ��
    %%�����ѭ��
    for i=1:NP
          fitnessP(i)=fun(P(i,:)',fun_flag);%������Ⱥ������Ӧֵ
    end
    while G<Gmax
        Scr=[];%��ʼ�ɹ��μӱ���ĸ���Ľ������Ϊ�ռ�
        Sf=[];%��ʼ�ɹ��μӱ���ĸ������������Ϊ�ռ�
        n1=1;%��¼Scr�е�Ԫ�ظ���
        n2=1;%��¼Sf�е�Ԫ�ظ���
 
        for i=1:NP
            %�Ե�G����ÿ��������㽻����ʺ���������
            CR(i)=normrnd(uCR,0.2);%������̬�ֲ�
            F(i)=cauchyrnd(uF,0.2);%���ӿ����ֲ�0.5��0.3
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
 
        %���ݸ�����Ӧֵ�������򣬵õ���õ�ǰtop������
        [~,indexSortP]=sort(fitnessP);
        for i=1:top
           bestTopP(i,:)=P(indexSortP(i),:); 
        end
 
        %%�������
        for i=1:NP   
            %��top�����������ѡ��һ����ΪXpbest
            k0=randperm(top,1);
            Xpbest=bestTopP(k0,:);
            %�ӵ�ǰ��ȺP�����ѡ��P1
            k1=randi(NP);
            P1=P(k1,:);
            while (k1==i||k1==k0)
                k1=randi(NP);
                P1=P(k1,:); 
            end
 
            %��P��A�����ѡ��P2
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
 
        %%�������
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
 
        %%�߽紦��
  %      for i=1:NP
           for j=1:Dim
              while (U(i,j)>Pmax||U(i,j)<Pmin)
                 U(i,j)=(Pmax-Pmin)*rand+Pmin; 
              end
           end
 %       end
 
        %%ѡ�����
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
 
        %�жϹ鵵��ȺA�Ĺ�ģ�Ƿ���NP֮�ڣ������ڣ�������Ƴ�����ʹ���ģ����NP
        [tA,~]=size(A);
        if tA>NP
            nRem=tA-NP;
            k4=randperm(tA,nRem);
            A(k4,:)=[]; 
            [tA,~]=size(A);
            t=tA+1;
        end
 
        %����Ӧ����������uCR��uF
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
