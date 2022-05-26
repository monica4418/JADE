  close all
  clear all
  clc
  %progressbar('1','2');
% mex cec14_func.cpp -DWINDOWS
  fobj=str2func('cec14_func');
  lb=-100;
  ub=100;
  format long%����С���㹻��
  func_NUM=30;
  iter=21;%ѭ��ִ�д���21  
  SearchAgents_no=500;
  dim=30;
  FEsmax=dim*10000;
  Max_iteration=FEsmax/SearchAgents_no;
  %����������������Ԥ�ȷ������ռ䡪��������������
%   PSO_cg_curve_all=zeros(iter,Max_iteration);
%   PSO_cg_curve_save=zeros(func_NUM*iter,Max_iteration);
%   PSO_best_pos_save=zeros(func_NUM*iter,dim);
%   SAVE_PSO_time_cost_save=zeros(func_NUM,iter);
%   SAVE_PSO_score_save=1000*ones(func_NUM,iter);
%   SAVE_PSO_mean_cg_curve=zeros(func_NUM,Max_iteration);
%   PSOscore_mean_best_Worst_std_time=zeros(func_NUM,5);
  %____________________end
 
for kk=1:30%30�����Ժ���
    func_num=kk; 
   
    parfor repeat=1:iter
        
        
        %������������������������mainGEDWOA;
        tic;                                 %GDSPSO_CEC 
        [Best_score,Best_pos,PSO_cg_curve]=JADE_demo(SearchAgents_no,FEsmax,lb,ub,dim,fobj,func_num);%FEsmax
        TotalTime=toc;
        Best_score=Best_score-kk*100;
%         if Best_score<10e-8
%             Best_score=0;
%         end
        display(['The running time of F',num2str(kk),' is: ',num2str(TotalTime),'    The best  score of F',num2str(kk),' is : ',num2str(Best_score)]);
        %������������������������mainGEDWOA��������end
        PSO_cg_curve_all(repeat,:)=PSO_cg_curve-kk*100;%�ظ�repeat�ε���������
       % PSO_best_pos_save(iter*(kk-1)+repeat,:)=Best_pos;
        SAVE_PSO_score_save(kk,repeat)=Best_score;
        SAVE_PSO_time_cost_save(kk,repeat)=TotalTime;
      %  p_ZHI_TOTAL(:,:,kk)=PPPP;
        
     
    end  
    SAVE_PSO_mean_cg_curve(kk,:)=mean(PSO_cg_curve_all,1);%�ظ�repeat�ε��������߾�ֵ baocun
  
end
    PSOscore_mean_best_Worst_std_time(:,1)=mean(SAVE_PSO_score_save,2);
    PSOscore_mean_best_Worst_std_time(:,2)=min(SAVE_PSO_score_save,[],2);
    PSOscore_mean_best_Worst_std_time(:,3)=max(SAVE_PSO_score_save,[],2);
    PSOscore_mean_best_Worst_std_time(:,4)=std(SAVE_PSO_score_save,[],2);%���Ž�ñ�׼�����s=sqrt([(x1-EX)^2+(x2-EX)^2....(xn-EX)^2]/(n-1))
    PSOscore_mean_best_Worst_std_time(:,5)=mean(SAVE_PSO_time_cost_save,2);

%save('PSO190509');
PSO_mean_b_w_sd_time= PSOscore_mean_best_Worst_std_time';

SAVE_date=reshape(PSO_mean_b_w_sd_time,30*5,1);

 %xlswrite('SSA_CEC.xlsx',PSO_date)
% MEANMEANMEAN(:,1)=mean(SAVE_PSO_score_save,2);
