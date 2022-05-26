function [Global_score,Global_pos,Convergence_curve] = IJADE(SearchAgents_no,FEsmax,lb,ub,dim,fobj,fnum)


%% 初始化参数
lb=lb.*ones(dim,1);ub=ub.*ones(dim,1);
NP = SearchAgents_no;
 uF = 0.8;
 uCR = 0.5;
Ar = 2.6;
p = 0.11;
NPmin = 50;
NPinit = NP;

	A		= [];
    fA=[];
	nA		= [];
D = dim;
counteval = 0;
countiter = 1;
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* rand(D, 1);
        fx(i) = feval(fobj,X(:,i),fnum);
        	counteval = counteval + 1;
            Convergence_curve(counteval)=fx(i);
    end
       curve(countiter)=min(fx);
    countiter = countiter + 1;
% Initialize archive
% if isempty(A)
 	Asize = round(Ar * NP);

if isempty(nA)
	nA = 0;
else
	nA = min(nA, Asize);
end
[fx, fidx] = sort(fx);
X = X(:, fidx);




% Initialize variables
V = X;
U = X;
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_df = zeros(1, NP);	% Set of df
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;

%%   进入循环
while counteval < FEsmax
	% Termination conditions	
	% Memory Indices
	% Crossover rates
	CR = uCR + 0.1 * randn(1, NP);	
	CR((CR < 0) ) = uCR + 0.1 * randn(1, numel(find(CR<0)));
	CR(CR > 1) = uCR + 0.1 * randn(1, numel(find(CR>1)));	
    %% 对CR进行排序
    CR=sort(CR);
	% Scaling factors
    %% 产生F
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = uF + Chy(iChy);
			iChy = mod(iChy, numel(Chy)) + 1;
		end
	end
	F(F > 1) = 1;	
	% pbest
   pbest = 1 + floor(max(2, round(p * NP)) * rand(1, NP));	

 	XA = [X, A];

	% Mutation 变异环节
	for i = 1 : NP		
		% Generate r1

		r1 = floor(1 + NP * rand);
		while i == r1
			r1 = floor(1 + NP * rand);
        end	
		% Generate r2
		r2 = floor(1 + (NP + nA) * rand);
		while i == r1 || r1 == r2
			r2 = floor(1 + (NP + nA) * rand);
        end

			V(:, i) = X(:, i) ...
				+ F(i) .* (X(:, pbest(i)) - X(:, i)) ...
				+ F(i) .* (X(:, r1) - XA(:, r2));


    end
		% Correction for outside of boundaries 边界处理
		for i = 1 : NP
			for j = 1 : D
				if V(j, i) < lb(j)
						V(j, i) = 0.5 * (lb(j) + X(j, i));

				elseif V(j, i) > ub(j)
						V(j, i) = 0.5 * (ub(j) + X(j, i));

				end
			end
        end
        %% 交叉环节
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		
			% Binomial Crossover
            m=0;
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
                    m=m+1;
				else
					U(j, i) = X(j, i);
				end
			end

	end
	
	% Evaluation 评价环节
    	nS = 0;
  	for i = 1 : NP  
	 fu(i)=feval(fobj,U(:,i),fnum);
    % Selection
		if fu(i) < fx(i)
			nS			= nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			S_df(nS)	= abs(fu(i) - fx(i));		
			if nA < Asize
				A(:, nA + 1)	= X(:, i);
                fA(nA + 1)=fx(i);
				nA				= nA + 1;
			else
				ri				= floor(1 + Asize * rand);
				A(:, ri)		= X(:, i);
                fA(ri)=fx(i);
            end		
               X(:, i)		= U(:, i);
               fx(i)		= fu(i);	
        elseif fu(i) == fx(i)			
			X(:, i)		= U(:, i);

        end
            counteval = counteval + 1;
            Convergence_curve(counteval)=min(fx);
    end	
	% Update uCR and uF 更新
	if nS > 0
        c=0.1;
            newSf=(sum(S_F(1:nS).^2))/(sum(S_F(1:nS)));
            uCR=(1-c)*uCR+c.*mean(S_CR(1:nS));
            uF=(1-c)*uF+c.*newSf;
    end	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
    [fA, fidA] = sort(fA);
	A = A(:, fidA);
	% Update NP and population 更新种群数
     NP = round(NPinit - (NPinit - NPmin) * counteval / FEsmax);
    fx = fx(1 : NP);
    X = X(:, 1 : NP);
	U = U(:, 1 : NP);
    Asize = round(Ar * NP);	
	if nA > Asize
		nA = Asize;
		A = A(:, 1 : Asize);
        fA=fA(1 : Asize);
	end
       curve(countiter)=fx(1);
    countiter = countiter + 1;

end

% Global_score,Global_pos,Convergence_curve
Global_score = fx(1);
Global_pos = X(:, 1);
    Convergence_curve=Convergence_curve(500:500:FEsmax);

	
