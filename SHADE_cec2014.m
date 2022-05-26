function [Global_score,Global_pos,Convergence_curve] = SHADE_cec2014(SearchAgents_no,FEsmax,lb,ub,dim,fobj,fnum)
% LSHADE_SPS_CEC15 L-SHADE algorithm with SPS framework (variant CEC15)
% LSHADE_SPS_CEC15(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% LSHADE_SPS_CEC15(..., options) minimize the function by solver options.
% if nargin <= 4
% 	options = [];
% end
rand('seed', sum(100 * clock));
lb=lb.*ones(dim,1);ub=ub.*ones(dim,1);
NP = SearchAgents_no;
 F = 0.5;
 CR = 0.5;
Ar = 2.6;
p = 0.11;
H = 6;
NPmin = 4;
NPinit = NP;
	A		= [];
	nA		= [];
D = dim;
counteval = 0;
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* rand(D, 1);
        fx(i) = feval(fobj,X(:,i),fnum);
                	counteval = counteval + 1;
            Convergence_curve(counteval)=fx(i);
    end
% Initialize archive
if isempty(A)
	Asize = round(Ar * NP);
	A = zeros(D, Asize);
	for i = 1 : Asize
		A(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
else
	[~, Asize] = size(A);
	if Asize > round(Ar * NP)
		Asize = round(Ar * NP);
		A = A(:, 1 : Asize);
	elseif Asize < round(Ar * NP)
		Asize = round(Ar * NP);
		A = zeros(D, Asize);
		for i = 1 : Asize
			A(:, i) = lb + (ub - lb) .* rand(D, 1);
		end
	end
end

if isempty(nA)
	nA = 0;
else
	nA = min(nA, Asize);
end
[fx, fidx] = sort(fx);
X = X(:, fidx);
% MF
	MF = F * ones(H, 1);
% MCR
	MCR = CR * ones(H, 1);
% iM
	iM = 1;



% Initialize variables
V = X;
U = X;
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_df = zeros(1, NP);	% Set of df
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;


while counteval < FEsmax
	% Termination conditions	
	% Memory Indices
	r = floor(1 + H * rand(1, NP));	
	% Crossover rates
	CR = MCR(r)' + 0.1 * randn(1, NP);	
	CR((CR < 0) | (MCR(r)' == -1)) = 0;
	CR(CR > 1) = 1;	
	% Scaling factors
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = MF(r(i)) + Chy(iChy);
			iChy = mod(iChy, numel(Chy)) + 1;
		end
	end
	F(F > 1) = 1;	
	% pbest
	pbest = 1 + floor(max(2, round(p * NP)) * rand(1, NP));	
	% Population + archive
	XA = [X, A];
	% Mutation
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
		% Correction for outside of boundaries
		for i = 1 : NP
			for j = 1 : D
				if V(j, i) < lb(j)
						V(j, i) = 0.5 * (lb(j) + X(j, i));

				elseif V(j, i) > ub(j)
						V(j, i) = 0.5 * (ub(j) + X(j, i));

				end
			end
        end
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		
			% Binomial Crossover
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end


	end
	
	% Evaluation
    	nS = 0;
  	for i = 1 : NP  
	 fu(i)=feval(fobj,U(:,i),fnum);
    % Selection

		if fu(i) < fx(i)
			nS			= nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			S_df(nS)	= abs(fu(i) - fx(i));
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);			
			if nA < Asize
				A(:, nA + 1)	= X(:, i);
				nA				= nA + 1;
			else
				ri				= floor(1 + Asize * rand);
				A(:, ri)		= X(:, i);
			end						
        elseif fu(i) == fx(i)			
			X(:, i)		= U(:, i);

        end
                    counteval = counteval + 1;
            Convergence_curve(counteval)=min(fx);
    end	
	% Update MCR and MF
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		
		if all(S_CR(1 : nS) == 0)
			MCR(iM) = -1;
		elseif MCR(iM) ~= -1
			MCR(iM) = sum(w .* S_CR(1 : nS) .* S_CR(1 : nS)) / sum(w .* S_CR(1 : nS));
		end
		
		MF(iM) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
		iM = mod(iM, H) + 1;
    end	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	
	% Update NP and population

    fx = fx(1 : NP);
    X = X(:, 1 : NP);
	U = U(:, 1 : NP);
    Asize = round(Ar * NP);	
	if nA > Asize
		nA = Asize;
		A = A(:, 1 : Asize);
	end
   
    
end
% Global_score,Global_pos,Convergence_curve
Global_score = fx(1);
Global_pos = X(:, 1);
    Convergence_curve=Convergence_curve(500:500:FEsmax);

	
