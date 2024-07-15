% O-infinity for autonomous AS stable:
% x^+ = Ax 
% subject to C*x \in Y
% where Y = {y: A_y y\leq b_y}
% kfin max iteration for which we compute O_{kfin}
function [AOi,bOi,tstar] = computeBasicMoas(A,C,Ay, by,kfin)

toll    = 1e-5; % tolerance parameter
epsVal = 0.00001;% epsBall taken outin the definition of \tilde Y_n: \tilde Y_n\oplu\epsilon \mathcal B \subseteq Y_{cl,n,\infty} where instead of \infty we take n+kfin
nCst = size(Ay,1);

byFin = by-epsVal;

% -------- Define some Matrices to ease the prediction of the state
A_pow_k = cell(kfin+1,1); A_pow_k{1} = eye(length(A));
for k = 2:kfin+1
    A_pow_k{k} = A*A_pow_k{k-1};
end


% --- Constraints on current state, i.e. k = 0

% build constraints on current state for O_{\infty,0}
AOi = Ay*C;
bOi = byFin;

% --- Constraints on predicted state at time k=1,2,..,kfin
for k = 1:1:kfin
% k
    AOi_pre = AOi;
    bOi_pre = bOi;
    % build constraints on current state for O_{\infty,0}
    AOi= [AOi; Ay*C*A_pow_k{k+1}];
    bOi =  [bOi;by];

    % delete constraints from time to time + check if we have finalized computation of the O_infty
    [AOi, bOi] = elimRedundant(AOi, bOi, toll);
    if ~isempty(AOi_pre),
        if issubset(AOi_pre, bOi_pre, AOi, bOi, toll)=='y'
            break;
        end
    end

end
tstar=k;
%check if we have completed the computation of the MOAS
if ~isempty(AOi_pre) && k ==kfin
    if issubset(AOi_pre, bOi_pre, AOi, bOi, toll)=='y',
        tstar = tstar-1;
    end
end
% [ min(bOi),tstar]
end