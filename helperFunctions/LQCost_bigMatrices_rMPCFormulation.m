function [Bhat, Ahat ,H ,Ac ,bc ,Cc] = LQCost_bigMatrices_rMPCFormulation(Q,R,P,nx,nu,N,A,B,uLimA,uLimb,xLimA,xLimb,xLimAf,xLimbf,xLimAe,xLimbe)
% creates the big matrices for LQ-formulation of rMPC as defined in Robust model predictive control of constrained linear systems with bounded disturbances D.Q. Maynea,∗, M.M. Seronb, S.V. Rakovi ́ ca
% LQCost x_NPx_N + \sum_{k=0}^{N-1}x_kQx_k + u_kRu_k + x_N'Px_N
% where matrices Q,R,A,B all depend on a parameter p
% Constraints: 
%(1)       x_i = Ax_{i-1} + Bu_i, i = 1,..., N
%(2)      Au u_i<= bu i = 0,...,N-1
%(3)       Ax x_i <= bx i = 0,...,N-1
%(4)       Af x_N <= bf 
%(5)        A_{x0} x_0 <= b_{x0} + C_{x0} x; Comes from the condition x\in {x_0}\oplus E
% Assuming E = {e: A_e e<= b_e} matrices in (5) are given by
%(5')    A_{x0} = -A_e,   b_{x0} = b_e,   C_{x0} = -A_e

% End result is J(v) = vHv s.t. Ac v \leq bc + Cc x where v' = [u', x0'];




% matrices ----- such that the cost is given by [z, x]' [H G; G' W][z; x]!!!!!!!!!!!!!!!!!!
Bhat = sparse((N+1)*nx,N*nu,nu*nx*N);
Ahat = sparse((N+1)*nx,nx,nx*nx*N);
Ahat(1:nx,:) = eye(nx);
for it = 1: N
    Bhat(it*nx+1:it*nx+nx, : ) = A*Bhat(it*nx-nx+1:it*nx,:);
    Bhat(it*nx+1:it*nx+nx,it*nu-nu+1:it*nu) = B;
    Ahat(it*nx+1:it*nx+nx, : ) = A*Ahat(it*nx-nx+1:it*nx,:);
end

Qbar =sparse( [kron(eye(N),Q) zeros(N*nx,nx);zeros(nx,N*nx) P]);
Rbar = sparse(kron(eye(N),R));
bar_BQB_bar = Bhat'*Qbar*Bhat;
bar_AQA_bar = Ahat'*Qbar*Ahat;
bar_AQB_bar = Ahat'*Qbar*Bhat;
bar_BQA_bar = Bhat'*Qbar*Ahat;



H =[bar_BQB_bar+Rbar, bar_BQA_bar;
    bar_AQB_bar, bar_AQA_bar];
H=(H+H')/2;


%% constrains --  Ac*U<bc + Cc x
%constraints from the input
delCstr = abs(uLimb)== Inf;
uLimA(delCstr,:) = [];
uLimb(delCstr,:) = [];

%Constraints on initial state
Ax0 = -xLimAe;
bx0 = xLimbe;
Cx0 = -xLimAe;


%constraints on the input
Au = kron(eye(N),sparse(uLimA));
bu = kron(ones(N,1),sparse(uLimb));
Cu = sparse(zeros(size(bu,1),nx));


%constraints on states i = 0,...,N-1 and terminal set
xLimSetA =  kron(eye(N),sparse(xLimA))  ;
xLimSetA = [xLimSetA, zeros(size(xLimSetA,1),nx);
    zeros(size(xLimAf,1), size(xLimSetA,2)), xLimAf];
xLimSetb = [kron(ones(N,1),sparse(xLimb));xLimbf];

Au_x = xLimSetA*Bhat;
bu_x = xLimSetb;
Cu_x = -xLimSetA*Ahat;

Au = [Au;Au_x];
bu = [bu;bu_x];
Cu = [Cu;Cu_x];

Ac = [Au, -Cu;
    zeros(length(bx0),nu*N), Ax0];
bc = [bu;bx0];
Cc = [Cu*0;Cx0];



delCstr = abs(bc) == Inf;
Ac(delCstr,:) = [];
bc(delCstr,:) = [];
Cc(delCstr,:) = [];

end
