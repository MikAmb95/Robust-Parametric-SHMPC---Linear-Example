function [Bhat, Ahat ,H, F ,Ac ,bc ,Cc,nLambda] = LQCost_bigMatrices_rMPCFormulation_zonotopex0(Q,R,P,nx,nu,N,A,B,uLimA,uLimb,xLimA,xLimb,xLimAf,xLimbf,x0Zonotope)
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
Huu = Bhat'*Qbar*Bhat +Rbar;
Hxx = Ahat'*Qbar*Ahat;
Hxu = Ahat'*Qbar*Bhat;
Hux = Bhat'*Qbar*Ahat;


%H =[Huu, Hux;
%   Hxu, Hxx];

% add zonotope change of var: x  = x_0 + G\lambda cost becomes: [z \lambda] H [z;\lambda] + x_0F [z;\lambda] + x_0 J x_0 <-- we don't worry about this guy 
GzTOP = x0Zonotope.G;
nLambda = size(GzTOP,2);

Huu = Huu;
Hlambdau = (GzTOP')*Hxu;
Hulambda = Hux*GzTOP;
Hlambdalambda = (GzTOP')*Hxx*GzTOP + .001*eye(nLambda);

H = [Huu Hulambda;
    Hlambdau Hlambdalambda];
Fu = 2*Hxu;
Flambda = 2*Hxx*GzTOP;

F = [Fu Flambda];

% J = Hxx; we don't really care about this lill guy

H=(H+H')/2;




%% constrains --  Ac*U<bc + Cc x
%constraints from the input
delCstr = abs(uLimb)== Inf;
uLimA(delCstr,:) = [];
uLimb(delCstr,:) = [];


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

Ac = sparse([Au -Cu*GzTOP;
    zeros(2*nLambda,N*nu),kron([-1;1],eye(nLambda))]);
bc = sparse([bu;ones(2*nLambda,1)]);
Cc = sparse([Cu;zeros(2*nLambda,nx)]);

delCstr = abs(bc) == Inf;
Ac(delCstr,:) = [];
bc(delCstr,:) = [];
Cc(delCstr,:) = [];



end
