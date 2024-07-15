clearvars
clc,close all


%% model parameters


% -- 1 -- System Parameters
mu = 398600.4; % gravitational constant for earth [km^3/s^2]
Re = 6378;     % Earth radius [km]
Ra = 500; % altitude [km]
m = 100;
r0 = Re+Ra;
n = sqrt(mu/r0^3);
Ts = 0.5;
TOrb  = 2* pi/ n; %orbital period, just to get an idea
tSim = 110; % simnulation time (start and stop) [s]

% dimensions of states, inputs, disturbance and reference.
nx = 6;
nu = 3;
nv = 3;
nw = 6;

figId = 0;
sysID = fullfile('results','20231918',strcat('rMPC_CWH_c',num2str(figId)));

% --2-- Linear and nonlinear models
% continuous linear system when considering the equilibrium theta =0

Ac = [zeros(3,3) eye(3);
    3*n^2 zeros(1,3) 2*n 0;
    zeros(1,3) -2*n 0 0;
    zeros(1,2) -n^2 zeros(1,3);
    ];

Bc = [zeros(3,3);eye(3)];
Bwc = [1*eye(3)*1e-1, zeros(3);zeros(3) eye(3)*1e-2];

% DTLTI

dtSys = c2d(ss(Ac,[Bc Bwc],[],[]),Ts);
Ad = dtSys.A;
Bwd = dtSys.B(:,nu+1:end);
Bd = dtSys.B(:,1:nu);

% -- 2.b-- References and equilibria

ic2r_func = @(ic) ic(1:3,1);
r2xSS_func = @(r) [r;zeros(3,size(r,2))];
r2uSS_func = @(r) [-3*n^2 0 0;zeros(1,3); 0 0 n^2]*r;

%desired end reference command -- in this inplementation We shift the problem by xRef
vRef = [0;4;0];
xRef = r2xSS_func(vRef);



% --3-- constraints
[xLim, uLim] =  CWH_generatePolyhedralConstraints(m,-xRef);

% --4-- ini Conditions
iniCon= [10;100;20;0;0;0]-xRef;

% --5-- Disturbance sets
dist.A = kron(eye(nx),[1;-1]);
dist.b = kron(Bwd*[ones(3,1)*1e-2;ones(3,1)],ones(2,1));

%uses the CORA toolbox for zonotopes. User manual: https://tumcps.github.io/CORA/data/archive/manual/Cora2024Manual.pdf
dist.intervalRep = interval(-Bwd*[ones(3,1)*1e-2;ones(3,1)],Bwd*[ones(3,1)*1e-2;ones(3,1)]);
dist.zonotope = dist.intervalRep.zonotope;

%% Controller Design
flg_rLQ = 0;

% --0-- Controller and IRG costs and parameters
MPC.Q = diag([10 10 10 1 1 1]);
MPC.R = eye(nu)*1e-4;
[LQR.K,MPC.P] = dlqr(Ad,Bd,MPC.Q,MPC.R);
LQR.K = -LQR.K;
MPC.N = 200;
MPC.maxIt = 2000;
MPC.eps = 1e-5;

% --1-- disturbance rejection feedback control law
LQR.agrQ = diag([10 10 10 1 1 1])*1e0;
LQR.agrR = eye(nu)*1e-4;
LQR.agrK = -dlqr(Ad,Bd,LQR.agrQ,LQR.agrR);

%--3-- error dynamics and error bounding sets
% error dynamics: e_+ = (A+BK)e + w
% error bounding sets: E_1 = A_cl E_0 + W; E_2 = A_cl^2 E_0 + A_cl W + W,... E_n = A_cl*E_{n-1} + W

Aerr = Ad+Bd*LQR.agrK;

% compute an epsilon-close outer approximation of the mRPI set, where epsilon =eps_mRPI below.
eps_mRPI = 1e-8;%5e-5;
[alphaVal, sVal] = paramsForEps_mRPISet(Aerr,dist.A,dist.b,eps_mRPI); %find the parameters for the approximation
alphaVal = max(0,alphaVal);

% compute {E_k}

vx0 = zeros(nx,1); % error uncertainty at time 0 (vertex rpz). We assume no initial uncertainty, i.e. the state is fully measured.
[Einf.seqzTOP,Einf.zTOP] = computeNStepReachableSetCORAToolBox(Aerr,eye(nx),Ts,alphaVal,sVal,dist.zonotope,vx0); %sequence of zonotopes. each element is a k-step reachable set. The last set is the outer approximation of the mRPI
KEinf.seqzTOP  = cellfun(@(y) mtimes(LQR.agrK,y), Einf.seqzTOP,'UniformOutput',false ); % set of possible feedback input values

%Constraint tightening

% x Constraints
seqPolytopeTemp = cellfun(@(y)minkDiff(polytope(xLim.A,xLim.b),y),Einf.seqzTOP,'UniformOutput',false );
seqPolyhedronTemp = cellfun(@(y) minHRep(Polyhedron(y.A,y.b)),seqPolytopeTemp,'UniformOutput',false ); % uses the toolbox MPT3
xLim.ABarSeq = cellfun(@(y)y.A,seqPolyhedronTemp,'UniformOutput',false );
xLim.bBarSeq = cellfun(@(y)y.b,seqPolyhedronTemp,'UniformOutput',false );
% u Constraints
seqPolytopeTemp = cellfun(@(y)minkDiff(polytope(uLim.A,uLim.b),y),KEinf.seqzTOP,'UniformOutput',false );
seqPolyhedronTemp = cellfun(@(y) minHRep(Polyhedron(y.A,y.b)),seqPolytopeTemp,'UniformOutput',false );
uLim.ABarSeq = cellfun(@(y)y.A,seqPolyhedronTemp,'UniformOutput',false );
uLim.bBarSeq = cellfun(@(y)y.b,seqPolyhedronTemp,'UniformOutput',false );

% we just care about the last element on the sequence in this one
xLim.ABar = xLim.ABarSeq{end};
xLim.bBar = xLim.bBarSeq{end};
uLim.ABar = uLim.ABarSeq{end};
uLim.bBar = uLim.bBarSeq{end};


%-- 4-- Terminal set - based on MOAS for a fixed reference (the origin of the shifted system)
kfin = 5000;
ABigLim = [xLim.ABar;uLim.ABar*LQR.K];
bBigLim = [xLim.bBar;uLim.bBar];
[xLim.AOi,xLim.bOi,xLim.tstar] = computeBasicMoas(Ad+Bd*LQR.K,diag(nx),ABigLim,bBigLim, kfin);

% --5-- build the MPC matrices with the constrained set

% if you don't want zonotopes you probably need to compute an H-representation of Einf.zTOP, should be
% given in the CORA manual. I'm guessing it'll be : Einf.H = polytope(Einf.zTOP). Or somerthing like that
% the file helperFunctions>LQCost_bigMatrices_rMPCFormulation computes the "big matrices" assuming an H representation of E0
% [MPC.Bhat, MPC.Ahat ,MPC.H,MPC.F ,MPC.Ac ,MPC.bc ,MPC.Cc,MPC.nLambda] =... 
% LQCost_bigMatrices_rMPCFormulation_zonotopex0(MPC.Q,MPC.R,MPC.P,nx,nu,MPC.N,Ad,Bd,uLim.ABar,uLim.bBar,xLim.ABar,xLim.bBar,xLim.AOi,xLim.bOi,Einf.zTOP);

for i = 1:MPC.N
    disp(i+" / "+MPC.N)
    iter = MPC.N - i + 1;
    [MPC.Bhat.el{i}, MPC.Ahat.el{i} ,MPC.H.el{i},MPC.F.el{i} ,MPC.Ac.el{i} ,MPC.bc.el{i} ,MPC.Cc.el{i},MPC.nLambda.el{i}] =... 
    LQCost_bigMatrices_rMPCFormulation_zonotopex0(MPC.Q,MPC.R,MPC.P,nx,nu,iter,Ad,Bd,uLim.ABar,uLim.bBar,xLim.ABar,xLim.bBar,xLim.AOi,xLim.bOi,Einf.zTOP);
    %LQCost_bigMatrices_rMPCFormulation_zonotopex0(MPC.Q,MPC.R,MPC.P,nx,nu,iter,Ad,Bd,uLim.ABar,uLim.ABar,xLim.ABar,xLim.bBar,xLim.AOi,xLim.bOi,Einf.zTOP);


end


% Put all the qp matrices inside the MPC.qp structure
MPC.qp.A = MPC.Ac;
MPC.qp.b = MPC.bc;
MPC.qp.H = MPC.H;
% the linear term in the objective function is defined inside the loop:
% MPC.qp.f = obj.Mint*x0;

opsYAL = sdpsettings('solver','mosek','verbose',0,'debug',0,'cachesolvers',1);
% uncomment following line to run quadprog but it might take a lot of time / be laggy
% opsYAL = sdpsettings('solver','quadprog','verbose',0,'debug',0,'cachesolvers',1);
%xYAL = sdpvar(nu*MPC.N + MPC.nLambda,1); %if not using zonotopes, MPC.nLambda should be replaced by nx


%% running the simulation
clear DATA
flg_svPlot = 1;
%NSamps = ceil(tSim/Ts);

NSamps = MPC.N;
tVec = (0:NSamps-1)*Ts;
commonInfo.uLim = uLim;
commonInfo.xLim = xLim;
commonInfo.dist = dist;
commonInfo.Einf = Einf;
commonInfo.KEinf = KEinf;
commonInfo.sVal = sVal;
commonInfo.alphaVal = alphaVal;
commonInfo.solverInfo = opsYAL;
commonInfo.MPC = MPC;
commonInfo.LQR = LQR;

DATA.x = zeros(nx,NSamps);
DATA.u = zeros(nu,NSamps);
DATA.w = zeros(nx,NSamps);
DATA.xNom = zeros(nx,NSamps);
DATA.uNom = zeros(nu,NSamps);
DATA.tComp = zeros(1,NSamps);

% generate the disturbance (I fixed how it is generated so I can better compare)
gNum = 32712938;
gType = 'twister';
rng(gNum,gType);
DATA.w  = dist.b([1:nx/2 3*nx/2+1:2*nx],1).*(rand(nx,NSamps)-0.5)*2;


DATA.x(:,1) = iniCon;
DATA.xNom(:,1) = DATA.x(:,1);
xSS_t = r2xSS_func(vRef);
uSS_t = r2uSS_func(vRef);


L = 2; %Length Interval
Intervals = [0:L:MPC.N]; %Define the sequence of intervals
Np = (ceil(MPC.N)/L); %initial number of paramenters
Np0 = Np;

H = zeros(MPC.N,Np);
for i = 0:(Np-1)
H((L*i)+1:L*(i+1),i+1) = ones(L,1);
end

nVar = zeros(NSamps,1);

start = tic;
for it = 1:NSamps
    disp(it+" / "+NSamps)
    
    
    % --1-- solve uMPC OPC for current reference
    MPC.qp.b = MPC.bc.el{it}  +  MPC.Cc.el{it}*(DATA.x(:,it)); % the Cc matrix will be different iin the formulation without zonotopes
    MPC.qp.f = MPC.F.el{it}'*DATA.x(:,it);
    
    Np = ceil((MPC.N-it+1)/L); %initial number of paramenters
    
    xYAL = sdpvar(nu*Np + MPC.nLambda.el{1},1); %if not using zonotopes, MPC.nLambda should be replaced by nx
    nVar(it) = nu*Np;
    %xYAL = sdpvar(nu*(MPC.N-it+1) + MPC.nLambda.el{1},1); %if not using zonotopes, MPC.nLambda should be replaced by nx
    t_xYAL = xYAL(1:nu*Np);
    tt_xYAL = reshape(t_xYAL,nu,Np);
    
    H_size = size(H,2);
    ut = H(1:end,(H_size-Np+1):end)*tt_xYAL';
    
    ut = ut';
    opt = [reshape(ut,(MPC.N-it+1)*nu,1);xYAL(nu*Np+1:end)]; 
    
    % define constraints for Yalmip
    CLax = MPC.qp.A.el{it}*opt<=MPC.qp.b;
    % define objective function for Yalmip
    OBJ = .5*opt'*MPC.qp.H.el{it}*opt + 0.5*MPC.qp.f'*opt;
    tic
    outInfo = optimize(CLax,OBJ,opsYAL);
    if outInfo.problem
        prob = 1
    end
    uOut = value(opt);
    
    J_val(it) = .5*uOut'*MPC.qp.H.el{it}*uOut;

    %separate x0 and input sequence
    DATA.xNom(:,it) =  DATA.x(:,it) + Einf.zTOP.G*uOut(end-MPC.nLambda.el{it}+1:end,1);
    uOpt = reshape(uOut(1:end-MPC.nLambda.el{it}),nu,[]);

    % --3-- aCompute current control input
    uNom = uOpt(:,1);
    uCur = uNom + LQR.agrK*(DATA.x(:,it) - DATA.xNom(:,it));
    %uCur = max(min(uCur,[0.1;0.1;0.1]),-[0.1;0.1;0.1]);
    DATA.tComp(1,it) = toc;
    DATA.tComp(2,it) = DATA.tComp(1,it)-outInfo.yalmiptime;
    DATA.x(:,it+1) = Ad*DATA.x(:,it) + Bd*uCur + DATA.w(:,it);
    
    % save to buffers
    DATA.u(:,it) = uCur;
    DATA.uNom(:,it) = uNom;
    
    H = H(2:end,:);
end
stop = toc(start)

DATA.x(:,end) = [];
DATA.xNom(:,end) = [];
% go back to original coordinates
DATA.x = DATA.x + xRef;
DATA.xNom = DATA.xNom + xRef;
extraInfo.rdnGen.type = gType;
extraInfo.rdnGen.num = gNum;

% save simulation data
DATA_w_Param = DATA; stop_w_P = stop; J_w_P = J_val; tVec_P = tVec; nVar_P = nVar;
save('Results/res_w_P_L2','DATA_w_Param','stop_w_P','tVec_P','J_w_P','nVar');

% %% basic plot
% 
% 
% figure;
% subplot(2,1,1);plot(tVec, DATA.x(1:3,:),'linewidth',2.5);grid on;
% subplot(2,1,2);plot(tVec, DATA.x(4:6,:),'linewidth',2.5);grid on;
% 
% figure;
% subplot(2,1,1);plot(tVec, DATA.u,'linewidth',2.5);grid on;
% 
% figure;
% subplot(2,1,1);plot(tVec, DATA.w(1:3,:),'linewidth',2.5);grid on;
% subplot(2,1,2);plot(tVec, DATA.w(4:6,:),'linewidth',2.5);grid on;




%% some functions

% for error bounding set
function [alphaVal,sVal] = paramsForEps_mRPISet(Ak,A_pol,b_pol,epsVal)
% given \epsilon determine value of s such that (1-\alpha)^{-1} E_s is an \epsilon outer approximation of E_{\infty}.
% See: Algoirthm 1 secIV of Rakovic, Sasa V., et al. "Invariant approximations of the minimal robust positively invariant set." IEEE Transactions on automatic control 50.3 (2005): 406-410.
flg_keepGoing = 1;
nx = size(Ak,1);
sCur = 0;
nCstr = length(b_pol);
x = sdpvar(nx,1);
CLax = A_pol*x<=b_pol;
yalmOpts = sdpsettings('solver','CDD','verbose',0,'cachesolvers',1);
hwFunc = @(vec) optimize(CLax,-(vec')*x,yalmOpts); % support function of the set
gi = b_pol';
fi = A_pol';
eBase = eye(nx);

for it = 1:nx
    hwFunc(eBase(:,it));
    hwAsT_ejp(it,1)=value(x)'*eBase(:,it);
    hwFunc(-eBase(:,it));
    hwAsT_ejm(it,1)=-value(x)'*eBase(:,it);
end

while flg_keepGoing
    sCur = sCur+1;
    As = Ak^sCur;
    AsT_fi = (As')*fi;
    % compute \alphs^o(s)
    for it = 1:nCstr
        hwFunc(AsT_fi(:,it));
        hwAsT_fi = value(x)'*AsT_fi(:,it);
        alphaCand(it) = hwAsT_fi/gi(it);
    end
    alphaCur = max(alphaCand(it));
    %compute M(s)
    for it = 1:nx
        hwFunc(As'*eBase(:,it));
        hwAsT_ejp(it,sCur+1)=value(x)'*(As'*eBase(:,it));
        hwFunc(-As'*eBase(:,it));
        hwAsT_ejm(it,sCur+1)=-value(x)'*(As'*eBase(:,it));
    end
    Ms = max(max(sum(hwAsT_ejp,2),sum(hwAsT_ejm,2)));
    flg_keepGoing = alphaCur>(epsVal/(epsVal+Ms));

end
sVal = sCur;
alphaVal = alphaCur;
end

function [seqOfReachableSets,mRPI_outerApprox] = computeNStepReachableSetCORAToolBox(ADyn,BDyn,Ts,alphaVal,Nval,uZonotope,vx0);

DTSys = linearSysDT(ADyn,BDyn,Ts);
params.tStart = 0;
params.tFinal= Nval*Ts;
params.R0 = zonotope(interval(vx0,vx0));
params.U = uZonotope;
options.zonotopeOrder = 5;
options.reductionTechnique = 'girard';
rop  = reach(DTSys,params,options);
seqOfReachableSets = rop.timePoint.set;
mRPI_outerApprox = (rop.timePoint.set{end});
mRPI_outerApprox = mtimes(1/(1-alphaVal),mRPI_outerApprox);
seqOfReachableSets{end+1} = mRPI_outerApprox;
seqOfReachableSets(1) = [];
end
