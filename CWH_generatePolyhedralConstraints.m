function [xLim, uLim] =  CWH_generatePolyhedralConstraints(m,xOff)
nu = 3;
nx = 6;


% input constraints
uLim.b = ones(nu*2,1)*1e1/m; %m/s^2
uLim.A = [-eye(nu); eye(nu)];

%state constraints

% min max values for state box constraints
xBox =  ([[-ones(3,1)*10000 ones(3,1)*10000]; ...
    [-ones(3,1)*3 ones(3,1)*3]])+xOff;    % adding some constraints to have a compact set
xBox(2,1) = 3+xOff(2);

%linear constraints derived from saturation values / box
xLim.A =kron([-1;1],eye(nx));
xLim.b= [-xBox(:,1); xBox(:,2)];

% generate Cone constraints
gamma = 15/180*pi;
N = 50; % number of inequality constraints in the cone approximation
normal = [0 1 0]';
v = [0 -1 0]'+xOff(1:3);
[Acone, bcone] = generateCone(gamma,N,normal,v);
Acone = [Acone zeros(size(Acone))];
xLim.A = [xLim.A ; Acone];
xLim.b = [xLim.b ; bcone];

end