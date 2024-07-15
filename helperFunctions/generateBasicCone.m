% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate a polyhedral approximation for a base cone aligned with y-axis 
% and vertex at the origin of x-y-z coordinate frame
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [Acone,bcone] = generateBasicCone(gamma,N)

% ++ generate a list of vectors
x = [];
y = [];
z = [];
n = [];
Acone = [];
bcone = [];
for (i=1:1:N),
    theta(i)= 2*pi/(N - 1)*(i - 1);
    x(i)    = tan(gamma)*cos(theta(i));
    z(i)    = tan(gamma)*sin(theta(i));
    y(i)    = 1;
    
    if i>1,
        n(:,i-1)=cross([x(i - 1),y(i - 1),z(i - 1)],[x(i),y(i),z(i)])';
        n(:,i-1)=n(:,i - 1)/norm(n(:,i - 1));
        
        Acone=[Acone; n(:,i - 1)'];
        bcone=[bcone; 0];
        
    end;
    
end



return

% 
% for m = 1:1:100000,
%     x=(rand([3,1])-0.5)*10;
%     if Acone*x<=bcone, plot3(x(1),x(2),x(3)); hold on; end;
% end;