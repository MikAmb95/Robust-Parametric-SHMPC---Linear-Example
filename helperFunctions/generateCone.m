% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate a polyhedral approximation for a rotated and translated cone 
%  aligned with n-vector and vertex at v of x-y-z coordinate frame
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [Acone,bcone] = generateCone(gamma,N,normal,v)

if norm([0,1,0]'-normal)<1e-6,
    [Acone, bcone] = generateBasicCone(gamma,N);
else,
    axisOfRotation = cross([0,1,0],normal);
    axisOfRotation = axisOfRotation/norm(axisOfRotation);
    
    angleOfRotation = acos(dot(normal,[0,1,0])/norm(normal));
    
    e = axisOfRotation;
    S = [0, -e(3), e(2); e(3), 0, -e(1); -e(2), e(1), 0];
    R = expm(-S*angleOfRotation);
    
    [Acone, bcone] = generateBasicCone(gamma,N);
    Acone = Acone*R;
end;

bcone = bcone + Acone*v(:);


return

% 
% for m = 1:1:100000,
%     x=(rand([3,1])-0.5)*10;
%     if Acone*x<=bcone, plot3(x(1),x(2),x(3),'r'); hold on; end;
% end;