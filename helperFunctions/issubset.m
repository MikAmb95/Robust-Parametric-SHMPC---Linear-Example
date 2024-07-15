% function issub=issubset(A,b,A1,b1,toll)
% Checks whether the set of x s.t. A1 x <=b1 is a subset
% of the set of x s.t. A x <= b
% issub takes values 'y' for yes or 'n' for no
% The parameter toll specifies the tolerance in determining set
% inclusion; 
% (c) Ilya V. Kolmanovsky

function issub=issubset(A,b,A1,b1,toll)

if nargin<5
    toll=1e-8;
end

if nargin<4
    disp('Error: Not enough input arguments specified'), return
end

szA1=size(A1);
m=szA1(1); % number of inequalitites which define A1

issub='y'; %

for it=1:1:m
    f=-A1(it,:);   %select the ith constraint, i=1,...,m, of A1, b1
    hx = lp(f,A,b,[],[],[],[],-1); %solve an lp problem to determine if redundant
    % hx=lp(f,A,b);
    % hx = linprog(f,A,b);
    hv = -f*hx-b1(it);
    %keyboard
    if hv>=toll, 
        issub='n'; break,
    end % if not redundant stop

end

return



%TEST EXAMPLE 1 (anser yes)

A=[1,0;-1,0;0,1;0,-1];
b=[1;1;1;1];

A1=[1,0;-1,0];
b1=[1;1];
issub=issubset(A,b,A1,b1)




%TEST EXAMPLE 2 (answer no)


A=[1,0;-1,0;0,1;0,-1];
b=[1;1;1;1];

A1=[1,0;-6,0];
b1=[1;1];
issub=issubset(A,b,A1,b1)

