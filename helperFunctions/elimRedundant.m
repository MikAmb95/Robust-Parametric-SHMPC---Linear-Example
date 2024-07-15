%   function [A1,b1]=elimRedundant(A,b,toll)
%   Eliminate all redundant constraints in Ax<=b
%   The output set is A1<=b1
%  (c) Ilya V. Kolmanovsky


function [A1,b1,IDX]=elimRedundant(A,b,toll)


if nargin==2
         toll=1e-6;
end

A1=A;
b1=b;
sizeA=size(A);

n=sizeA(1);
options=optimset('Display','off');

% for it=1:1:n,
%     % disp(i);
%     f=-A1(1,:);
%     Aa=A1(2:n,:);
%     ba=b1(2:n);
%     hx=lp(f,Aa,ba,[],[],[],[],-1);
%     % hx=linprog(f,Aa,ba,[],[],[],[],[],options);
%     hxval=-f*hx;
%     % disp(hxval-b(i))
%     if hxval<b(it)+toll
%         A1=Aa;
%         b1=ba;n=n-1; % disp('Constraint Eliminated')
%     else
%         A1=[Aa;-f];
%         b1=[ba;b1(1)];
%     end
%     if n==1
%         break; disp('Warning: only one active constraint remains');
%     end
% end

  IDX = [2:n]; % active constraints
  n_del = 0;
%   tic1 = [];
%   tic2 = [];
  
for it=1:1:n
    % disp(i);   
%     tic
    hx=lp(-A(it,:),A(IDX,:),b(IDX),[],[],[],[],-1);
%     tic1 = [tic1 toc];
%     tic
%     hx=linprog(-A(it,:),A(IDX,:),b(IDX),[],[],[],[],[],options);
    hxval=A(it,:)*hx;
    % disp(hxval-b(i))
    IDX(it-n_del) =  it;
    if hxval<b(it)+toll
        IDX(it-n_del) = [];
        n_del = n_del+1;% disp('Constraint Eliminated')
    end
%     tic2 = [tic2 toc];
    if n_del+1 == n
        IDX = n;
        break; disp('Warning: only one active constraint remains');
    end
end

A1 = A(IDX,:);
b1 = b(IDX);

return

% TEST EXAMPLE 1
%A=[1,0;-1,0;0,1;0,-1;1,0;-1,0;0,1;0,-1];
%b=[1;1;1;1;1;1;1-1e-16;1-1e-16];
%[A1,b1]=elimRedundant(A,b,0)