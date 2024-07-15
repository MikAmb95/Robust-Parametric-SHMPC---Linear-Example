function [x,lambda,how]=lp(f,A,B,vlb,vub,x,neqcstr,verbosity)
%LP     Linear programming.                   
%       X=LP(f,A,b) solves the linear programming problem:
%        
%            min f'x    subject to:   Ax <= b 
%             x
%   
%       [x,LAMBDA]=LP(f,A,b) returns the set of Lagrangian multipliers,
%       LAMBDA, at the solution. 
%
%       X=LP(f,A,b,VLB,VUB) defines a set of lower and upper
%       bounds on the design variables, X, so that the solution is always in
%       the range VLB < X < VUB.
%
%       X=LP(f,A,b,VLB,VUB,X0) sets the initial starting point to X0.
%
%       X=LP(f,A,b,VLB,VUB,X0,N) indicates that the first N constraints defined
%       by A and b are equality constraints.
%
%       LP produces warning messages when the solution is either unbounded
%       or infeasible. Warning messages can be turned off with the calling
%       syntax: X=LP(f,A,b,VLB,VUB,X0,N,-1).

%       Copyright (c) 1990 by the MathWorks, Inc.
%       Andy Grace 7-9-90.
         
if nargin<8, verbosity = 0; 
        if nargin<7, neqcstr=0; 
                if nargin<6, x=[]; 
                        if nargin<5, vub=[];
                                if nargin<4, vlb=[];
end, end, end, end, end
[x,lambda,how]=qp([],f(:),A,B(:),vlb, vub, x(:),neqcstr,verbosity);                                
          
