function new = mapping_parameters(original,Params_LowerBound,...
Params_UpperBound, parameterNames)
% To map the parameters inside the parameter boundary
%
% Inputs:
%   original: m*n matrix - n is the number of parameters, each column is a 
%       parameter vector
%   Params_LowerBound: a structure contains parameter names and lower
%       boudanry values. e.g.Params_LowerBound.x1=1; or
%       Params_LowerBound.x1=-inf;
%   Params_UpperBound: a structure contains parameter names and Upper
%       boudanry values. e.g.Params_UpperBound.x1=10; or
%       Params_UpperBound.x1=inf;
%   parameterNames: cell format contains the parameter names. 
%       e.g.parameterNames = {'x1'};
% Output:
%   new: mapped parameter matrix (m*n)
%
% Example:
%   original=[-10:1:16]';
%   LowerBound.x1 = 3;
%   UpperBound.x1 = 6;
%   ParaName = {'x1'};
%   new = mapping_parameters(original,LowerBound,UpperBound,ParaName);
%
%
% QJ Wang & Jie Jian, the University of Melbourne, March 2019

nparams = length(parameterNames);

% check the format of original parameter matrix
if size(original,2)~=nparams
    error('Each parameter vector should be in a column.')
end

new = original;

% do for each parameter
for i=1:nparams
    LowerBound = Params_LowerBound.(parameterNames{i});% get the lower bound
    UpperBound = Params_UpperBound.(parameterNames{i});% get the upper bound
    
    if isinf(LowerBound)% lower boundary is -inf
        new(:,i) = UpperBound-abs( original(:,i)-UpperBound );
    elseif isinf(UpperBound)% upper boundary is +inf
        new(:,i) = LowerBound+abs( original(:,i)-LowerBound );
    else
    % Step 1: flip values smaller than lower boundary
    new(original(:,i)<LowerBound,i)=...
        2*LowerBound-original(original(:,i)<LowerBound,i);
    % Step 2: create periodic pattern
    new(:,i)=new(:,i)-2.*(UpperBound-LowerBound).*...
        floor( (new(:,i)-LowerBound)./(2.*(UpperBound-LowerBound)) );
    % Step 3: flip values larger than upper boundary
    new(new(:,i)>UpperBound,i)=2.*UpperBound-new(new(:,i)>UpperBound,i);
    end
end
end