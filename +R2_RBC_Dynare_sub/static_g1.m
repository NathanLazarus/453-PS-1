function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = R2_RBC_Dynare_sub.static_g1_tt(T, y, x, params);
end
g1 = zeros(4, 4);
g1(1,1)=(-exp((-y(1))))-(1+exp(y(3))*params(4)*exp((params(4)-1)*y(4))*exp((1-params(4))*y(2))-params(3))*params(2)*(-exp((-y(1))));
g1(1,2)=(-(exp((-y(1)))*params(2)*exp(y(3))*params(4)*exp((params(4)-1)*y(4))*(1-params(4))*exp((1-params(4))*y(2))));
g1(1,3)=(-(exp((-y(1)))*params(2)*exp(y(3))*params(4)*exp((params(4)-1)*y(4))*exp((1-params(4))*y(2))));
g1(1,4)=(-(exp((-y(1)))*params(2)*exp((1-params(4))*y(2))*exp(y(3))*params(4)*(params(4)-1)*exp((params(4)-1)*y(4))));
g1(2,1)=exp(y(3))*(1-params(4))*exp(params(4)*y(4))*(-exp((-y(1))));
g1(2,2)=(-((params(4)+params(5))*exp(y(2)*(params(4)+params(5)))));
g1(2,3)=exp((-y(1)))*exp(y(3))*(1-params(4))*exp(params(4)*y(4));
g1(2,4)=exp((-y(1)))*exp(y(3))*(1-params(4))*params(4)*exp(params(4)*y(4));
g1(3,1)=exp(y(1));
g1(3,2)=(-((1-params(4))*exp((1-params(4))*y(2)+y(3)+params(4)*y(4))));
g1(3,3)=(-exp((1-params(4))*y(2)+y(3)+params(4)*y(4)));
g1(3,4)=exp(y(4))-(exp(y(4))*(1-params(3))+params(4)*exp((1-params(4))*y(2)+y(3)+params(4)*y(4)));
g1(4,3)=1-params(8);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
