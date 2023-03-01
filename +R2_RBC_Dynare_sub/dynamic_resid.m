function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = R2_RBC_Dynare_sub.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(4, 1);
lhs = exp((-y(3)));
rhs = params(2)*exp((-y(7)))*(1+exp(y(8))*params(4)*exp((params(4)-1)*y(6))*exp((1-params(4))*y(4))-params(3));
residual(1) = lhs - rhs;
lhs = exp((-y(3)))*(1-params(4))*exp(y(5))*exp(params(4)*y(2));
rhs = exp(y(4)*(params(4)+params(5)));
residual(2) = lhs - rhs;
lhs = exp(y(3))+exp(y(6));
rhs = exp((1-params(4))*y(4)+y(5)+params(4)*y(2))+(1-params(3))*exp(y(2));
residual(3) = lhs - rhs;
lhs = y(5);
rhs = params(8)*y(1)+x(it_, 1);
residual(4) = lhs - rhs;

end
