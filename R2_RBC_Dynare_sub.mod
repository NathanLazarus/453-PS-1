// Solev RBC with labor
// Recitation 2
% endogenous variables
var
lnc, lnn, lnz, lnk;

% exogenous variables
varexo
e_z;

// Model Parameters
parameters
sigma, beta, delta, alpha, eps, c_y, k_y, rho_z, sigma_z;

% set parameter values
load param_RBC;
set_param_value('sigma',sigma)
set_param_value('beta',beta)
set_param_value('delta',delta)
set_param_value('alpha',alpha)
set_param_value('eps',eps)
set_param_value('rho_z',rho_z)
set_param_value('sigma_z',sigma_z)

K_SS	= (alpha/(1/beta - (1-delta)))^(1/(1-alpha));
Y_SS	= K_SS^alpha;
k_y = K_SS/Y_SS;
i_y	= delta * K_SS/Y_SS;
c_y	= 1 - i_y;

// Model Solution
model;
% Euler
exp(-lnc) = beta*exp(-lnc(+1))*(exp(lnz(+1))*alpha*exp((alpha-1)*lnk)*exp((1-alpha)*lnn)+1-delta);

% Labor Supply
(1-alpha)*exp(lnz)*exp(alpha*lnk(-1))*exp(-lnc)=exp((eps+alpha)*lnn);

% Resource const
exp(lnc)+exp(lnk)=exp(lnz+alpha*lnk(-1)+(1-alpha)*lnn)+(1-delta)*exp(lnk(-1));

% shocks
lnz=rho_z*lnz(-1)+e_z;
end;

initval;
lnk = log(6.366837);
lnc = log(1.310525);
lnn = 0;
lnz = 0;
end;

shocks;

var e_z; stderr 1;

end;

steady;

stoch_simul(order=1,irf=100,nograph) lnz lnc lnn lnk;