%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBC model
% Solved by value function iteration
% For 14.453, Recitation 2
%
% Name: R2_RBC_VFI.m
% This version: 02/16/2023
% Originally by Adrien Auclert
% Modified by Shinnosuke Kikuchi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
profile on

%% Provide parameter values
run R2_load_param_RBC.m


%% Discretize the shock AR(1) into Markov process, S states
S = 10;  % Number of states - set to 1 to see the policy functions of the neoclassical growth model with labor supply
if S>1
    [log_z,Pr] = tauchen(S,0,rho_z,shock,sigma_z);
    z = A*exp(log_z);

    % Find the stationary distribution of the chain
    pr = ones(1,S)/S;                       
    diff=1;
    while diff>tol
        pri=pr*Pr;
        diff=max(abs(pri-pr));
        pr=pri;
    end
    zss=pr*z;
else
    z=1;
    zss=1;
    Pr=1;
    pr=1;
end

% maximum labor (just to compute upper bound on consumption)
nmax  = 100;
nmin  =.01;

%% Get Steady states
kss     = (((1/beta - 1 + delta)/alpha)/A)^(1/(alpha-1)); % compute steady state
css     = A*kss^alpha-delta*kss; % steady state consumption
yss     = A*kss^alpha;
xss     = yss-css;

% "new" steady state (for permanent shock case)
kss2    = (((1/beta - 1 + delta)/alpha)/A)^(1/(alpha-1)); 

% normalize coefficient psi so that ss labor equals 1...
psi = (1-alpha)*A*kss^alpha*css^(-sigma);

%% Compute c+ and n+ functions over (k,k',z)
tic
i=1; crit = 1;

% Grid for capital today
kt = linspace(0.5*kss,2*kss, gridsize-1);
ikss=find(kt<kss, 1, 'last');
k=[kt(1:ikss), kss, kt(ikss+1:end)];
ikss=ikss+1;

% grid for capital tomorrow
h = k;

% Fill in the vectors c_val(k,k',z), n_val(k,k',z) with values for
% consumption and labor
nfun=@(cz, zz, kz, hz) max(nmin,min(nmax,(((1-alpha)*zz*kz^(alpha)*cz^(-sigma))/psi)^(1/(eps+alpha))));
cfun=@(cz, zz, kz, hz) (cz>0)*(cz+hz-(1-delta)*kz-zz*kz^alpha*(nfun(cz,zz,kz,hz))^(1-alpha))+(-50)*(cz<=0);

n_val=zeros(gridsize, gridsize, S);
c_val=zeros(gridsize, gridsize, S);
for s=1:S
    s
    % initiate
        fprintf('-');
        
        cfunzhk=@(cz) cfun(cz,z(s), k(1), h(1));
        ctest=fzero(cfunzhk, css); % steady-state is our best guess here
        c_val(1,1,s)=ctest;
        n_val(1,1,s)=nfun(ctest,z(s), k(1), h(1));
        for j=2:gridsize  
            cfunzhk=@(cz) cfun(cz,z(s), k(1), h(j));
            ctest=fzero(cfunzhk, c_val(1,j-1,s));
            c_val(1,j,s)=ctest;
            n_val(1,j,s)=nfun(ctest,z(s), k(1), h(j));
        end
        
    for i=2:gridsize
        fprintf('-');
        if mod(i,50)==0
            fprintf('\n');
        end
        
        cfunzhk=@(cz) cfun(cz,z(s), k(i), h(1));
        ctest=fzero(cfunzhk, css);
        c_val(i,1,s)=ctest;
        n_val(i,1,s)=nfun(ctest,z(s), k(i), h(1));
        for j=2:gridsize  
            cfunzhk=@(cz) cfun(cz,z(s), k(i), h(j));
            ctest=fzero(cfunzhk, c_val(i,j-1,s));
            c_val(i,j,s)=ctest;
            n_val(i,j,s)=nfun(ctest,z(s), k(i), h(j));
        end
    end
end

save cnval c_val n_val


toc
%% Use Value function iteration
tic

diff=1; 
polC = zeros(gridsize,S); polN = zeros(gridsize,S); polH=zeros(gridsize,S);
% Initial guess is: value of being forced to stay at k forever
V=zeros(gridsize, S);

for s=1:S
    for i=1:gridsize
        V(i,s)=1/(1-beta)*(log(c_val(i,i,s))-psi*n_val(i,i,s)^(1+eps)/(1+eps));
    end
end

EV=V;

while diff>tol
for s=1:S
    % form objective function
    F=log(c_val(:,:,s))-psi*n_val(:,:,s).^(1+eps)/(1+eps)+beta*repmat(EV(:,s)',gridsize,1);
    % Find policy for future capital
    [V(:,s), polH(:,s)] = max(F,[],2);  
end
% Update V
EVi = V*Pr';
diff = max(max(abs(EVi-EV))) %display(diff)
EV = EVi;
end
toc
%% Recover policy functions
for s=1:S
    for i=1:gridsize
        polC(i,s)=c_val(i,polH(i,s),s);
        polN(i,s)=n_val(i,polH(i,s),s);
    end
end

%% Figures
% Value Function
figure
plot(k,V(:,1),'LineWidth',3,'MarkerSize',10)
hold on
plot(k,V(:,end), 'r','LineWidth',3,'MarkerSize',10)
plot([kss kss], ylim, 'k-.','LineWidth',3,'MarkerSize',10)
xlabel('k','Fontsize',16)
ylabel('V','Fontsize',16)
legend(['V(k, z=',num2str(z(1)),')'], ['V(k, z=',num2str(z(end)),')'], 'Non-stoch. SS', 'location', 'SouthEast')
set(gca,'FontSize',16);
filename=['fig_R2_RBC_VFI_vfunc', '.png'];
saveas(gcf,filename,'png')

% Policy Function for Capital
figure
plot(k,k(polH(:,1)),'LineWidth',3,'MarkerSize',10)
hold on
plot(k,k(polH(:,end)), 'r','LineWidth',3,'MarkerSize',10)
plot(k, k,'k','LineWidth',3,'MarkerSize',10)
plot([kss kss], ylim, 'k-.','LineWidth',2,'MarkerSize',10)
plot(xlim,[kss kss], 'k-.','LineWidth',2,'MarkerSize',10)
xlabel('k','Fontsize',16)
ylabel('k''','Fontsize',16)
legend(['k''(k, z=',num2str(z(1)),')'], ['k''(k, z=',num2str(z(end)),')'], '45° line', 'Non-stoch. SS', 'location', 'SouthEast')
set(gca,'FontSize',16);
filename=['fig_R2_RBC_VFI_pfunc_kprime', '.png'];
saveas(gcf,filename,'png')

% Policy Function for Consumption
figure
plot(k,polC(:,1),'LineWidth',3,'MarkerSize',10)
hold on
plot(k,polC(:,end), 'r','LineWidth',3,'MarkerSize',10)
plot([kss kss], ylim, 'k-.','LineWidth',2,'MarkerSize',10)
plot(xlim, [css css],  'k-.','LineWidth',2,'MarkerSize',10)
xlabel('k','Fontsize',16)
ylabel('c','Fontsize',16)
legend(['c(k, z=',num2str(z(1)),')'], ['c(k, z=',num2str(z(end)),')'], 'Non-stoch. SS', 'location', 'SouthEast')
set(gca,'FontSize',16);
filename=['fig_R2_RBC_VFI_pfunc_c', '.png'];
saveas(gcf,filename,'png')


% Policy Function for Labor
figure
plot(k,polN(:,1),'LineWidth',3,'MarkerSize',10)
hold on
plot(k,polN(:,end), 'r','LineWidth',3,'MarkerSize',10)
plot([kss kss], ylim, 'k-.','LineWidth',2,'MarkerSize',10)
plot(xlim, [1 1],  'k-.','LineWidth',2,'MarkerSize',10)
xlabel('k','Fontsize',16)
ylabel('n','Fontsize',16)
legend(['n(k, z=',num2str(z(1)),')'], ['n(k, z=',num2str(z(end)),')'], 'Non-stoch. SS', 'location', 'NorthEast')
set(gca,'FontSize',16);
filename=['fig_R2_RBC_VFI_pfunc_n', '.png'];
saveas(gcf,filename,'png')