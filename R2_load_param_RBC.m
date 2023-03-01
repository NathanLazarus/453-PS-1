%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Parameters
% For 14.453, Recitation 2
%
% Name: R2_load_param_RBC.m
% This version: 02/16/2023
% Modified by Shinnosuke Kikuchi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha =  .333; % share of capital
delta =  .025; %.08; % depreciation
beta  =  .984; % discount factor
sigma =  1;  % 1/IES
A     =  1; % original SS TFP (use A instead of z to clarify)

% separate disutility of labor:  n^(1+eps)/(1+eps)
% eps = 1/ frisch elasticity of labor
eps   = .25; 

% TFP shock process
shock   = 0.05; % .01; % percent shock
rho_z     = 0.9; % 0.97;
sigma_z =  1.2;

% Tol
tol=1e-7;

% Grid size. Decreasing this may make policy functions ugly
gridsize=500;

save param_RBC alpha delta beta sigma A eps rho_z sigma_z tol gridsize