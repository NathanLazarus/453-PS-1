%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'R2_RBC_Dynare_sub';
M_.dynare_version = '5.3';
oo_.dynare_version = '5.3';
options_.dynare_version = '5.3';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e_z'};
M_.exo_names_tex(1) = {'e\_z'};
M_.exo_names_long(1) = {'e_z'};
M_.endo_names = cell(4,1);
M_.endo_names_tex = cell(4,1);
M_.endo_names_long = cell(4,1);
M_.endo_names(1) = {'lnc'};
M_.endo_names_tex(1) = {'lnc'};
M_.endo_names_long(1) = {'lnc'};
M_.endo_names(2) = {'lnn'};
M_.endo_names_tex(2) = {'lnn'};
M_.endo_names_long(2) = {'lnn'};
M_.endo_names(3) = {'lnz'};
M_.endo_names_tex(3) = {'lnz'};
M_.endo_names_long(3) = {'lnz'};
M_.endo_names(4) = {'lnk'};
M_.endo_names_tex(4) = {'lnk'};
M_.endo_names_long(4) = {'lnk'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'sigma'};
M_.param_names_tex(1) = {'sigma'};
M_.param_names_long(1) = {'sigma'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'eps'};
M_.param_names_tex(5) = {'eps'};
M_.param_names_long(5) = {'eps'};
M_.param_names(6) = {'c_y'};
M_.param_names_tex(6) = {'c\_y'};
M_.param_names_long(6) = {'c_y'};
M_.param_names(7) = {'k_y'};
M_.param_names_tex(7) = {'k\_y'};
M_.param_names_long(7) = {'k_y'};
M_.param_names(8) = {'rho_z'};
M_.param_names_tex(8) = {'rho\_z'};
M_.param_names_long(8) = {'rho_z'};
M_.param_names(9) = {'sigma_z'};
M_.param_names_tex(9) = {'sigma\_z'};
M_.param_names_long(9) = {'sigma_z'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 9;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 4;
M_.eq_nbr = 4;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 3 7;
 0 4 0;
 1 5 8;
 2 6 0;]';
M_.nstatic = 1;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , '3' ;
  4 , 'name' , 'lnz' ;
};
M_.mapping.lnc.eqidx = [1 2 3 ];
M_.mapping.lnn.eqidx = [1 2 3 ];
M_.mapping.lnz.eqidx = [1 2 3 4 ];
M_.mapping.lnk.eqidx = [1 2 3 ];
M_.mapping.e_z.eqidx = [4 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 4 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(4, 1), 'log_deflator', cell(4, 1), 'growth_factor', cell(4, 1), 'log_growth_factor', cell(4, 1));
M_.NNZDerivatives = [17; -1; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
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
M_.params(7) = K_SS/Y_SS;
k_y = M_.params(7);
i_y	= delta * K_SS/Y_SS;
M_.params(6) = 1-i_y;
c_y = M_.params(6);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 1.851102799953296;
oo_.steady_state(1) = 0.2704278202876532;
oo_.steady_state(2) = 0;
oo_.steady_state(3) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
steady;
options_.irf = 100;
options_.nograph = true;
options_.order = 1;
var_list_ = {'lnz';'lnc';'lnn';'lnk'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'R2_RBC_Dynare_sub_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
