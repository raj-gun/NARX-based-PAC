function [Y_est,Y_est_sub] = model_simulation_clstr(model,u,y,match_ind)

a1=model{1};a2=model{2};b1=model{3};b2=model{4};theta=model{5};trm_chsn_lin=model{6};trm_chsn_lin_org=model{7};
n_lin_trms_org=model{8};nl_ord_max=model{9};trm_chsn_nl=model{10};bias=model{11};n_inpts=model{12};
%% Form the linear regressors/monomials/lagged-terms/model terms
min_dyn_ord_u = ones(1,n_inpts).*b1;
max_dyn_ord_u = ones(1,n_inpts).*b2;
min_dyn_ord_y = a1;
max_dyn_ord_y = a2;
inpt0=ones(1,n_inpts).*0;
n_terms_yi = max_dyn_ord_y-min_dyn_ord_y+1; %No. of lagged terms from each output
n_terms_ui = max_dyn_ord_u-min_dyn_ord_u+1; %No. of lgged terms from each input
n_terms_y = sum(n_terms_yi);
n_terms_u = sum(n_terms_ui);
%---------------------------------
[U_delay_mat_sim,~,~] = info_mat_sysID_i(min_dyn_ord_u,max_dyn_ord_u,u,y);
%---------------------------------
%---------------------------
y_lag_lin_srt = []; % Retained for sim_model_clstr interface compatibility
%---------------------------
[Y_est,Y_est_sub] = sim_model_clstr(theta,trm_chsn_lin,trm_chsn_lin_org,U_delay_mat_sim,n_lin_trms_org,y_lag_lin_srt,nl_ord_max,trm_chsn_nl,bias,match_ind);
%---------------------------
% error_sim = Y_sim - y_est;
% dat_len_y_sim = length(Y_sim);
% msse = (error_sim'*error_sim)/dat_len_y_sim;
% %%
% sse = [msse, mspe, mskpe];
% y_hat = [y_est, y_pred, Y_kSA];
% error = [error_sim, error_pred, error_kSA];
end