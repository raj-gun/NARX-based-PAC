function [Y_est,Y_est_sub] = sim_model_clstr(theta,trm_chsn_lin,trm_chsn_lin_org,U_delay_mat,n_lin_trms,y_lag_lin_srt,nl_ord,trm_chsn_nl,bias,mod_sub_cmp)

% Simulate ploynomial ARX and NARX type different equation models
% n_lin_trms - number of linear terms 
% trm_chsn_lin - term indecies of chosen linear terms from the superset of terms
% max_dyn_ord - maximum dynamic order of the original superset of terms 
% U_delay_mat - Matrix containing all the linear input lagged terms (superset)
% mod_sub_cmp - Sub-component of the model to be extracted
%% test
% n_lin_trms = length(trm_chsn);
% max_dyn_ord = 3;
% U_delay_mat = diff_eq_mat(nb,U,'u1');

%%
size_U_d_m = size(U_delay_mat);
dat_len = size_U_d_m(1);

[unq_nl_comb] = nl_term_comb(nl_ord,n_lin_trms); % Evalutate nonlinear combinations 
nl_terms = sum([unq_nl_comb{2:end,2}]); % Total no. of nonlinear terms
X_nl = zeros(dat_len,nl_terms);

X_sup_lin = U_delay_mat;

X_sub_lin = X_sup_lin(:,trm_chsn_lin); % Linear subset of terms 
X_sub_lin_org = X_sup_lin(:,trm_chsn_lin_org); % Initial linear subset of terms 
% ---------- Nonlinear regressors are formed from X_sub_lin -----------
if nl_ord >= 2
    nl_ind_end = 0;
    for n = 2:nl_ord
        unq_comb = unq_nl_comb{(n),1};
        [X_comp] = nl_reg_data_mat(X_sub_lin_org,unq_comb);
        
        nl_ind = nl_ind_end + unq_nl_comb{(n),2};
        
        X_nl(:,nl_ind_end+1:nl_ind) = X_comp;
        
        nl_ind_end = nl_ind ; 
    end
    
    X_main_sup = [X_sub_lin_org , X_nl];
    X_main = [X_sub_lin , X_main_sup(:,trm_chsn_nl)];
else
    X_main = X_sub_lin;           
end
% ---------------------------------------------------------------------

mod_sub_cmp = logical(mod_sub_cmp);

if isempty(X_main)
    Y_est = bias.*ones(dat_len,1); % Simulate model
    Y_est_sub = zeros(dat_len,1); % Simulation of just the extracted sub-component of the model
else
    Y_est = X_main*theta + bias.*ones(dat_len,1); % Simulate model
    Y_est_sub = X_main(:,mod_sub_cmp)*theta(mod_sub_cmp); % Simulate model sub-component
end

end
