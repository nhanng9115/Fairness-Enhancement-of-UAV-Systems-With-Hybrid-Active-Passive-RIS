function [minrate1,minrate_vec] = optimize_alg(RIS,opt_TX_BF,opt_location,sys,cc)

% global N Na
% sys.pmax = db2pow(pmax_dBm);

%% load small scale-fading channels
if RIS == 0
    file_name = strcat('./channels/',num2str(32),'/small_scale',num2str(cc));
else
    file_name = strcat('./channels/',num2str(sys.N),'/small_scale',num2str(cc));
end
data = load(file_name);
chan.g0 = data.g0; chan.G1 = data.G1; chan.g2 = data.g2;

%% Initialization
[w_old,Ups_old,v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,abar_old,gamma_old,Ntilde_old,sys,chan,cons] = initialize(sys,chan);

run_RIS = 1; run_TX = 1; run_loc = 1;
minrate = 0; ii = 0;
minrate0 = -Inf;

while (abs(minrate - minrate0) >= 1e-3 && ii <= 50) || (ii <= 5)
    
    minrate0 = minrate;  ii = ii + 1;
    
    %% UAV TX beamforming -------------------------------------------------------------------
    if opt_TX_BF == 1 && run_TX == 1% && ii <= 5
        [obj_TXBF,w_new,gamma_new] = opt_tx_BF(w_old,gamma_old,Ups_old,sys,chan);
        %w_new
        w_old = w_new; gamma_old = gamma_new;
        if sys.N > 0
            [cons.Q,cons.Qtilde,cons.t,cons.ttilde,cons.e,cons.etilde,cons.Xi,~] = update_Q(w_old,Ups_old,sys,chan);
        end
        [cons.c0,cons.c1,cons.c2] = update_constant(Ups_old,w_old,sys,chan);
    end
    
    
    %% trajectory --------------------------------------------------------------------------------------------
    if opt_location == 1 && run_loc == 1
        [obj_loc,v_new,v0hat_new,v1hat_new,v0bar_new,v1bar_new,abar_new] = ...
            opt_loc(v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,Ups_old,abar_old,sys,cons);
        v_old = v_new;
        if v_new(3) < 50
            disp('*********INFEASIBILITY************');
            ii
            break;
        end
        abar_old = abar_new;
        [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
        [cons.c0,cons.c1,cons.c2] = update_constant(Ups_old,w_old,sys,chan);
        [v0hat_old,v1hat_old,v0bar_old,v1bar_old] = update_v(v_old,sys,cons);
    end
    
    %% RIS --------------------------------------------------------------------------------------------
    if RIS > 0 && run_RIS == 1
        [obj_RIS,Ntilde_new,Ups_new] = opt_RIS(Ups_old,Ntilde_old,sys,cons);
        Ups_old = Ups_new;
        %abs(diag(Ups_old))
        Ntilde_old = Ntilde_new;
        [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
        [cons.Q,cons.Qtilde,cons.t,cons.ttilde,cons.e,cons.etilde,cons.Xi,~] = update_Q(w_old,Ups_old,sys,chan);
    end
    
    
    %% Compute min rate
    minrate = compute_rate(w_old,sys,chan);
    minrate_vec(ii) = minrate;
    %minrate_vec.'
end

%% recompute min-rate based on true channel
eps = 0;
if eps > 0
    [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = get_true_channel(v_old,Ups_old,sys,chan,eps);
    minrate = compute_rate(w_old,sys,chan);
end

minrate1 = minrate;

end % EOF