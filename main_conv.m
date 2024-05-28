
% run addpathyalmip.m
close all
clear% all
clc


%% Load simulation parameters
UE_fix = 1;% = 1: use the generated locations
sys = config(UE_fix);
n_chan = 1;
n_iter = 30; % number of iterations
sys.ops_soc = sdpsettings('solver','mosek','verbose',0,'debug',1) ;

opt_location = 1; opt_TX_BF = 1;

% RIS = 1: passive RIS, RIS = 2: hybrid active-passive RIS, RIS = 3: full passive RIS
obj_TXBF = zeros(n_iter,n_chan,3); obj_RIS = zeros(n_iter,n_chan,3); obj_loc = zeros(n_iter,n_chan,3);
minrate = zeros(n_iter,n_chan,3);  minrate_1 = zeros(n_iter,n_chan,3);
linestyle = {':ko','-g^','-rs','--r','-r^','-cs','--m'};
G = db2pow(0);
RIS_vec = [2];

V = zeros(3,n_chan,length(RIS_vec));
for cc = 1:n_chan
    cc
    
    for rr = 1:length(RIS_vec)
        RIS = RIS_vec(rr);
        
        % configurations of RIS
        [sys.N,sys.Na,sys.amax,sys.AA,sys.ONE] = config_RIS(RIS,sys.N0);
        
        % load small scale-fading channels
        file_name = strcat('./channels/',num2str(sys.N0),'/small_scale',num2str(cc)); data = load(file_name);
        chan.g0 = G^2*data.g0; chan.G1 = G*data.G1; chan.g2 = G*data.g2;
        
        % Initialization
        [w_old,Ups_old,v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,abar_old,gamma_old,btilde_old,sys,chan,cons] = initialize(sys,chan);
        
        %% Start iterations
        pmax0 = sys.pmax/5;
        
        run_RIS = 1; run_TX = 1; run_loc = 1; objRIS = 0; objBF = 0; stop_TXBF = 0;
        for ii = 1:n_iter
            ii
            if ii <= 5
                sys.pmax = pmax0*ii;
                sys.pmax
            end
            %% UAV TX beamforming -------------------------------------------------------------------
            if opt_TX_BF == 1 && run_TX == 1
                [obj_TXBF(ii,cc,RIS+1),w_new,gamma_new] = opt_tx_BF(w_old,gamma_old,Ups_old,sys,chan);
                %w_new
                w_old = w_new; gamma_old = gamma_new;
                if sys.N > 0
                    [cons.Q,cons.Qtilde,cons.t,cons.ttilde,cons.e,cons.etilde,cons.Xi,~] = update_Q(w_old,Ups_old,sys,chan);
                end
                [cons.c0,cons.c1,cons.c2] = update_constant(Ups_old,w_old,sys,chan);
                objBF = obj_TXBF(ii,cc,RIS+1);
            end
            
            
            %% trajectory --------------------------------------------------------------------------------------------
            if opt_location == 1 && run_loc == 1
                [obj_loc(ii,cc,RIS+1),v_new,v0hat_new,v1hat_new,v0bar_new,v1bar_new,abar_new] = ...
                    opt_loc(v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,Ups_old,abar_old,sys,cons);
                v_old = v_new;
                
                abar_old = abar_new;
                [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
                [cons.c0,cons.c1,cons.c2] = update_constant(Ups_old,w_old,sys,chan);
                [v0hat_old,v1hat_old,v0bar_old,v1bar_old] = update_v(v_old,sys,cons);
            end
            
            %% RIS --------------------------------------------------------------------------------------------
            if RIS > 0 && run_RIS == 1
                [obj_RIS(ii,cc,RIS+1),Ntilde_new,Ups_new] = opt_RIS(Ups_old,btilde_old,sys,cons);
                Ups_old = Ups_new;
                %abs(diag(Ups_old))
                btilde_old = Ntilde_new;
                [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
                [cons.Q,cons.Qtilde,cons.t,cons.ttilde,cons.e,cons.etilde,cons.Xi,~] = update_Q(w_old,Ups_old,sys,chan);
                objRIS = obj_RIS(ii,cc,RIS+1);
            end
            
            
            %% Compute min rate
            minrate(ii,cc,RIS+1) = compute_rate(w_old,sys,chan);
            
        end
        V(:,cc,rr) = v_old;
        
        %if RIS == 0
        figure
        %end
        %plot(1:n_iter,obj_loc(:,cc,RIS+1),'-m^','LineWidth',1); hold on
        %plot(1:n_iter,obj_RIS(:,cc,RIS+1),'-b^','LineWidth',1); hold on
        %plot(1:n_iter,obj_TXBF(:,cc,RIS+1),'-co','LineWidth',1); hold on
        plot(1:n_iter,minrate(:,cc,RIS+1),linestyle{RIS+1},'LineWidth',2); hold on
    end
end
% plot_UE_location(sys,V,RIS_vec);

xlabel('iteration numbers','Interpreter','latex','FontSize',12)
ylabel('Minimum rate [nats/s/Hz]','Interpreter','latex','FontSize',12)

% for RIS = 2:3
%     figure
%     minrate_mean = mean(minrate(:,:,RIS+1),2);
%     plot(1:n_iter,minrate_mean,linestyle{RIS+1},'LineWidth',2); hold on
%     xlabel('iteration numbers')
%     ylabel('Minimum rate [nats/s/Hz]')
% end
%% Plot locations
% plot_UE_location(V)



