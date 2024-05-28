
% run addpathyalmip.m
close all
clear all
clc

% load simulation parameters
sys = config(1);
sys.ops_soc = sdpsettings('solver','mosek','verbose',0,'debug',1) ;

n_chan = 1;
pmax_ris_vec = -10:10:20; % number of iterations
linestyle = {':k','--g','-.bs','-rd',':bo',':rs','-r*'};
opt_TX = 1; opt_loc = 1;

RIS_vec = [2];

for rr = 1:length(RIS_vec)
    RIS = RIS_vec(rr);
    
    % configurations of RIS
    [sys.N,sys.Na,sys.amax,sys.AA,sys.ONE] = config_RIS(RIS,sys.N0);
    
    minrate = zeros(length(pmax_ris_vec),n_chan);
    
    
    for pp = 1:length(pmax_ris_vec)
        pp
        pmax_ris_dBm = pmax_ris_vec(pp);
        sys.pmax_r = db2pow(pmax_ris_dBm);
        for cc = 1:n_chan
            %cc
            [minrate(pp,cc),minrate_vec] = optimize_alg(RIS,opt_TX,opt_loc,sys,cc);
        end
        
    end
    minrate_mean = mean(minrate,2);
    plot(pmax_ris_vec,minrate_mean, linestyle{RIS+1}, 'LineWidth', 2); hold on
end

xlabel('Maximum transmit power of the RIS $(p_{\mathrm{max}}^{\mathrm{RIS}})$ [dBm]','Interpreter','latex','FontSize',12)
ylabel('Minimum rate [nats/s/Hz]','Interpreter','latex','FontSize',12)
legend('No RIS','Passive RIS','Hybrid RIS, 4 active elements')