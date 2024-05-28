clc
% clear all
% config

sys = config(1);
sys.N0 = 128;

data_location = strcat('./channels/',num2str(sys.N0));
if ~exist(data_location, 'dir')
    mkdir(data_location)
end

n_chan = 20;

for ii = 1:n_chan
    
    % UAV-UE small-scale channel
    g0_LoS = sqrt(sys.kappa0/(sys.kappa0+1))*ones(sys.Nt,sys.K);
    g0_NLoS = sqrt(1/(sys.kappa0+1))*(1/sqrt(2)*(randn(sys.Nt,sys.K) + 1i*randn(sys.Nt,sys.K)));
    g0 = g0_LoS + g0_NLoS;% Rician fading
    
    % UAV-RIS small-scale channel
    g1_NLoS = sqrt(1/(sys.kappa1+1))*1/sqrt(2)*(randn(sys.N0,sys.Nt) + 1i*randn(sys.N0,sys.Nt));
    g1_LoS = sqrt(sys.kappa1/(sys.kappa1+1))*LoS_channel(sys.Nt, sys.N0, 'UAV-RIS');
    G1 = g1_LoS + g1_NLoS;% Rician fading
    
    % RIS-UE small-scale channel
    g2 = zeros(sys.N0,sys.K);
    for k = 1:sys.K
        g2_NLoS = sqrt(1/(sys.kappa2+1))*1/sqrt(2)*(randn(1,sys.N0) + 1i*randn(1,sys.N0));
        g2_LoS = sqrt(sys.kappa2/(sys.kappa2+1))*LoS_channel(sys.N0, 1, 'RIS-UE');
        g2(:,k) = (g2_LoS + g2_NLoS).';% Rician fading
    end
    
    % save large scale params
    file_name = strcat(data_location,'/small_scale',num2str(ii));
    save(file_name,'g0','G1','g2');
    %save(file_name,'g0','G1','g2','err0','Err1','err2');
    
end

% end % EOF