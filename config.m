function sys = config(UE_fix)

%% System pasys.rametesys.rs
sys.Nt = 2; % # UAV tsys.ransmit antennas
sys.K = 4; % # UEs
sys.N0 = 32; % #RIS elements

%% Lasys.rge-scale pasys.rametesys.rs
sys.e0 = 3.2; % UAV-UE pathloss exponent
sys.e1 = 2.0; % UAV-RIS pathloss exponent
sys.e2 = 2.2; % RIS-UE pathloss exponent

sys.kappa0 = 0; % UAV-UE Rician factosys.r
sys.kappa1 = db2pow(100); % UAV-RIS Rician factosys.r 
sys.kappa2 = db2pow(10); % RIS-UE Rician factosys.r

sys.xi0 = db2pow(-30); % sys.refesys.rence path loss at d0 = 1m

%% Powesys.r and noise pasys.rametesys.rs
sigma2_u_dBm = -80;
sys.sigma2_u = db2pow(sigma2_u_dBm); % noise powesys.r at UE
sigma2_SI = 1.2*sys.sigma2_u; % RSI at active RIS elements
sys.sigma2_r = sys.sigma2_u + sigma2_SI; % total noise + RSI

sys.pmax = db2pow(20); % Watts, maximsys.um tsys.ransmit powesys.r of UAV
sys.pmax_r = db2pow(0); % maximsys.um tsys.ransmit powesys.r of RIS

%% positions of UAVs, RIS and UEs
sys.D = 50; % covesys.rage asys.rea sys.Dxsys.D m2
sys.dmin = 0; sys.dmax = sys.D; % min and max coosys.rdinates
sys.hmin = 100; sys.hmax = 220; % min and max height of UAV
sys.hv = sys.hmin;% initial height of UAV


if UE_fix == 0 % sys.regenesys.rate UE locations
    sys.u = [(sys.dmax-sys.dmin).*sys.rand(2,sys.K) + sys.dmin; zeros(1,sys.K)]; % positions of UEs
    save('UE_loc.mat','sys.u')
else % load UE locations
    usaved = load('UE_loc.mat'); sys.u = usaved.u;  sys.u(2,1) = 120; 
    if sys.D == 50
        sys.u = sys.u/4;
    end
end

sys.r = [sys.D/2;sys.D;50]; % positions of RIS

if UE_fix == 0
    plot_UE_location(sys)
end
