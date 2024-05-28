function H = LoS_channel(Nt, Nr, channel)
Nx = 8;
H = zeros(Nr,Nt);
rand_angle = rand;
if strcmp(channel, 'RIS-UE')
    AoD = 2*pi * rand_angle;
    AoA_az = 2*pi * rand_angle;
    AoA_el = 2*pi * rand_angle;%pi*rand_angle - pi/2;
    
    a_rx = exp(1i * [0:Nr-1] .* pi * sin(AoD)).';
    
    sinsin = sin(AoA_az)*sin(AoA_el);
    sincos = sin(AoA_az)*cos(AoA_el);
    n = [0:Nt-1];
    a_tx = exp(1i * pi * ( floor(n./Nx).*sinsin + (n-floor(n./Nx)*Nx) * sincos )).';
    H = a_rx*a_tx';
    
elseif strcmp(channel, 'UAV-RIS')
    
    AoD = 2*pi * rand_angle;
    a_tx = exp(1i * [0:Nt-1] .* pi * sin(AoD)).';
    
    AoA_az = 2*pi * rand_angle;
    AoA_el = pi*rand_angle - pi/2;
    sinsin = sin(AoA_az)*sin(AoA_el);
    sincos = sin(AoA_az)*cos(AoA_el);
    n = [0:Nr-1];
    a_rx = exp(1i * pi * ( floor(n./Nx).*sinsin + (n-floor(n./Nx)*Nx) * sincos )).';
    
    H = a_rx*a_tx';
    
else % 'UE-AP'
    AoD = 2*pi * rand_angle;
    AoA = 2*pi * rand_angle;
    
    a_rx = exp(1i * [0:Nr-1] .* pi * sin(AoA)).';
    a_tx = exp(1i * [0:Nt-1] .* pi * sin(AoD)).';
    
    H = a_rx*a_tx';
end
end