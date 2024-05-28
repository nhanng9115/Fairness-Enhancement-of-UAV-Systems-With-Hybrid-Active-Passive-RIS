function [N,Na,amax,AA,ONE] = config_RIS(RIS_type,Ntmp)

if RIS_type == 0 % No RIS
    N = 0; Na = 0;
elseif RIS_type == 1 % Passive RIS
    N = Ntmp; Na = 0;
elseif RIS_type == 2 % Passive RIS
    N = Ntmp; Na = 2;
elseif RIS_type == 3 % Hybrid RIS
    N = Ntmp; Na = 4;
elseif RIS_type == 4 % Hybrid RIS
    N = Ntmp; Na = 8;
elseif RIS_type == 5 % Hybrid RIS
    N = Ntmp; Na = 16;
else %  Active RIS
    N = Ntmp; Na = N;
end
if Na > 0
    amax = db2pow(40);
else
    amax = 1;
end

% %% determine active elements
% if fix_RIS == 1
    AA = 1:Na; elements = zeros(N,1); elements(AA) = 1; ONE = diag(elements);
% else
%     norm_h1 = vecnorm(H1);
%     AA = maxk(norm_h1,Na);
% end

end % EOF