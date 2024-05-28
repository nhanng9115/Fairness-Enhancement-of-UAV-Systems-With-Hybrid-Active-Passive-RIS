function array_res = gen_array_res(angle, M)
% M = so antennas

array_res = 1/sqrt(M) * exp(-1i * [0:M-1] .* pi * sin(angle));
array_res = array_res.';
% array_res = 1/sqrt(M) * exp(-1i * [0:M-1] .* pi * sin(angle) * randn(M,1));

end