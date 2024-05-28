%% Result Proposed method
function [obj_value,p_new,w_new] = opt_power(p_old,Ups_old,pmax,w)

global K pmax_r sigma2_r Nt ops_soc sigma2 h0 H1 Na h
% pmax
yalmip('clear')
p = sdpvar(K,1,'full','real');
tau = sdpvar(1,1,'full','real');

obj = tau;
F = [];

F = [F, tau >= 0];
F = [F, p >= 0, sum(p) == pmax];

if Na > 0
    pris = 0;
    for n = 1:Na
        h_1n = H1(n,:); a_n = Ups_old(n,n);
        xi_n = sigma2_r + norm(h_1n)^2 * sum(p);
        pris = pris + abs(a_n)^2 * xi_n;
    end
    F = [F, pris <= pmax_r]; % (11)
end

%% (59), another way, based on Q. Wu paper
for k = 1:K

    % compute first rate term
    sum_power = sigma2(k) + p(k)*abs(h(:,k)'*w(:,k))^2;
    for j = 1:K
        if j~= k
            sum_power = sum_power + p(j)*abs(h(:,k)'*w(:,j))^2;
        end
    end
    rate_1 = log(sum_power);
    
    % compute second rate term
    denum_Dmk = sigma2(k);
    for j = 1:K
        if j ~= k
            h_kj = abs(h(:,k)'*w(:,j))^2;
            denum_Dmk = denum_Dmk + p_old(j)*h_kj;
        end
    end
    
    % compute sum D_kj
    sum_Dmj = 0;
    for j = 1:K
        if j ~= k
            h_kj = abs(h(:,k)'*w(:,j))^2;
            D_mj = h_kj/denum_Dmk;
            sum_Dmj = sum_Dmj + D_mj*(p(j) - p_old(j));
        end
    end
    
    rate_2 = sum_Dmj + log(denum_Dmk);
    
    F = [F,  rate_1 >= tau + rate_2];
end

%% Get results
disp('*********** POWER **************************');

optimize(F,-obj,ops_soc)
obj_value = double(obj)
p_new = double(p);
% gamma_new = double(gamma);


w_new = zeros(Nt,K);
for k = 1:K
    w_new(:,k) = sqrt(p_old(k))*h0(:,k)/norm(h0(:,k));
end
end % EOF




