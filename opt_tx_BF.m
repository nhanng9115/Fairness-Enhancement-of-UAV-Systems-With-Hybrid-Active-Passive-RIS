%% Optimize variable w
function [obj_value,w_new,gamma_new] = opt_tx_BF(w_old_mat,gamma_old,Ups_old,sys,chan)

K = sys.K; Nt = sys.Nt; N = sys.N; pmax_r = sys.pmax_r; Na = sys.Na; pmax = sys.pmax;
sigma2_u = sys.sigma2_u; sigma2_r = sys.sigma2_r; sigma2 = sys.sigma2;
h = chan.h; H1 = chan.H1; h0 = chan.h0;

ops_soc = sys.ops_soc; 

yalmip('clear')

w = sdpvar(Nt*K,1,'full','complex');
gamma = sdpvar(K,1,'full','real');
tau = sdpvar(1,1,'full','real');

obj = tau;
F = [];

F = [F, tau >= 0];

F = [F, norm(w)^2 <= pmax]; % (5c)

for k = 1:K
    F = [F, log(1 + gamma(k)) >= tau]; % (12a)
    %F = [F, log(sigma2(k)/sigma2_u + gamma(k)) >= tau]; % (12a)
    
    % convert w_old_mat (Nr x K) to w_old (NrK x 1)
    w_old = w_old_mat(:); % (13)

    % construct H_khat, H_kbar
    if N == 0
        hk = h0(:,k);
    else
        hk = h(:,k);
    end
    
    Hhat = zeros(K*Nt,K*Nt);
    Hhat(Nt*(k-1)+1:Nt*k,Nt*(k-1)+1:Nt*k) = hk*hk'/sigma2_u; % (14)
    Hbar = kron(eye(K),hk*hk')/sigma2_u;
    Hbar(Nt*(k-1)+1:Nt*k,Nt*(k-1)+1:Nt*k) = zeros(Nt,Nt); % (15)
    
    F_qol = w_old'*Hhat*w_old/gamma_old(k)^2 * gamma(k) - 2*real(w_old'*Hhat*w)/gamma_old(k); % (17)
    %F = [F, real(w_old'*Hhat*w) >= 0];
    F = [F, w'*Hbar*w + sigma2(k)/sigma2_u + F_qol <= 0]; % (18)
    
end


%% (19)
if Na > 0
    pris = 0;
    for nn = 1:Na
        h_1n = H1(nn,:); a_n = Ups_old(nn,nn);
        xi_n = sigma2_r + norm(h_1n)^2 * norm(w)^2;
        pris = pris + abs(a_n)^2 * xi_n;
    end
    F = [F, pris <= pmax_r]; % (11)
end


%% Get results
% disp('*********** TX BF **************************');

optimize(F,-obj,ops_soc);
obj_value = double(obj);
gamma_new = double(gamma);
w_new_mat = double(w);
% convert w back to matrix (Nt x K)
w_new = reshape(w_new_mat, Nt, []);
end % EOF




