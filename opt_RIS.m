%% Result Proposed method
function [obj_value,btilde_new,Ups_new] = opt_RIS(Ups_old,btilde_old,sys,cons)

ops_soc = sys.ops_soc; K = sys.K; N = sys.N; Na = sys.Na;
amax = sys.amax; pmax_r = sys.pmax_r;% sigma2_r = sys.sigma2_r; AA = sys.AA;
Q = cons.Q; Qtilde = cons.Qtilde; t = cons.t; ttilde = cons.ttilde; Xi = cons.Xi; e = cons.e; etilde = cons.etilde;

ZERO = 0;
yalmip('clear') % clearing YALMIPs internal database

%% Variables
alpha = sdpvar(N,1,'full','complex');
tau = sdpvar(1,1,'full','real');
% b = sdpvar(N,1,'full','real');
btilde = sdpvar(K,1,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];
% F = [F, tau <= objBF];
% F = [F, tau >= objBF];

F = [F, tau >= 0, btilde >= ZERO];


%% (72)
F = [F, abs(alpha(1:Na)) <= amax];
F = [F, abs(alpha(Na+1:N)) <= 1];
if Na > 0
    F = [F, alpha'*Xi*alpha <= pmax_r];
end


%% (71)
alpha_old = diag(Ups_old(1:N,1:N));
for k = 1:K
    Qbar = Q(:,:,k) + Qtilde(:,:,k);
    qbar = t(:,k) + ttilde(:,k);
    ebar = e(k) + etilde(k);
    %rate1 = log(alpha'*Qbar*alpha + 2*real(alpha'*qbar) + ebar);
    
    Qtmp = sqrtm(Qbar);
    Fqua = -get_Fqua(Qtmp*alpha,zeros(N,1),Qtmp*alpha_old);
    rate1 = log(Fqua + 2*real(alpha'*qbar) + ebar);
    
    F = [F, btilde(k) >= alpha'*Qtilde(:,:,k)*alpha + 2*real(alpha'*ttilde(:,k))];
    D = real(btilde_old(k) + etilde(k));
    r_ub = log(D) + real(btilde(k) - btilde_old(k))/D;
    rate2 = real(r_ub);
    F = [F, rate1 >= tau + rate2];
end

%% Get results
% disp('*********** RIS **************************');
optimize(F,-obj,ops_soc);
obj_value = double(tau);
alpha_new = double(alpha);
% abs(alpha_new)
Ups_new = diag(alpha_new);
btilde_new = double(btilde);
end % EOF
