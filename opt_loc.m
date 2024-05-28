%% Optimize variable v
function [obj_value,v_new,v0hat_new,v1hat_new,v0bar_new,v1bar_new,abar_new] = ...
    opt_loc(v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,Ups_old,abar_old,sys,cons)

ZERO = 0;

ops_soc = sys.ops_soc;
N = sys.N; K = sys.K; Na = sys.Na; pmax_r = sys.pmax_r;
u = sys.u; r = sys.r; AA = sys.AA; e0 = sys.e0; e1 = sys.e1; 
dmin = sys.dmin; dmax = sys.dmax; hmin = sys.hmin; hmax = sys.hmax;
c0 = cons.c0; c1 = cons.c1; c2 = cons.c2; c3 = cons.c3; 
sigma2_r = sys.sigma2_r; sigma2 = sys.sigma2; sigma2_u = sys.sigma2_u;

yalmip('clear') % clearing YALMIPs internal database

%% Variables
v = sdpvar(3,1,'full','real');
tau = sdpvar(1,1,'full','real');
v0hat = sdpvar(K,1,'full','real');
v1hat = sdpvar(1,1,'full','real');
v0bar = sdpvar(K,1,'full','real');
v1bar = sdpvar(1,1,'full','real');

x0 = sdpvar(K,1,'full','real');
x1 = sdpvar(1,1,'full','real');

ahat = sdpvar(K,1,'full','real');
abar = sdpvar(K,K,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];

F = [F, tau >= ZERO, ahat >= ZERO, abar >= ZERO];
F = [F, v0hat >= ZERO, v1hat >= ZERO, v0bar >= ZERO, v1bar >= ZERO];
F = [F, x0 >= ZERO, x1 >= ZERO];
F = [F, v([1,2]) <= dmax, v([1,2]) >=  dmin];
% F = [F, hmin <= v(3,:) <= hmax]; % (9b)
F = [F, v(3) == hmin]; % (9b)

% if N > 0 && UBtau_BF == 1
%     F = [F, tau <= objRIS];
% end

%% (32)
for k = 1:K
    % (32a)
    F00_pow = get_Fpow(v0hat(k),-2/e0,v0hat_old(k)); % (32a)
    F = [F, norm(v - u(:,k)) + F00_pow <= 0]; % (32a)
    
    % (32b) ********************************************* 3 dòng dưới này
    % là biến đổi của (32b) nhưng em không nhớ tại sao lại làm thế nữa, chỉ
    % nhớ là đặt vế trái của (32b) lớn hơn hoặc bằng x0(k). Nay xem lại em
    % lại không hiểu thế có đúng không ??????????????
    F0_qua = get_Fqua(v,u(:,k),v_old);
    F = [F, F0_qua <= 0];
    F01_pow = get_Fpow(v0bar(k),4/e0,v0bar_old(k)); % for (46)
    F = [F, F0_qua <= x0(k)];
    F = [F, cone([1,0.5*(-F01_pow - x0(k))],0.5*(-F01_pow + x0(k)))];
    
    % (32e)
    Fquav0hat = get_Fqua(v0hat(k),0,v0hat_old(k));
    Fquav1hat = get_Fqua(v1hat,0,v1hat_old);
    if c2(k,k) > 0
        Fbilv0v1hat = get_Fbil(v0hat(k),v1hat,c2(k,k),v0hat_old(k),v1hat_old);
    else
        Fbilv0v1hat = get_Fbil(v0bar(k),v1bar,c2(k,k),v0bar_old(k),v1bar_old);
    end
    F = [F, ahat(k) + c0(k,k)*Fquav0hat + c1(k,k)*Fquav1hat + abs(c2(k,k))*Fbilv0v1hat <= 0];
    
    % (32f)
    for j = 1:K
        if j ~= k
            if c2(k,j) > 0
                Fbilv0v1bar = get_Fbil(v0bar(k),v1bar,c2(k,j),v0bar_old(k),v1bar_old);
            else
                Fbilv0v1bar = get_Fbil(v0hat(k),v1hat,c2(k,j),v0hat_old(k),v1hat_old);
            end
            F = [F, abar(k,j) >= c0(k,j)*v0bar(k)^2 + c1(k,j)*v1bar^2 + abs(c2(k,j))*Fbilv0v1bar];
        end
    end
end

%% (32c,d) do not have index k
% (32c)
F10_pow = get_Fpow(v1hat,-2/e1,v1hat_old);
F = [F, norm(v - r) + F10_pow <= 0];

% (32d)
F1_qua = get_Fqua(v,r,v_old); % for (46)
F = [F, F1_qua <= 0];
F11_pow = get_Fpow(v1bar,4/e1,v1bar_old); % for (47)
F = [F, F1_qua <= x1];
F = [F, cone([1,0.5*(-F11_pow - x1)],0.5*(-F11_pow + x1))];

%% (27)
for k = 1:K
    set_j = [1:k-1,k+1:K];
    
    % compute first rate term
    rate_1 = log(sigma2(k)/sigma2_u + ahat(k) + sum(abar(k,set_j))); % RHS of (27)
    
    % compute rtilde_ub in (27)
    denum_Dmk = sigma2(k)/sigma2_u + sum(abar_old(k,set_j));

    % gradient 
    sum_Dmj = sum(abar(k,set_j) - abar_old(k,set_j));

    rate_2 = log(denum_Dmk) + sum_Dmj/denum_Dmk;
    F = [F,  rate_1  >= tau + rate_2];
end
%% ---------------------------------------------------------------------------------


%% (33)
if Na > 0
    p_ris = 0;
    for n = 1:Na
        aa = AA(n);
        a_n = Ups_old(aa,aa);
        xi_n = sigma2_r + c3(aa)*(v1bar)^2;
        p_ris = p_ris + abs(a_n)^2*xi_n;
    end
    F = [F, p_ris <= pmax_r];
end


%% Get results
% disp('*********** TRAJECTORY **************************');
optimize(F,-obj,ops_soc);
obj_value = double(obj);
v_new = double(v);
v0hat_new = double(v0hat);
v1hat_new = double(v1hat);
v0bar_new = double(v0bar);
v1bar_new = double(v1bar);
abar_new = double(abar);

end % EOF




