%% Optimize variable v
function [obj_value,v_new,v0hat_new,v1hat_new,v0bar_new,v1bar_new,abar_new] = ...
    opt_loc_noRIS(v_old,v0hat_old,v1hat_old,v0bar_old,v1bar_old,Ups_old,abar_old)

ZERO = 0;

global ops_soc K Na pmax_r e0 e1 u r AA dmin dmax c0 c1 c2 c3 hmin hmax sigma2_r sigma2 sigma2_u

yalmip('clear') % clearing YALMIPs internal database

%% Variables
v = sdpvar(3,1,'full','real');
tau = sdpvar(1,1,'full','real');
v0hat = sdpvar(K,1,'full','real');
% v1hat = sdpvar(1,1,'full','real');
v0bar = sdpvar(K,1,'full','real');
% v1bar = sdpvar(1,1,'full','real');

x0 = sdpvar(K,1,'full','real');

ahat = sdpvar(K,1,'full','real');
abar = sdpvar(K,K,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];

F = [F, tau >= ZERO, ahat >= ZERO, abar >= ZERO];
F = [F, v0hat >= ZERO, v0bar >= ZERO];
% F = [F, , x0 >= ZERO, x1 >= ZERO];
F = [F, v([1,2]) <= dmax, v([1,2]) >=  dmin];
% F = [F, hmin <= v(3,:) <= hmax]; % (9b)
F = [F, v(3) == hmin]; % (9b)


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
    F01_pow = get_Fpow(v0bar(k),4/e0,v0bar_old(k)); % for (46)
    F = [F, F0_qua <= x0(k)];
    F = [F, cone([1,0.5*(-F01_pow - x0(k))],0.5*(-F01_pow + x0(k)))];
    
    % (32e)
    Fquav0hat = get_Fqua(v0hat(k),0,v0hat_old(k));
    F = [F, ahat(k) + c0(k,k)*Fquav0hat <= 0];
    
    % (32f)
    for j = 1:K
        if j ~= k
            F = [F, abar(k,j) >= c0(k,j)*v0bar(k)^2];
        end
    end
end

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

%% Get results
disp('*********** TRAJECTORY **************************');
optimize(F,-obj,ops_soc)
obj_value = double(obj)
v_new = double(v);
v0hat_new = double(v0hat);
v1hat_new = 0;
v0bar_new = double(v0bar);
v1bar_new = 0;
abar_new = double(abar);

end % EOF




