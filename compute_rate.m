
function rate_min = compute_rate(w,sys,chan)

K = sys.K; sigma2 = sys.sigma2; h = chan.h; 
% update_channel(v,Ups);

rate = zeros(K,1);
for k = 1:K
    S_mk = abs(h(:,k)'*w(:,k))^2;
    I_mk = sigma2(k);
    for j = 1:K
        if j ~= k
            I_mk = I_mk + abs(h(:,k)'*w(:,j))^2;
        end
    end
    gamma(k) = S_mk/I_mk;
    rate(k) = log(1 + gamma(k));
    
end
rate_min = min(rate);

