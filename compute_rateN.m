
function rate_min = compute_rateN(w,sys,chan,N,Ups)

K = sys.K; sigma2 = sys.sigma2; Nt = sys.Nt;
H1 = chan.H1; h0 = chan.h0; h2 = chan.h2;

%% effective channel
h = zeros(Nt,K);
for k = 1:K
%     if N > 0
%         h(:,k) = h0(:,k) + (h2(1:N,k)'*Ups(1:N,1:N)*H1(1:N,:))';
%     else
        h(:,k) = h0(:,k);
%     end
end

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

