function [sigma2,h0,H1,h2,h,PL2] = update_channel(v,Ups,sys,chan)

K = sys.K; N = sys.N; Na = sys.Na; Nt = sys.Nt;
u = sys.u; r = sys.r; 
g0 = chan.g0; G1 = chan.G1; g2 = chan.g2;
e0 = sys.e0; e1 = sys.e1; e2 = sys.e2; xi0 = sys.xi0; ONE = sys.ONE;
sigma2_r = sys.sigma2_r; sigma2_u = sys.sigma2_u; 

%% Channels h0 h1 h2
h0 = zeros(Nt,K); H1 = zeros(N,Nt); h2 = zeros(N,K); PL2 = zeros(K,1);
for k = 1:K
    d0 = norm(u(:,k) - v);
    PL0 = xi0*d0^(-e0);
    h0(:,k) = sqrt(PL0)*g0(:,k);
    
    if N > 0
        d1 = norm(v - r);
        PL1 = xi0*d1^(-e1);
        H1 = sqrt(PL1)*G1;
        
        d2 = norm(u(:,k) - r);
        PL2(k) = xi0*d2^(-e2);
        h2(:,k) = sqrt(PL2(k))*g2(:,k);
    end
end

%% effective channel
h = zeros(Nt,K);
for k = 1:K
    h(:,k) = h0(:,k) + (h2(:,k)'*Ups*H1)';
end

%% effective noise
sigma2 = zeros(K,1);
for k = 1:K
    if N > 0 && Na > 0
        sigma2(k) = sigma2_u + sigma2_r*norm(h2(:,k)'*Ups*ONE)^2;
    else
        sigma2(k) = sigma2_u;
    end
end

end % EOF