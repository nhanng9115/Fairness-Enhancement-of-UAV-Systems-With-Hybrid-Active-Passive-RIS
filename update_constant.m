function [c0,c1,c2,c3] = update_constant(Ups,w,sys,chan)

K = sys.K; N = sys.N; Na = sys.Na;
g0 = chan.g0; G1 = chan.G1; xi0 = sys.xi0; h2 = chan.h2;
sigma2_u = sys.sigma2_r; 

c0 = zeros(K,K); c1 = zeros(K,K); c2 = zeros(K,K);

for k = 1:K
    scale = 1/sigma2_u;
    c0(k,k) = scale * abs(g0(:,k)'*w(:,k))^2*xi0;
    if N > 0
        c1(k,k) = scale * abs(h2(:,k)'*Ups*G1*w(:,k))^2*xi0;
        c2(k,k) = scale * 2*xi0*real(w(:,k)'*g0(:,k)*h2(:,k)'*Ups*G1*w(:,k));
    end
    
    for j = 1:K
        if j ~= k
            c0(k,j) = scale * abs(g0(:,k)'*w(:,j))^2*xi0;
            if N > 0
                c1(k,j) = scale * abs(h2(:,k)'*Ups*G1*w(:,j))^2*xi0;
                c2(k,j) = scale * 2*xi0*real(w(:,j)'*g0(:,k)*h2(:,k)'*Ups*G1*w(:,j));
            end
        end
    end
end


c3 = zeros(N,1);
if Na > 0
    for nn = 1:Na
        c3(nn) = xi0*norm(G1(nn,:))^2*norm(w,'fro')^2;
    end
end

end % EOF