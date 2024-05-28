function [v0hat0,v1hat0,v0bar0,v1bar0,abar0] = update_v(v0,sys,cons)


K = sys.K; N = sys.N;
u = sys.u; r = sys.r; 
e0 = sys.e0; e1 = sys.e1;

v0hat0 = zeros(K,1); v0bar0 = zeros(K,1); abar0 = zeros(K,K);

v1hat0 = norm(v0 - r)^(-e1/2);
v1bar0 = norm(v0 - r)^(-e1/2);

for k = 1:K
    v0hat0(k) = norm(v0 - u(:,k))^(-e0/2);
    v0bar0(k) = norm(v0 - u(:,k))^(-e0/2);
    
    for j = 1:K
        if j ~= k
            if N == 0
                abar0(k,j) = cons.c0(k,j)*v0bar0(k)^2;
            else
                abar0(k,j) = cons.c0(k,j)*v0bar0(k)^2 + cons.c1(k,j)*v1bar0^2 + cons.c2(k,j)*v0bar0(k)*v1bar0;
            end
        end
    end
end


end %EOF