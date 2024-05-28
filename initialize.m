function [w0,Ups0,v0,v0hat0,v1hat0,v0bar0,v1bar0,abar0,gamma0,btilde0,sys,chan,cons] = initialize(sys,chan)

K = sys.K; N = sys.N; Na = sys.Na; Nt = sys.Nt; pmax = sys.pmax; 
u = sys.u; hv = sys.hv; r = sys.r; D = sys.D; g0 = chan.g0;% G1 g2

init_loc = 2; init_w = 1; init_Ups = 1; % circular initalization

%% UAV location
if init_loc == 1 % initialize UAV right above the UE closest to the RIS
    dis_RIS2UE = vecnorm(repmat(r,1,K) - u(:,[1:K]));
    [~,k_min] = min(dis_RIS2UE);
    v0 = u(:,k_min) + [0;0;hv];
else % initialize UAV at the center of the area
    %v0 = r + [0;0;hmin - r(3)];
    hmin = sys.hmin;
    v0 = [D/2;D/2;hmin];
end



%% search for the best Ups and w
if init_w == 1
    w0 = zeros(Nt,K);
    for k = 1:K
        w0(:,k) = sqrt(pmax/K)*g0(:,k)/norm(g0(:,k));
        %w0(:,k) = pmax/K*ones(Nt,1);
    end
else
    Wrand = 500; % number of random precoding vector
    W0 = randn(K*Nt,Wrand)+ 1j*randn(K*Nt,Wrand);
    [~,minpos] = min(abs(diag(W0'*W0) - pmax));
    w0tmp = W0(:,minpos);
    w0 = zeros(Nt,K);
    for k = 1:K
        w0(:,k) = w0tmp(Nt*(k-1)+1:Nt*k);
    end
end

%% Ups0, Psi0
if init_Ups == 1
    Ups0 = 0*diag(exp(1i*2*pi.*randn(N,1)));
else
    Ups0 = zeros(N,N);
    gamma0 = zeros(K,1); rate_min0 = 0;% w = 100*randn(Nt,K) + 1j*randn(Nt,K);
    for m = 1:300
        Ups = diag(exp(1i*2*pi.*randn(N,1)));
        %w = randn(Nt,K) + 1j*randn(Nt,K);
        %w = pmax*w/norm(w,'fro');
        %norm(w,'fro')
        update_channel(v0,Ups)
        
        % compute SINR
        gamma = zeros(K,1); rate = zeros(K,1);
        for k = 1:K
            S_mk = abs(h(:,k)'*w0(:,k))^2;
            I_mk = sigma2(k);
            for j = 1:K
                if j ~= k
                    I_mk = I_mk + abs(h(:,k)'*w0(:,j))^2;
                end
            end
            gamma(k) = S_mk/I_mk;
            rate(k) = log(1 + gamma(k));
        end
        rate_min = min(rate);
        if rate_min > rate_min0
            Ups0 = Ups;
            %w0 = w;
            gamma0 = gamma;
            rate_min0 = rate_min;
        end
    end
end

Psi0 = zeros(N,N);
for n = 1:Na
    alpha_active = 1;
    Ups0(n,n) = sqrt(alpha_active)*Ups0(n,n);
    Psi0(n,n) = Ups0(n,n);
end

%% Initialize channels h0 h1 h2, effective channel h, and effective noise
[sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v0,Ups0,sys,chan);


%% Initialize Ntilde
btilde0 = zeros(K,1);
if N > 0
    [cons.Q,cons.Qtilde,cons.t,cons.ttilde,cons.e,cons.etilde,cons.Xi,btilde0] = update_Q(w0,Ups0,sys,chan);
end

%% Initialize constants c0,c1,c2,c3 and v0hat0,v1hat0,v0bar0,v1bar0,abar0
[cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups0,w0,sys,chan);
[v0hat0,v1hat0,v0bar0,v1bar0,abar0] = update_v(v0,sys,cons);

%% Initialize SINR values gamma
gamma0 = zeros(K,1);
for k = 1:K
    S_mk = abs(chan.h(:,k)'*w0(:,k))^2;
    I_mk = sys.sigma2(k);
    for j = 1:K
        if j ~= k
            I_mk = I_mk + abs(chan.h(:,k)'*w0(:,j))^2;
        end
    end
    gamma0(k) = S_mk/I_mk;
end


end % EOF