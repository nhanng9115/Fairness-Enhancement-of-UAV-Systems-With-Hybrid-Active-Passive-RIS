function [Q,Qtilde,t,ttilde,e,etilde,Xi,btilde0] = update_Q(w0,Ups0,sys,chan)


K = sys.K; N = sys.N; Na = sys.Na; ONE = sys.ONE;
sigma2_r = sys.sigma2_r; sigma2_u = sys.sigma2_r; 
h0 = chan.h0; H1 = chan.H1; h2 = chan.h2;

test_derivation = 0;

%% compute H2tilde and h12tilde
H2tilde = zeros(N,N,K); h12tilde = zeros(N,K,K); h0bar = zeros(K,K); h1bar = zeros(N,K);

for k = 1:K
    h1bar(:,k) = H1*w0(:,k);
    H2tilde(:,:,k) = diag(h2(:,k)');
end

for k = 1:K
    for j = 1:K
        h0bar(k,j) = h0(:,k)'*w0(:,j);
        h12tilde(:,k,j) = H2tilde(:,:,k)*h1bar(:,j);
    end
end

%% compute Q,q,c4 for all m, k
Q = zeros(N,N,K); t = zeros(N,K);
Qtilde = zeros(N,N,K); ttilde = zeros(N,K);
e = zeros(K,1); etilde = zeros(K,1);
Ntilde0 = zeros(K,1);

alpha0 = diag(Ups0);
scale = sigma2_u; %%%%%%%%%%%%%%%%%%%%%%%%% chỗ này em phải chia cả tử số và mẫu số của SINR ở (39) cho sigma2_u, k thì nó k chạy được ????
for k = 1:K
    %scale = norm(w0(:,k))^2;
    %% compute Q and q
    e(k) = abs(h0bar(k,k))^2/scale;
    Q(:,:,k) = conj(h12tilde(:,k,k))*h12tilde(:,k,k).'/scale;
    t(:,k) = conj(h12tilde(:,k,k))*h0bar(k,k)/scale;
    
    %% compute Qtilde, qtilde
    Qtilde(:,:,k) = sigma2_r/scale * ONE*conj(H2tilde(:,:,k))*H2tilde(:,:,k).'*ONE;
    ttilde(:,k) = 0;
    etilde(k) = sigma2_u/scale;
    
    for j = 1:K
        if j ~= k
            etilde(k) = etilde(k)  + abs(h0bar(k,j))^2/scale;
            Qtilde(:,:,k) = Qtilde(:,:,k) + conj(h12tilde(:,k,j))*h12tilde(:,k,j).'/scale;
            ttilde(:,k) = ttilde(:,k) + conj(h12tilde(:,k,j))*h0bar(k,j)/scale;
        end
    end
    btilde0(k) = alpha0'*Qtilde(:,:,k)*alpha0 + 2*real(alpha0'*ttilde(:,k));
    
    % compute Ntilde
    Ntilde0(k) = alpha0'*Q(:,:,k)*alpha0 + 2*real(alpha0'*t(:,k)) + e(k);

    % compute gammatilde
    %num = alpha0'*Q(:,:,k)*alpha0 + 2*real(alpha0'*q(:,k)) + c4(k);
    %denom = alpha0'*Qtilde(:,:,k)*alpha0 + 2*real(alpha0'*qtilde(:,k)) + c5(k);
    %gammatilde0(k) = num/denom;
end

%% compute Xi
Xi = zeros(N,N);
if Na > 0
    for nn = 1:Na
        xi_n = sigma2_r + norm(H1(nn,:))^2 * norm(w0,'fro')^2;
        Xi(nn,nn) = xi_n;
    end
end

%% ************** phần dưới này chỉ để check cái form (39) thôi chứ không liên quan đến optimization đâu anh.
%% ==================================================================================
%% check the derivations in (68) and (69)
if test_derivation == 1
    Psi0 = Ups0.*ONE;
    alpha = diag(Ups0);
    for k = 1:K
        
        % test numerator of SINR ==> ok
        Num = abs( h0(:,k)'*w0(:,k) + h2(:,k)'*Ups0*H1*w0(:,k) )^2
        Num = abs( h0bar(k,k) + h2(:,k)'*Ups0*h1bar(:,k) )^2
        Num1 = abs( h0bar(k,k) + alpha.'*h12tilde(:,k,k) )^2
        Num2 = real(abs(h0bar(k,k))^2 + alpha'*Q(:,:,k)*alpha + alpha'*t(:,k) + t(:,k)'*alpha)
        Num_ana = real(alpha'*Q(:,:,k)*alpha + 2*real(alpha'*t(:,k)) + e(k))
        
        % test denominator of SINR ==>
        A0 = sigma2_r*norm(h2(:,k)'*Psi0)^2;
        Denom0 = sigma2_u + A0;
        Denom1 = sigma2_u + A0;
        Denom2 = sigma2_u + A0;
        for j = 1:K
            if j ~= k
                Denom0 = Denom0 + abs(h0(:,k)'*w0(:,j) + h2(:,k)'*Ups0*H1*w0(:,j))^2;
                Denom1 = Denom1 + abs(h0bar(k,j) + h2(:,k)'*Ups0*h1bar(:,j))^2;
                Denom2 = Denom2 + abs(h0bar(k,j) + alpha.'*h12tilde(:,k,j))^2;
            end
        end
        Denom0
        Denom1
        Denom2
        Denom_ana = real(etilde(k) + alpha'*Qtilde(:,:,k)*alpha + 2*real(alpha'*ttilde(:,k)))
    end
end

end % EOF