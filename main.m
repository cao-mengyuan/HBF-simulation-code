clc
clear

global N N_RF K P_t
global H H_t P
global V_RF

K = 8;
N = 64;
N_RF = 9;

P_t = K;
beta = ones(1,K);% beta_k = 1; %均匀分配权重
sigma2 = K;

SNR = P_t/sigma2;% SNR = P/sigma^2;
SNR_dB = 10* log10(SNR);

H = channel(N,K); % 生成信道

% 生成随机可行解
% 从可行的一个解开始
P = eye(K);
% V_RF = ones(N,N_RF);
tt = 2*pi*rand(1,N*N_RF);
ttt = exp(1j*tt); 
V_RF = reshape(ttt,N,N_RF);
% 生成随机数，才能保证V_RF满秩，否则后面的A_j不满秩

H_t = P^(-0.5) * H; 

% 预分配空间
for loop = 1:200
    V_RF_last = V_RF;
    V_RF = change_V_RF();
%     flag = 1;
%     for jj = 1:1:N_RF
%         for ii = 1:1:N
%             ttt_r = real(V_RF(ii,jj))/real(V_RF_last(ii,jj))
%             ttt_i = imag(V_RF(ii,jj))/imag(V_RF_last(ii,jj))
%             if(ttt_r >= 1.05 || ttt_r <= 0.95 || ttt_i <= 0.95 || ttt_i >= 1.05) 
%                 flag = 0; 
%                 break;
%             end
%         end
%         if(flag == 0) break; end
%     end
%     if(flag == 1) 
%         disp("successfully converged"); 
%         V_RF
%         break; 
%     end
end

max_ = 0;
for jj = 1:N_RF
    for ii = 1:N
        ttt_r = abs(real(V_RF(ii,jj))/real(V_RF_last(ii,jj))-1);
        ttt_i = abs(imag(V_RF(ii,jj))/imag(V_RF_last(ii,jj))-1);
        max_ = max([ttt_r,ttt_i,max_]);
    end
end
% 生成功率分配矩阵
V_D_t = (V_RF') * (H')/ ( H * V_RF * (V_RF') * (H'));
Q_t = (V_D_t') * (V_RF') * V_RF * V_D_t;

% 迭代求出lamda
lamda = 1;
while 1
    initPower = 0;
    posi = 0;
    for k = 1:1:K
        tttt = (beta(k)/lamda) - Q_t(k,k)*sigma2;
        if(tttt > 0 )
            initPower = initPower + tttt;
            posi = posi + 1;
        end
    end
    if( abs(initPower / P_t -1) <= 0.05 )
        disp("find P");
        break;
    end
    if(posi > 0) lamda = lamda + 0.5*(initPower - P_t)/posi;
    else lamda = lamda/4;
    end
end
% 求出P
P = zeros(K,K);
for kk = 1:1:K
    P(kk,kk) = max([(beta(kk)/lamda) - Q_t(kk,kk)*sigma2, 0]) / Q_t(kk,kk);
end

% V_D
V_D = (V_RF') * (H') / ( H *V_RF * (V_RF') *(H'));

% 计算 R_k 求和
R = zeros(1,K);
for k = 1:1:K
    R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
end
Sum_R = sum(R)