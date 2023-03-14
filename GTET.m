function [Wx,TFx,GD,t,f,C] = GTET(x,fs,s,k,gamma)
% % Example1：
% % Parameters
% N = 700;
% fs = 200;
% t = (0:N-1)/fs;
% f = (0:N/2)*fs/N;  
% % Signal
% a = 2;
% b = 2;
% c = 0.015;
% d = 0.01;
% e = 0.1;
% GD_t4 = a + b*exp(-c*f+d).*(c*cos(e*pi*f)+e*pi*sin(e*pi*f)); 
% Phi_f4 = a*f-b*exp(-c*f+d).*cos(e*pi*f);
% A_f4 = 1.5.*ones(size(f));
% X4 = A_f4.*exp(-1i*2*pi*Phi_f4);
% X4(end) = -A_f4(end);%
% Y4 = [X4  conj(fliplr(X4(2:end-1)))];    
% y = ifft(Y4);
% figure;
% plot(t,real(y),'color','b','linewidth',0.2);
% xlabel({'Time(s)','(b)'},'FontSize',10);
% ylabel('Amplitude','FontSize',10);
% % TFA
% k=3;
% s = 0.12;
% gamma = 0;
% [Wx,TFx,GD,t,f,C] = GTET(y,fs,s,k,gamma);
% figure('color',[1 1 1]);
% imagesc(t,f,abs(Wx'));
% axis xy
% xlabel('Time(s)','FontSize',8);
% ylabel('Frequency(Hz)','FontSize',8);
% set(gca,'FontSize',8);
% figure('color',[1 1 1]);
% imagesc(t,f,abs(TFx'));
% axis xy
% xlabel('Time(s)','FontSize',8);
% ylabel('Frequency(Hz)','FontSize',8);
%%%%%%%
% If you use this code, please refer to the following paper:
% 1.Generalized Transient-Extracting Transform and Its Accurate Signal
% Reconstruction;2.Generalized Transient-Squeezing Transform: Algorithm and Applications
% https://www.researchgate.net/profile/Wenjie-Bao/research
% I hope this code can help you with your research work.
% I wish you all the best in your research!
% Welcome to communicate with me through the following ways:
% Email: baowenjie@sjtu.edu.cn; baowenjie_mail@163.com;
% WeChat: baowenjie0001
% Copyright (c) 2021 BAO WENJIE. All rights reserved.
%% 参数数目检查
if (nargin > 5)
    error('输入参数过多！');
elseif(nargin == 4)    
    gamma =0.000;
elseif(nargin < 4)
    error('缺少输入参数！');
end
%%
[xrow,~] = size(x);
if (xrow~=1)
    x = x';
end
N = length(x);
t = (0:N-1)/fs;
tao = (0:N-1)/fs;
if mod(N,2)==0
    L = N/2+1;
else
    L = (N+1)/2;
end
f = fs/N*(0:L-1);
TFx = zeros(N,L);
dt = 1/fs;
df = f(2)-f(1);
%     a = (pi*s^2)^(-1/4);
%     b = -1/s^2;
%     c = sqrt(2*s)*pi^(1/4);
    d = -s^2;
%%
if k == 1
    S1 = S_wkG(s,L,N,t,tao,x,1);
    S2 = S_wkG(s,L,N,t,tao,x,2);
    Wx = S1;
    Denominator = S1;
    Numerator = d*S2;
    
    C = cell(1,2*k-1);
    for i = 1:2*k-1
        C{1,i} = S_wkG(s,L,N,t,tao,x,i);
    end    
else
    C = cell(1,2*k-1);
    for i = 1:2*k-1
        C{1,i} = S_wkG(s,L,N,t,tao,x,i);
    end
    Wx = C{1,1};
    sum_A_1k = zeros(size(Wx));
    A_det = zeros(size(Wx));
    C_n = STRU_CELL(C,k);
    C_1 = C_n;
    C_nn = C_1(1,:);
    C_1(1,:) = [];
    for i = 1:k
        C_2 = C_1;
        C_2(:,i) = [];
        C_2_det = DET_CELL_D(C_2);
        if i>1
        sum_A_1k = sum_A_1k + (-1)^(i+1)*(i-1)*C_2_det.*C_n{1,i-1};
        end
        A_det = A_det + (-1)^(i+1)*C_2_det.*C_n{1,i};
    end
    Denominator = A_det;
    Numerator = sum_A_1k;
end
    p = Numerator./Denominator;
    for ptr = 1:N
        p(ptr,:) = p(ptr,:) - 1i*t(ptr);
    end
    GD = -imag(p);
    GD(abs(Denominator)<gamma)=0;
    %重排
    for k = 1:N
        for m = 1:L
            time = min(max(1 + round((real(GD(k,m))-t(1))/dt),1),N);
            if (time == k)
                TFx(time, m) = Wx(k,m);
            end
        end
    end
end
function [S] = S_wkG(s,L,N,t,tao,x,n)
S = zeros(N,L);
for ptr = 1:N
    gh = (1j)^(n-1)*fg_k(t-tao(ptr),n-1,s);
    gh = conj(gh);
    xcpsi = fft(gh .* x);
    S(ptr,:) = xcpsi(1:L);
end 
end
function [G] = fg_k(w,k,s)
a = (pi*s^2)^(-1/4);
b = -1/s^2;
if k == 0 
    G = a.*exp(b*w.^2/2);
elseif k == 1 
    G = b.*w.*a.*exp(b*w.^2/2);
else
    G = b*((k-1)*fg_k(w,k-2,s)+w.*fg_k(w,k-1,s));
end
end
function [C_N] = STRU_CELL(C,N)
    l = length(C);
    if (2*N-1)>l
        error('cell长度不够！');
    end
    C_N = cell(N);
    
    for i = 1:N
       for j = 1:N
          C_N{i,j}=C{i+j-1}; 
       end
    end

end
function [C_det] = DET_CELL_D(C)
N = length(C);
if N>=2
    C_1 = C;
    C_1(1,:)=[];
    sum = zeros(size(C{1,1}));
    for i=1:N
        C_2 = C_1;
        C_2(:,i) = [];
        sum = sum + (-1)^(i+1)*C{1,i}.*DET_CELL_D(C_2); 
    end
    C_det = sum;
else 
  C_det = C{1,1};      
end      
end