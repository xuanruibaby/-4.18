clear;
close all;
clc;

%产生观测信号Un
M=4;                        %滤波器抽头数
N=1000;                     %样本数
f=[0.1 0.25 0.27];          %归一化频率
SNR=[30 30 27];             %信噪比
sigma=1;                    %方差
Am=sqrt(sigma*10.^(SNR/10));%信号幅度
t=linspace(0,1,N);
phi=2*pi*rand(size(f));     %随机相位
vn=sqrt(sigma/2)*rand(size(t))+1i*sqrt(sigma/2)*rand(size(t));%加性高斯白噪声
Un=vn;
for k=1:length(f)
    s=Am(k)*exp(1i*2*pi*N*f(k).*t+1i*phi(k));
    Un=Un+s;
end
Un=Un.';

%构建矩阵
A=zeros(M,N-M+1);%构建观测矩阵
for n=1:N-M+1
    A(:,n)=Un(M+n-1:-1:n);
end
[U,S,V]=svd(A');
invphi=V*inv(S'*S)*V';%构建矩阵phi

%构建向量
p=1024;
f=linspace(-0.5,0.5,p);
omega=2*pi*f;
a=zeros(M,p);%每一列均为一个向量
for k=1:p
    for m=1:M
        a(m,k)=exp(-1i*omega(k)*(m-1));
    end
end

%计算MVDR谱
Pmvdr=zeros(size(omega));
for k=1:p
    Pmvdr(k)=1/(a(:,k)'*invphi*a(:,k));
end
%Pmvdr=abs(Pmvdr);
Pmvdr=abs(Pmvdr/max(abs(Pmvdr)));
Pmvdr=10*log10(Pmvdr);

h=1:p;
h=(h-512)/1024;
plot(h,Pmvdr);
xlabel('M=4的MVDR谱');
ylabel('归一化MVDR谱');