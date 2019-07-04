clear;
close all;
clc;
%产生阵列接收信号x(n)
N=100; %信号样本长度
M=10;  %阵元个数
K=2;   %信源个数
theta=[-10;40]*pi/180;   %信源入射角度（单位：度）
SNR=[10;20];             %信噪比（单位：dB)
Am=[sqrt(10.^(SNR/10))];
S=Am*ones(1,N);
S(2,:)=S(2,:).*exp(1i*2*pi*rand(1,N)); %加入随机相位
for a=1:M
    for b=1:K
        A(a,b)=exp(-1i*(a-1)*pi*sin(theta(b)));
    end
end
V=zeros(M,N); %产生方向矩阵
for m=1:M
    v=wgn(1,N,0,'complex'); %产生复高斯噪声
    v=v-mean(v);  %零均值
    v=v/std(v);   %方差归一
    V(m,:)=v;
end
X=A*S+V;          %阵列接收信号
%利用接收数据估计信号的空间相关矩阵R
R=zeros(M,M);  
for i=1:N
    R=R+X(:,i)*X(:,i)';
end
R=R/N;
%MUSIC算法
[VR,D]=eig(R);       %特征值分解
[B,IX]=sort(diag(D));
G=VR(:,IX(M-K:-1:1));%得到矩阵G
P=[];
for n=-pi/2:pi/180:pi/2
    a=exp(-1i*[0:M-1]'*pi*sin(n));%方向向量
    P=[P,1/(a'*G*G'*a)];          %MUSIC谱
end
t=-90:1:90;
figure(1);
plot(t,10*log10(abs(P)/max(abs(P))))
title('MUSIC算法')
ylabel('归一化功率谱/dB')
xlabel('空间角度')
%RootMUSIC算法
syms z
pz=z.^([0:M-1]');
pz1=(z^(-1)).^([0:M-1]);
fz=z^(M-1)*pz1*G*G'*pz; %构造多项式
a=sym2poly(fz);
r=roots(a);             %求根
r1=abs(r);
for i=1:2*K
    [Y,I(i)]=min(abs(r1-1));
    r1(I(i))=inf;
end
%确定靠近单位圆的2K个解
for i=1:2*K   
    theta_esti1(i)=asin(-angle(r(I(i)))/pi)*180/pi;
end
display('按照RootMUSIC算法进行100次蒙特卡洛实验得到的平均结果为:');
sort(theta_esti1)

%ESPRIT算法
S=VR(:,IX(M:-1:M-K+1));%构造信号子空间
S1=S(1:M-1,:);
S2=S(2:M,:);
fai=S1\S2;
[U_fai,V_fai]=eig(fai);%特征值分解
for i=1:K
    theta_esti2(i)=asin(-angle(V_fai(i,i))/pi)*180/pi;
end
display('按照ESPRIT算法进行100次蒙特卡洛实验得到的平均结果为:');
sort(theta_esti2)
%MVDR算法 
P=[];
for n=-pi/2:pi/180:pi/2
    a=exp(-1i*[0:M-1]'*pi*sin(n));%方向向量
    P=[P,1/(a'*inv(R)*a)];        %MVDR谱
end
P=abs(P/max(abs(P)));
P=10*log10(P);
figure(2);
plot(t,P);
title('MVDR算法')
ylabel('归一化功率谱/dB')
xlabel('空间角度')
%F-SAPES算法 
P=6;        %子阵数目
L=M+1-P;    %子阵阵元数目
Rf=zeros(L,L);
for i=1:P
    Rf=Rf+X(i:i+L-1,:)*X(i:i+L-1,:)'/N;
end
Rf=Rf/P;    %子阵平滑后的空间相关矩阵
n1=0:P-1;
n2=0:L-1;
cc=[1 zeros(1,L-1)];
for n3=-90:.5:90
    fy=exp(1i*pi*sin(n3/180*pi));
    tt=[(fy.^(n1')).' zeros(1,M-P)];
    Tfy=toeplitz(cc,tt);        %构造toeplitz矩阵
    GfTheta=1./(P^2)*Tfy*R*Tfy';%获得GF矩阵
    Qf=Rf-Rf-GfTheta;           %获得Qf矩阵
    aTheta=fy.^(-n2');
    Wof=((inv(Qf))*aTheta)./(aTheta'*inv(Qf)*aTheta);%计算最优权向量
    sigma2Theta(((n3+90)/.5+1))=Wof'*GfTheta*Wof;
end
sigma2Theta=abs(sigma2Theta/max(abs(sigma2Theta)));
sigma2Theta=10*log10(sigma2Theta);
figure(3);
t1=-90:.5:90;
plot(t1,sigma2Theta);
title('F-SAPES算法')
ylabel('归一化功率谱/dB')
xlabel('空间角度')