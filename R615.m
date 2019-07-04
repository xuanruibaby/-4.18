clear;
close all;
clc;

a1=0.99;    %AR模型系数
sigma=0.995;%白噪声方差
N=1000;     %数据个数
trials=500; %实验次数
vn=sqrt(sigma)*randn(N,1,trials);%产生白噪声样本
nume=1;     %分子系数
deno=[1 a1];%分母系数
u0=zeros(length(deno)-1,1);%初始数据
n0=1;       %需实现n0步线性预测
M=2;        %滤波器阶数
L=N-n0;
w=zeros(M,L+1,trials);%存储权向量
epsilon=zeros(L,1,trials);
W_mean=zeros(M,L+1);
MSE=zeros(L,1);

for m=1:trials
    xic=filtic(nume,deno,u0);          %初始条件
    un=filter(nume,deno,vn(:,1,m),xic);%产生数据
    b=un(n0+1:N);                      %预测的期望响应
    un1=[zeros(M-1,1).',un.'].';       %扩展数据
    A=zeros(M,L);
    for k=1:L                          %构建观测数据矩阵
        A(:,k )= un1(M-1+k:-1:k);
    end
    delta=0.004;                       %调整参数
    lambda=0.98;                       %遗忘因子
    P1=eye(M)/delta;
    for k=1:L                          %RLS算法迭代过程
        PIn=P1*A(:,k);
        denok=lambda+A(:,k)'*PIn;
        kn=PIn/denok;
        epsilon(k,1,m)=b(k)-w(:,k,m)'*A(:,k);
        w(:,k+1,m)=w(:,k,m)+kn*conj(epsilon(k,1,m));
        P1=P1/lambda-kn*A(:,k)'*P1/lambda;
    end
    MSE=MSE+abs(epsilon(:,1,m)).^2;
    W_mean=W_mean+w(:,:,m);
end
figure(1);
plot(1:L,1/trials*MSE);
xlabel('迭代次数');
ylabel('MSE');
figure(2);
plot(1:L+1,1/trials*W_mean);
xlabel('迭代次数');
ylabel('权值');