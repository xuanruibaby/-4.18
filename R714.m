clear;
close all;
clc;
%产生给定方差的高斯白噪声 v(n)
N=3000;                       %噪声序列长度
gv=0.0332;                    %噪声方差
trials=100;                   %迭代次数
v=randn(N,1,trials)*sqrt(gv); %产生高斯白噪声序列
%根据给定的AR模型产生u(n)序列
a1=-2;                %模型中u(n-1)的系数
a2=1.4;               %模型中u(n-2)的系数
a3=-0.4;              %模型中u(n-3)的系数
a4=0.0384;            %模型中u(n-4)的系数
u1=zeros(N,1,trials); 
for m=1:trials
    %产生u(n)序列
 for i=1:N-4
   u1(i+4,1,m)=-a1*u1(i+3,1,m)-a2*u1(i+2,1,m)-a3*u1(i+1,1,m)...
   -a4*u1(i,1,m)+v(i+4,1,m);
 end
end
%卡尔曼滤波
N2=2000;                   %卡尔曼滤波点数
Jmin=0.005;                %测量噪声方差
W_esti=zeros(4,N2,trials); %状态向量初始化
W=zeros(4,N2+1);
e=zeros(N2,1);
for m=1:trials
    P_esti=eye(4);%状态误差自相关矩阵初始化
  for n=5:N2   %构成观测矩阵
    P_pre=P_esti;%预测误差自相关矩阵
    U(:,n,m)=[u1(n-1,1,m);u1(n-2,1,m);u1(n-3,1,m);u1(n-4,1,m)];
    A=(U(:,n,m))'*P_pre*U(:,n,m)+Jmin;%新息过程自相关矩阵
    K=P_pre*U(:,n,m)/A;%卡尔曼增益
    alpha(n,m)=u1(n,1,m)-(U(:,n,m))'*W_esti(:,n,m);%新息过程
    W_esti(:,n+1,m)=W_esti(:,n,m)+K*alpha(n,m);%状态估计
    P_esti=P_pre-K*(U(:,n,m))'*P_pre;%状态估计误差自相关矩阵
  end
   W=W+W_esti(:,:,m);
   e=e+abs(alpha(:,m)).^2;
end
W=1/trials*W;
e=1/trials*e;
h=1:N2+1;
figure(1);
plot(h,W);
xlabel('迭代次数');
ylabel('权值');
gtext('w1');
gtext('w2');
gtext('w3');
gtext('w4');
figure(2);
plot(1:N2,e);
xlabel('迭代次数');
ylabel('MSE');
