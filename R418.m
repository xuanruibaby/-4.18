clc;
clear;
close all;
%产生512点样本序列
data_len=512;    %设置样本序列的长度   
trials=100;      %随机实验的次数 
nvector=1:data_len;
a1=-0.975;
a2=0.95;
sigma_v_2=0.0731;
v=sqrt(sigma_v_2)*randn(data_len,1,trials);
u0=[0 0];
num=1;
den=[1 a1 a2];
Zi=filtic(num,den,u0);  %滤波器的初始条件   
u=filter(num,den,v,Zi); %产生样本序列u(n)   
%画出样本序列
figure(1);
plot(u(:,:,1)) %画出其中一个样本序列
title('一次实验中的样本序列u(n)')
xlabel('n')
ylabel('u(n)')
%LMS 迭代算法
%初始化
mu1=0.05; %步长
mu2=0.005;%步长
w1=zeros(2,data_len);%2*512*100
w2=zeros(2,data_len);
mean1=w1;
mean2=w2;
e1=zeros(data_len,1);
e2=zeros(data_len,1);
mean_e1=e1;
mean_e2=e2;
d1=zeros(data_len,1);
d2=zeros(data_len,1);
m=1;%实验次数控制
%LMS迭代过程
for m=1:trials  
    for n=3:data_len-1 %u(2) u(1)估计w(4)
        w1(:,n+1)=w1(:,n)...
            +mu1*u(n-1:-1:n-2,:,m)*conj(e1(n));%线性预测器 u(n)=[u(n-1) u(n-2)]
        w2(:,n+1)=w2(:,n)...
            +mu2*u(n-1:-1:n-2,:,m)*conj(e2(n));
        d1(n+1)=w1(:,n+1)'*u(n:-1:n-1,:,m);
        d2(n+1)=w2(:,n+1)'*u(n:-1:n-1,:,m);
        e1(n+1)=u(n+1,:,m)-d1(n+1);
        e2(n+1)=u(n+1,:,m)-d2(n+1);
        
        mean1(:,n+1)=mean1(:,n+1)+w1(:,n+1);%多次试验平均结果
        mean2(:,n+1)=mean2(:,n+1)+w2(:,n+1);
        
        mean_e1(n+1)=mean_e1(n+1)+e1(n+1)^2;%均方误差
        mean_e2(n+1)=mean_e2(n+1)+e2(n+1)^2;
        
    end
end
%画出滤波器权系数
%步长0.05
figure(2);
plot(nvector,w1',nvector,mean1'/trials )
legend('单次实验w(1)','单次实验w(2)','100次平均w(1)','100次平均w(2)')
xlabel('n')
ylabel('权系数')
%步长0.005
figure(3);
plot(nvector,w2',nvector,mean2'/trials )
legend('单次实验w(1)','单次实验w(2)','100次平均w(1)','100次平均w(2)')
xlabel('n')
ylabel('权系数')   
%计算剩余误差和失调参数
%初始化
wopt=zeros(2,trials);
Jmin=zeros(1,trials);
sum_eig=zeros(trials,1);
%通过维纳-霍夫方程计算最小均方误差
for m=1:trials
    rm1=xcorr(u(:,:,m),'biased');
    R=[rm1(512),rm1(513);rm1(511),rm1(512)];
    p=[rm1(511);rm1(510)];
    wopt(:,m)=R\p;                %单次最佳权值
    [v,d]=eig(R);
    Jmin(m)=rm1(512)-p'*wopt(:,m);%单次维纳误差
    sum_eig(m)=d(1,1)+d(2,2);     %单次特征值之和
end
sJmin=sum(Jmin)/trials;           %100次平均误差
Jex1=mean_e1-sJmin;      %剩余均方误差Jex1 
Jex2=mean_e2-sJmin;      %剩余均方误差Jex2
sum_eig_100trials=sum(sum_eig)/100;
Jexfin=mu1*sJmin*(sum_eig_100trials...
    /(2-mu1*sum_eig_100trials));
Jexfin2=mu2*sJmin*(sum_eig_100trials...
    /(2-mu2*sum_eig_100trials));
M1=Jexfin/sJmin;     %失调参数M1
M2=Jexfin2/sJmin;    %失调参数M2
%学习曲线
figure(4);
plot(nvector,mean_e1/trials,nvector,mean_e2/trials )
legend('步长0.05','步长0.005')
xlabel('迭代次数n')
ylabel('均方误差')
