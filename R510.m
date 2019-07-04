clc;
clear ；
close all;
%产生实验数据
data_len=10000;
sigma_v1_2=0.93627;
mu1=0.0005;%步长
trials=100;%试验次数
u=zeros(data_len,1);
M=3;%滤波器阶数
mean1=zeros(M,data_len);%多次试验平均结果
e1=zeros(data_len,1);
mean_e1=e1;
d1=zeros(data_len,1);
for m=1:trials  
    v=sqrt(sigma_v1_2)*randn(data_len,1);
    %生成AR模型信号u(n)
    a1=-0.99; %AR模型系数
    u(1)=v(1);
        for k=2:data_len
            u(k)=-a1*u(k-1)+v(k);
        end                
        w1=zeros(M,data_len);
    for n=(M+1):data_len-1 
        w1(:,n+1)=w1(:,n)+mu1*u(n-1:-1:n-M)*conj(e1(n));%线性预测器 u(n)=[u(n-1) u(n-2)]
        
        d1(n+1)=w1(:,n+1)'*u(n:-1:n-M+1);
        
        e1(n+1)=u(n+1)-d1(n+1);        
        mean1(:,n+1)=mean1(:,n+1)+w1(:,n+1);%多次试验平均结果        
        mean_e1(n+1)=mean_e1(n+1)+e1(n+1)^2;%均方误差        
    end
end
nvector=1:data_len;
figure(1)
plot(nvector,w1',nvector,mean1'/trials )
% legend('单次实验w(1)','单次实验w(2)','100次平均w(1)','100次平均w(2)')
xlabel('迭代次数n')
ylabel('权系数')
%学习曲线
figure(2);
plot(nvector,mean_e1/trials)
legend('步长0.0005')
xlabel('迭代次数n')
ylabel('均方误差')
