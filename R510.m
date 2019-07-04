clc;
clear ��
close all;
%����ʵ������
data_len=10000;
sigma_v1_2=0.93627;
mu1=0.0005;%����
trials=100;%�������
u=zeros(data_len,1);
M=3;%�˲�������
mean1=zeros(M,data_len);%�������ƽ�����
e1=zeros(data_len,1);
mean_e1=e1;
d1=zeros(data_len,1);
for m=1:trials  
    v=sqrt(sigma_v1_2)*randn(data_len,1);
    %����ARģ���ź�u(n)
    a1=-0.99; %ARģ��ϵ��
    u(1)=v(1);
        for k=2:data_len
            u(k)=-a1*u(k-1)+v(k);
        end                
        w1=zeros(M,data_len);
    for n=(M+1):data_len-1 
        w1(:,n+1)=w1(:,n)+mu1*u(n-1:-1:n-M)*conj(e1(n));%����Ԥ���� u(n)=[u(n-1) u(n-2)]
        
        d1(n+1)=w1(:,n+1)'*u(n:-1:n-M+1);
        
        e1(n+1)=u(n+1)-d1(n+1);        
        mean1(:,n+1)=mean1(:,n+1)+w1(:,n+1);%�������ƽ�����        
        mean_e1(n+1)=mean_e1(n+1)+e1(n+1)^2;%�������        
    end
end
nvector=1:data_len;
figure(1)
plot(nvector,w1',nvector,mean1'/trials )
% legend('����ʵ��w(1)','����ʵ��w(2)','100��ƽ��w(1)','100��ƽ��w(2)')
xlabel('��������n')
ylabel('Ȩϵ��')
%ѧϰ����
figure(2);
plot(nvector,mean_e1/trials)
legend('����0.0005')
xlabel('��������n')
ylabel('�������')
