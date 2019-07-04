clc;
clear;
close all;
%����512����������
data_len=512;    %�����������еĳ���   
trials=100;      %���ʵ��Ĵ��� 
nvector=1:data_len;
a1=-0.975;
a2=0.95;
sigma_v_2=0.0731;
v=sqrt(sigma_v_2)*randn(data_len,1,trials);
u0=[0 0];
num=1;
den=[1 a1 a2];
Zi=filtic(num,den,u0);  %�˲����ĳ�ʼ����   
u=filter(num,den,v,Zi); %������������u(n)   
%������������
figure(1);
plot(u(:,:,1)) %��������һ����������
title('һ��ʵ���е���������u(n)')
xlabel('n')
ylabel('u(n)')
%LMS �����㷨
%��ʼ��
mu1=0.05; %����
mu2=0.005;%����
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
m=1;%ʵ���������
%LMS��������
for m=1:trials  
    for n=3:data_len-1 %u(2) u(1)����w(4)
        w1(:,n+1)=w1(:,n)...
            +mu1*u(n-1:-1:n-2,:,m)*conj(e1(n));%����Ԥ���� u(n)=[u(n-1) u(n-2)]
        w2(:,n+1)=w2(:,n)...
            +mu2*u(n-1:-1:n-2,:,m)*conj(e2(n));
        d1(n+1)=w1(:,n+1)'*u(n:-1:n-1,:,m);
        d2(n+1)=w2(:,n+1)'*u(n:-1:n-1,:,m);
        e1(n+1)=u(n+1,:,m)-d1(n+1);
        e2(n+1)=u(n+1,:,m)-d2(n+1);
        
        mean1(:,n+1)=mean1(:,n+1)+w1(:,n+1);%�������ƽ�����
        mean2(:,n+1)=mean2(:,n+1)+w2(:,n+1);
        
        mean_e1(n+1)=mean_e1(n+1)+e1(n+1)^2;%�������
        mean_e2(n+1)=mean_e2(n+1)+e2(n+1)^2;
        
    end
end
%�����˲���Ȩϵ��
%����0.05
figure(2);
plot(nvector,w1',nvector,mean1'/trials )
legend('����ʵ��w(1)','����ʵ��w(2)','100��ƽ��w(1)','100��ƽ��w(2)')
xlabel('n')
ylabel('Ȩϵ��')
%����0.005
figure(3);
plot(nvector,w2',nvector,mean2'/trials )
legend('����ʵ��w(1)','����ʵ��w(2)','100��ƽ��w(1)','100��ƽ��w(2)')
xlabel('n')
ylabel('Ȩϵ��')   
%����ʣ������ʧ������
%��ʼ��
wopt=zeros(2,trials);
Jmin=zeros(1,trials);
sum_eig=zeros(trials,1);
%ͨ��ά��-���򷽳̼�����С�������
for m=1:trials
    rm1=xcorr(u(:,:,m),'biased');
    R=[rm1(512),rm1(513);rm1(511),rm1(512)];
    p=[rm1(511);rm1(510)];
    wopt(:,m)=R\p;                %�������Ȩֵ
    [v,d]=eig(R);
    Jmin(m)=rm1(512)-p'*wopt(:,m);%����ά�����
    sum_eig(m)=d(1,1)+d(2,2);     %��������ֵ֮��
end
sJmin=sum(Jmin)/trials;           %100��ƽ�����
Jex1=mean_e1-sJmin;      %ʣ��������Jex1 
Jex2=mean_e2-sJmin;      %ʣ��������Jex2
sum_eig_100trials=sum(sum_eig)/100;
Jexfin=mu1*sJmin*(sum_eig_100trials...
    /(2-mu1*sum_eig_100trials));
Jexfin2=mu2*sJmin*(sum_eig_100trials...
    /(2-mu2*sum_eig_100trials));
M1=Jexfin/sJmin;     %ʧ������M1
M2=Jexfin2/sJmin;    %ʧ������M2
%ѧϰ����
figure(4);
plot(nvector,mean_e1/trials,nvector,mean_e2/trials )
legend('����0.05','����0.005')
xlabel('��������n')
ylabel('�������')
