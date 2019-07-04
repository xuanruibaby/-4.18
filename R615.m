clear;
close all;
clc;

a1=0.99;    %ARģ��ϵ��
sigma=0.995;%����������
N=1000;     %���ݸ���
trials=500; %ʵ�����
vn=sqrt(sigma)*randn(N,1,trials);%��������������
nume=1;     %����ϵ��
deno=[1 a1];%��ĸϵ��
u0=zeros(length(deno)-1,1);%��ʼ����
n0=1;       %��ʵ��n0������Ԥ��
M=2;        %�˲�������
L=N-n0;
w=zeros(M,L+1,trials);%�洢Ȩ����
epsilon=zeros(L,1,trials);
W_mean=zeros(M,L+1);
MSE=zeros(L,1);

for m=1:trials
    xic=filtic(nume,deno,u0);          %��ʼ����
    un=filter(nume,deno,vn(:,1,m),xic);%��������
    b=un(n0+1:N);                      %Ԥ���������Ӧ
    un1=[zeros(M-1,1).',un.'].';       %��չ����
    A=zeros(M,L);
    for k=1:L                          %�����۲����ݾ���
        A(:,k )= un1(M-1+k:-1:k);
    end
    delta=0.004;                       %��������
    lambda=0.98;                       %��������
    P1=eye(M)/delta;
    for k=1:L                          %RLS�㷨��������
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
xlabel('��������');
ylabel('MSE');
figure(2);
plot(1:L+1,1/trials*W_mean);
xlabel('��������');
ylabel('Ȩֵ');