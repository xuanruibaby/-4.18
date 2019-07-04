clear;
close all;
clc;
%������������ĸ�˹������ v(n)
N=3000;                       %�������г���
gv=0.0332;                    %��������
trials=100;                   %��������
v=randn(N,1,trials)*sqrt(gv); %������˹����������
%���ݸ�����ARģ�Ͳ���u(n)����
a1=-2;                %ģ����u(n-1)��ϵ��
a2=1.4;               %ģ����u(n-2)��ϵ��
a3=-0.4;              %ģ����u(n-3)��ϵ��
a4=0.0384;            %ģ����u(n-4)��ϵ��
u1=zeros(N,1,trials); 
for m=1:trials
    %����u(n)����
 for i=1:N-4
   u1(i+4,1,m)=-a1*u1(i+3,1,m)-a2*u1(i+2,1,m)-a3*u1(i+1,1,m)...
   -a4*u1(i,1,m)+v(i+4,1,m);
 end
end
%�������˲�
N2=2000;                   %�������˲�����
Jmin=0.005;                %������������
W_esti=zeros(4,N2,trials); %״̬������ʼ��
W=zeros(4,N2+1);
e=zeros(N2,1);
for m=1:trials
    P_esti=eye(4);%״̬�������ؾ����ʼ��
  for n=5:N2   %���ɹ۲����
    P_pre=P_esti;%Ԥ���������ؾ���
    U(:,n,m)=[u1(n-1,1,m);u1(n-2,1,m);u1(n-3,1,m);u1(n-4,1,m)];
    A=(U(:,n,m))'*P_pre*U(:,n,m)+Jmin;%��Ϣ��������ؾ���
    K=P_pre*U(:,n,m)/A;%����������
    alpha(n,m)=u1(n,1,m)-(U(:,n,m))'*W_esti(:,n,m);%��Ϣ����
    W_esti(:,n+1,m)=W_esti(:,n,m)+K*alpha(n,m);%״̬����
    P_esti=P_pre-K*(U(:,n,m))'*P_pre;%״̬�����������ؾ���
  end
   W=W+W_esti(:,:,m);
   e=e+abs(alpha(:,m)).^2;
end
W=1/trials*W;
e=1/trials*e;
h=1:N2+1;
figure(1);
plot(h,W);
xlabel('��������');
ylabel('Ȩֵ');
gtext('w1');
gtext('w2');
gtext('w3');
gtext('w4');
figure(2);
plot(1:N2,e);
xlabel('��������');
ylabel('MSE');
