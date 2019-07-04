clear;
close all;
clc;
%�������н����ź�x(n)
N=100; %�ź���������
M=10;  %��Ԫ����
K=2;   %��Դ����
theta=[-10;40]*pi/180;   %��Դ����Ƕȣ���λ���ȣ�
SNR=[10;20];             %����ȣ���λ��dB)
Am=[sqrt(10.^(SNR/10))];
S=Am*ones(1,N);
S(2,:)=S(2,:).*exp(1i*2*pi*rand(1,N)); %���������λ
for a=1:M
    for b=1:K
        A(a,b)=exp(-1i*(a-1)*pi*sin(theta(b)));
    end
end
V=zeros(M,N); %�����������
for m=1:M
    v=wgn(1,N,0,'complex'); %��������˹����
    v=v-mean(v);  %���ֵ
    v=v/std(v);   %�����һ
    V(m,:)=v;
end
X=A*S+V;          %���н����ź�
%���ý������ݹ����źŵĿռ���ؾ���R
R=zeros(M,M);  
for i=1:N
    R=R+X(:,i)*X(:,i)';
end
R=R/N;
%MUSIC�㷨
[VR,D]=eig(R);       %����ֵ�ֽ�
[B,IX]=sort(diag(D));
G=VR(:,IX(M-K:-1:1));%�õ�����G
P=[];
for n=-pi/2:pi/180:pi/2
    a=exp(-1i*[0:M-1]'*pi*sin(n));%��������
    P=[P,1/(a'*G*G'*a)];          %MUSIC��
end
t=-90:1:90;
figure(1);
plot(t,10*log10(abs(P)/max(abs(P))))
title('MUSIC�㷨')
ylabel('��һ��������/dB')
xlabel('�ռ�Ƕ�')
%RootMUSIC�㷨
syms z
pz=z.^([0:M-1]');
pz1=(z^(-1)).^([0:M-1]);
fz=z^(M-1)*pz1*G*G'*pz; %�������ʽ
a=sym2poly(fz);
r=roots(a);             %���
r1=abs(r);
for i=1:2*K
    [Y,I(i)]=min(abs(r1-1));
    r1(I(i))=inf;
end
%ȷ��������λԲ��2K����
for i=1:2*K   
    theta_esti1(i)=asin(-angle(r(I(i)))/pi)*180/pi;
end
display('����RootMUSIC�㷨����100�����ؿ���ʵ��õ���ƽ�����Ϊ:');
sort(theta_esti1)

%ESPRIT�㷨
S=VR(:,IX(M:-1:M-K+1));%�����ź��ӿռ�
S1=S(1:M-1,:);
S2=S(2:M,:);
fai=S1\S2;
[U_fai,V_fai]=eig(fai);%����ֵ�ֽ�
for i=1:K
    theta_esti2(i)=asin(-angle(V_fai(i,i))/pi)*180/pi;
end
display('����ESPRIT�㷨����100�����ؿ���ʵ��õ���ƽ�����Ϊ:');
sort(theta_esti2)
%MVDR�㷨 
P=[];
for n=-pi/2:pi/180:pi/2
    a=exp(-1i*[0:M-1]'*pi*sin(n));%��������
    P=[P,1/(a'*inv(R)*a)];        %MVDR��
end
P=abs(P/max(abs(P)));
P=10*log10(P);
figure(2);
plot(t,P);
title('MVDR�㷨')
ylabel('��һ��������/dB')
xlabel('�ռ�Ƕ�')
%F-SAPES�㷨 
P=6;        %������Ŀ
L=M+1-P;    %������Ԫ��Ŀ
Rf=zeros(L,L);
for i=1:P
    Rf=Rf+X(i:i+L-1,:)*X(i:i+L-1,:)'/N;
end
Rf=Rf/P;    %����ƽ����Ŀռ���ؾ���
n1=0:P-1;
n2=0:L-1;
cc=[1 zeros(1,L-1)];
for n3=-90:.5:90
    fy=exp(1i*pi*sin(n3/180*pi));
    tt=[(fy.^(n1')).' zeros(1,M-P)];
    Tfy=toeplitz(cc,tt);        %����toeplitz����
    GfTheta=1./(P^2)*Tfy*R*Tfy';%���GF����
    Qf=Rf-Rf-GfTheta;           %���Qf����
    aTheta=fy.^(-n2');
    Wof=((inv(Qf))*aTheta)./(aTheta'*inv(Qf)*aTheta);%��������Ȩ����
    sigma2Theta(((n3+90)/.5+1))=Wof'*GfTheta*Wof;
end
sigma2Theta=abs(sigma2Theta/max(abs(sigma2Theta)));
sigma2Theta=10*log10(sigma2Theta);
figure(3);
t1=-90:.5:90;
plot(t1,sigma2Theta);
title('F-SAPES�㷨')
ylabel('��һ��������/dB')
xlabel('�ռ�Ƕ�')