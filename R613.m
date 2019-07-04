clear;
close all;
clc;

%�����۲��ź�Un
M=4;                        %�˲�����ͷ��
N=1000;                     %������
f=[0.1 0.25 0.27];          %��һ��Ƶ��
SNR=[30 30 27];             %�����
sigma=1;                    %����
Am=sqrt(sigma*10.^(SNR/10));%�źŷ���
t=linspace(0,1,N);
phi=2*pi*rand(size(f));     %�����λ
vn=sqrt(sigma/2)*rand(size(t))+1i*sqrt(sigma/2)*rand(size(t));%���Ը�˹������
Un=vn;
for k=1:length(f)
    s=Am(k)*exp(1i*2*pi*N*f(k).*t+1i*phi(k));
    Un=Un+s;
end
Un=Un.';

%��������
A=zeros(M,N-M+1);%�����۲����
for n=1:N-M+1
    A(:,n)=Un(M+n-1:-1:n);
end
[U,S,V]=svd(A');
invphi=V*inv(S'*S)*V';%��������phi

%��������
p=1024;
f=linspace(-0.5,0.5,p);
omega=2*pi*f;
a=zeros(M,p);%ÿһ�о�Ϊһ������
for k=1:p
    for m=1:M
        a(m,k)=exp(-1i*omega(k)*(m-1));
    end
end

%����MVDR��
Pmvdr=zeros(size(omega));
for k=1:p
    Pmvdr(k)=1/(a(:,k)'*invphi*a(:,k));
end
%Pmvdr=abs(Pmvdr);
Pmvdr=abs(Pmvdr/max(abs(Pmvdr)));
Pmvdr=10*log10(Pmvdr);

h=1:p;
h=(h-512)/1024;
plot(h,Pmvdr);
xlabel('M=4��MVDR��');
ylabel('��һ��MVDR��');