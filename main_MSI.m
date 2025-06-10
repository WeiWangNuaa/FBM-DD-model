
clc
clear all


T=120;
dt=0.1;
N=T/dt;
x0=0;

H=0.2
a=round(N.^(1));
V=[1:N];
n=10000;

parfor i=1:n  

    X(i,:)=generate_sample_LE(T,dt,H,x0,V);
    Y(i,:)=generate_sample_MN(T,dt,H,x0,a,V);
    Z(i,:)=generate_sample_levy(T,dt,H,x0,V);

end



% filename=['MSI_H=' num2str(H) '_T=' num2str(T) '_dt=' num2str(dt) '.mat']
% save(filename,'X', 'Y','Z','n','T','dt','V')


t0=5;
n0=t0/dt;

for k=1:N-n0
   MSI_X(k)= mean((X(:,k+n0)-X(:,n0)).^2,1);
   MSI_Y(k)= mean((Y(:,k+n0)-Y(:,n0)).^2,1);
   MSI_Z(k)= mean((Z(:,k+n0)-Z(:,n0)).^2,1);
end


%H=0.8
% V=floor(logspace(-1,log10(length(MSI_X)*dt),12)/dt)
% d=[V]*dt;
% figure
% semilogx(d,MSI_X(V)./d.^(2*H),'rs','markersize',5,'LineWidth',1.5)
% hold on
% loglog(d,MSI_Y(V)./d.^(2*H),'bo','markersize',5,'LineWidth',1.5)
% loglog(d,MSI_Z(V)./d.^(2*H),'g^','markersize',5,'LineWidth',1.5)
% 
% xlabel('$t$','Interpreter','latex','Fontsize',11)
% ylabel('$\left< x_\Delta^2(t)\right>/\Delta^{2H}$','Interpreter','latex','Fontsize',11) 
% set(gca,'FontSize',10);
% yline(1,'k-.','LineWidth',1.5)
% yline(2/pi,'k-.','LineWidth',1.5)
% yline(gamma(H+1/2)^2/gamma(2*H)/sin(pi*H),'k-.','LineWidth',1.5)
% ylim([0.4 1.8])
% % 
% 
% figure
% V1=floor(logspace(-1,log10(length(MSI_X)*dt),20)/dt)
% d1=[V1]*dt;
% loglog(d1,MSI_X(V1),'rs','markersize',8,'LineWidth',1.5)
% hold on
% loglog(d1,MSI_Y(V1),'bo','markersize',8,'LineWidth',1.5)
% loglog(d1,MSI_Z(V1),'g^','markersize',8,'LineWidth',1.5)
% loglog(d1,d1.^(2*H),'k--','LineWidth',1.5)
% xlabel('$t$','Interpreter','latex','Fontsize',16)
% ylabel('$\left< x_\Delta^2(t)\right>$','Interpreter','latex','Fontsize',16) 
% legend({'LE-FBM-DD','MN-FBM-DD','RL-FBM-DD'},'Interpreter','latex','Fontsize',16)    
% legend('boxoff')
% set(gca,'FontSize',16);
% xlim([0.1 1e2])

%%%H=0.2
figure
V1=floor(logspace(-1,log10(length(MSI_X)*dt),20)/dt)
d1=[V1]*dt;
loglog(d1,MSI_X(V1),'rs','markersize',8,'LineWidth',1.5)
hold on
loglog(d1,MSI_Y(V1),'bo','markersize',8,'LineWidth',1.5)
loglog(d1,MSI_Z(V1),'g^','markersize',8,'LineWidth',1.5)
loglog(d1,d1.^(2*H),'k--','LineWidth',1.5)
xlabel('$\Delta$','Interpreter','latex','Fontsize',16)
ylabel('$\left< x_\Delta^2(t)\right>$','Interpreter','latex','Fontsize',16) 
legend({'LE-FBM-DD','MN-FBM-DD','RL-FBM-DD'},'Interpreter','latex','Fontsize',16)    
legend('boxoff')
set(gca,'FontSize',16);
xlim([0.1 1e2])
%%numerical
h=0.01;
for j=1:length(V)
f=@(s) 2*(V(j)*dt-s)*(h)^(-2).*(abs(s-h).^(2*H)+abs(s+h).^(2*H)-2*abs(s).^(2*H))...
    .*(sqrt(1-exp(-2*s))+exp(-s).*atan(exp(-s)./sqrt(1-exp(-2*s))))/pi;
w(j)=integral(f,0,V(j)*dt);
end
plot(V*dt,w,'r-')


V=floor(logspace(-1,log10(length(MSI_X)*dt),12)/dt)
d=[V]*dt;
figure
semilogx(d,MSI_X(V)./d.^(2*H),'rs','markersize',5,'LineWidth',1.5)
hold on
semilogx(d,MSI_Y(V)./d.^(2*H),'bo','markersize',5,'LineWidth',1.5)
hold on
loglog(d,MSI_Z(V)./d.^(2*H),'g^','markersize',5,'LineWidth',1.5)

xlabel('$\Delta$','Interpreter','latex','Fontsize',11)
ylabel('$\left< x_\Delta^2(t)\right>/\Delta^{2H}$','Interpreter','latex','Fontsize',11) 
set(gca,'FontSize',10);
yline(1,'k-.','LineWidth',1.5)
yline(2/pi,'k-.','LineWidth',1.5)
yline(gamma(H+1/2)^2/gamma(2*H)/sin(pi*H),'k-.','LineWidth',1.5)
ylim([0.4 1.8])
% 
% 
V=[1:1:10 10:10:100 100:10:10000];
for j=1:length(V)
f=@(s) ((1+s).^(H-1/2)-s.^(H-1/2)).^2;
v(j)=(integral(f,0,t0/V(j)/dt)+1/2/H)*2*H;
end
plot(V*dt,v,'g-')
xlim([0.1 1e3])

% 
% h=0.001;
% for j=1:length(V)
% f=@(s) 2*(V(j)*dt-s)*(h)^(-2).*(abs(s-h).^(2*H)+abs(s+h).^(2*H)-2*abs(s).^(2*H))...
%     .*(sqrt(1-exp(-2*s))+exp(-s).*atan(exp(-s)./sqrt(1-exp(-2*s))))/pi;
% w(j)=integral(f,0,V(j)*dt);
% end
% plot(V*dt,w./(V*dt).^(2*H),'r-')

H=0.2324
gamma(H+1/2)^2/gamma(2*H)/sin(pi*H)


