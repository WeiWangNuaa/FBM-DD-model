

clc
clear all

t0=100;
dt=0.1;
n0=t0/dt;
delta=[1:1:100];
V=[n0 delta+n0];

T=ceil(max(V)*dt);

N=T/dt;
x0=0;

H=0.8
a=round(N.^(3/2));

n=10000;

parfor i=1:n    

    X(i,:)=generate_sample_LE(T,dt,H,x0,V);
    Y(i,:)=generate_sample_MN(T,dt,H,x0,a,V);
    Z(i,:)=generate_sample_levy(T,dt,H,x0,V);

end

% filename=['VACF_H=' num2str(H) '_t0=' num2str(t0) '_T=' num2str(T) '_dt=' num2str(dt) '.mat']
% save(filename,'X', 'Y','Z','n','T','dt','t0')

for k=1:length(V)-1
   VACF_X(k)= mean((X(:,k+1)-X(:,k)).*(X(:,2)-X(:,1))/dt^2,1);
   VACF_Y(k)= mean((Y(:,k+1)-Y(:,k)).*(Y(:,2)-Y(:,1))/dt^2,1);
   VACF_Z(k)= mean((Z(:,k+1)-Z(:,k)).*(Z(:,2)-Z(:,1))/dt^2,1);
end



%%H=0.8
figure
d=(delta-1)*dt;
s=unique(round(logspace(-1,1,20)/dt));
loglog(d(s), VACF_X(s),'rs','markersize',8,'LineWidth',1.5)
hold on
loglog(d(s), VACF_Y(s),'bo','markersize',8,'LineWidth',1.5)
loglog(d(s), VACF_Z(s),'g^','markersize',8,'LineWidth',1.5)

d3=[1 10];
loglog(d3, 0.7*d3.^(2*H-2),'k--','markersize',8,'LineWidth',1.5)
loglog(d3, 0.2*d3.^(H-3/2),'k--','markersize',8,'LineWidth',1.5)
xlabel('$\Delta$','Interpreter','latex','Fontsize',13)
ylabel('$C^\delta(t,\Delta)$','Interpreter','latex','Fontsize',13) 
set(gca,'FontSize',16);
legend({'LE-FBM-DD','MN-FBM-DD','RL-FBM-DD'},'Interpreter','latex','Fontsize',14)    
legend('boxoff')
xlim([0.1 T])
ylim([1e-2 1e1])

%%%%inset
figure
semilogx(d(s), VACF_X(s)./d(s).^(2*H-2),'rs','markersize',8,'LineWidth',1.5)
hold on
loglog(d(s), VACF_Y(s)./d(s).^(2*H-2),'bo','markersize',8,'LineWidth',1.5)
loglog(d(s), VACF_Z(s)./d(s).^(2*H-2),'g^','markersize',8,'LineWidth',1.5)
xlabel('$\Delta$','Interpreter','latex','Fontsize',13)
ylabel('$\langle \overline{\delta^2(\Delta)}\rangle$','Interpreter','latex','Fontsize',13) 

yline(H*(2*H-1),'b-.','LineWidth',1.5)
yline(2/pi*H*(2*H-1),'r-.','LineWidth',1.5)
% yline(2*H*(2*H-1)*dt^(H-1/2)/(2*H+1),'g-.','LineWidth',1.5)
yline(H*(2*H-1)*gamma(H+1/2)^2/gamma(2*H)/sin(pi*H),'g-.','LineWidth',1.5)

legend({'LE-FBM-DD','MN-FBM-DD','RL-FBM-DD'},'Interpreter','latex','Fontsize',14)    
legend('boxoff')
set(gca,'FontSize',16);

