function [X]=generate_sample_levy(T,dt,H,x0,V)

N=T/dt;
v=randn(1,N);
x(1)=x0;
D=diffusing_diffusivity(N+1,dt);

for j=1:length(V)
   
      i=[1:1:V(j)];    
      r=H-1/2;   
      s=2*H;
      % w=((V(j)-i+1).^(r+1)-(V(j)-i).^(r+1)).*dt.^(r+1)./(r+1)./dt.*sqrt(4*H*D(i)); 
      w=(((V(j)-i+1).^s-(V(j)-i).^s)*dt^(s)/(s*dt)).^(1/2).*sqrt(4*H*D(i)); 
      mfbm(j)=sum(w(1:V(j)).*v(1:V(j))*sqrt(dt));
      
end

X=mfbm+x0;





end