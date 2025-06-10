function X=generate_sample_MN(T,dt,H,x0,a,V)


N=T/dt;
v1=randn(1,a);
v2=randn(1,N);
D1=diffusing_diffusivity(N+1,dt);
D2=diffusing_diffusivity(a,dt);
for j=1:length(V)
       
    r=H-1/2;
    s=2*H;
    i=[1:1:V(j)];        
    % w=((V(j)-i+1).^(r+1)-(V(j)-i).^(r+1))*dt^(r+1)/(r+1)/dt.*sqrt(2*D1(i));
    w=(((V(j)-i+1).^s-(V(j)-i).^s)*dt^(s)/(s*dt)).^(1/2).*sqrt(2*D1(i)); 
    
    k=[1:1:a];  
    % m=(((V(j)+k).^(r+1)-(V(j)+k-1).^(r+1))*dt^(r+1)/(r+1)/dt-((k).^(r+1)-(k-1).^(r+1))*dt^(r+1)/(r+1)/dt).*sqrt(2*D2);
    m=((((V(j)+k).^s-(V(j)+k-1).^s)*dt^(s)/(s*dt)).^(1/2)-((k).^(r+1)-(k-1).^(r+1))*dt^(r+1)/(r+1)/dt).*sqrt(2*D2);

 
AH=sqrt(gamma(2*H+1)*sin(pi*H))/gamma(H+1/2);  
mfbm(j)=(sum(m.*v1)+sum(w(1:V(j)).*v2(1:V(j))))*sqrt(dt)*AH;
      
end

X=mfbm;

end


