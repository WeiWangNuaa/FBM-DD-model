function p=generate_sample_LE(T,dt,H,x0,V)

N=T/dt;   %采样点数
dBH=fGN(N+1,H,dt);
r(1)=x0;
D=diffusing_diffusivity(N+1,dt);
 
for i=1:N+1
    
 r(i+1)=r(i)+sqrt(2*D(i))*dBH(i);

end

p=r(V+1);



end