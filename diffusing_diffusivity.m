function p=diffusing_diffusivity(N,h)


t=h+[0:N]*h;

x=zeros(1,N+1);

x(1)=normrnd(0,sqrt(1/2));

for i=1:N+1
    

 x(i+1)=x(i)-x(i)*h+sqrt(h)*randn;
 V(i)=x(i)^2;
 
end

p=V(2:end);

end