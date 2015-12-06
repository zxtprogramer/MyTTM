function Fun3()
tTot=100;
N=100;
dt=tTot/100000;
dr=1/N;

s=dt/dr/dr;

t=0;
T0=zeros(N+1,1);

for j=1:N+1
    T0(j,1)=1200;
end

T1=T0;

while t<tTot
    A=zeros(N+1); 
    C=zeros(N+1,1);
   
   
    T0=T1;
    for j=2:N
        A(j,j+1)=0.5*s+0.25/(j-1);
        A(j,j)=-1-s;
        A(j,j-1)=0.5*s-0.25/(j-1);
        
        C(j,1)=(-0.5*s-0.25/(j-1))*T0(j+1,1) + (s-1)*T0(j,1) + (-0.5*s+0.25/(j-1))*T0(j-1,1);
            
    end
    A(1,2)=s; A(1,1)=-s-1;
    C(1,1)=-s*T0(2,1) + (s-1)*T0(1,1);
    A(N+1,N+1)=1;
    C(N+1,1)=300;
        

    T1=inv(A)*C;
   
    plot(T1);
    pause
    t=t+dt;
end

end






