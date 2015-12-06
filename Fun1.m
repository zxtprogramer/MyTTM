function Fun1()
tTot=100;
N=10;
dt=tTot/10000;
dx=1/N;

s=dt/2/dx/dx;

t=0;
T0=zeros(N+1,1);

for j=1:N+1
    T0(j,1)=2000;
end
T1=T0
pause
while t<tTot
    A=zeros(N+1); B=zeros(N+1);
    A(1,1)=1;A(N+1,N+1)=1;
    B(1,1)=1;B(N+1,N+1)=1;
    
    for j=2:N
        A(j,j-1)=-s;
        A(j,j)=1+2*s;
        A(j,j+1)=-s;
        
        B(j,j-1)=s;
        B(j,j)=1-2*s;
        B(j,j+1)=s;
        
    end
        
    T0=T1;
    T1=inv(A)*B*T0;
    T1(1)=300;
    T1(N+1)=1200;
    plot(T1);
    pause
    t=t+dt;
end

end

