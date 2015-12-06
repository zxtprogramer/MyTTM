function Fun2()
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
    A=zeros(N+1); B=zeros(N+1);
    A(N+1,N+1)=1;
    B(N+1,N+1)=1;
    
    A(1,1)=-s-1; A(1,2)=s;
    B(1,1)=s-1; B(1,2)=-s;
    
    for j=2:N
        A(j,j+1)=0.5*s+0.25/j;
        A(j,j)=-1-s;
        A(j,j-1)=0.5*s-0.25/j;
        
        B(j,j+1)=-0.5*s-0.25/j;
        B(j,j)=s-1;
        B(j,j-1)=-0.5*s+0.25/j;
        
    end
        
    T0=T1;
    T1=inv(A)*B*T0;
    T1(N+1)=300;
    plot(T1);
    pause
    t=t+dt;
end

end






