function MyTTM()

Global();

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global EM0 AM0 M0 T0 EM1 AM1 M1 T1;
global Temp0;

t=tBegin;

while t<tEnd
    M1=zeros(rNum*zNum*2);
    M0=zeros(rNum*zNum*2);
    A0=zeros(rNum*zNum*2,1);
    
%%%%%boundary condition%%%%%
    for i=1
        for j=2:zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;

            M1(num, num)=Ce/dt + KeV/dr/dr + KeH/dz/dz;
            M1(num, num+zNum)=-KeV/dr/dr;
            M1(num, num+1)=-KeH/2/dz/dz;
            M1(num, num-1)=-KeH/2/dz/dz;

            M0(num, num)=Ce/dt - KeV/dr/dr - KeH/dz/dz -g;
            M0(num, num+zNum)=KeV/dr/dr;
            M0(num, num+1)=KeH/2/dz/dz;
            M0(num, num-1)=KeH/2/dz/dz;
            M0(numA, num)=g;

            M1(numA, numA)=Ca/dt + KaV/dr/dr + KaH/dz/dz;
            M1(numA, numA+zNum)=-KaV/dr/dr;
            M1(numA, numA+1)=-KaH/2/dz/dz;
            M1(numA, numA-1)=-KaH/2/dz/dz;

            M0(numA, numA)=Ca/dt - KaV/dr/dr - KaH/dz/dz - g;
            M0(numA, numA+zNum)=KaV/dr/dr;
            M0(numA, numA+1)=KaH/2/dz/dz;
            M0(numA, numA-1)=KaH/2/dz/dz;
            M0(numA, num)=g;
        end
    end

    for i=rNum
        for j=1:zNum
            num=(i-1)*zNum + (j-1) + 1;
            M1(num,num)=1;
            M0(num,num)=1;
        end
    end

    for i=1:rNum
        for j=1
            num=(i-1)*zNum + (j-1) + 1;
            M1(num,num)=1;
            M0(num,num)=1;
        end
    end

    for i=1:rNum
        for j=zNum
            num=(i-1)*zNum + (j-1) + 1;
            M1(num,num)=1;
            M0(num,num)=1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:rNum
        for j=1:zNum

            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
     
            r=(i-1)*dr; z=(j-1)*dz;
            A0(num,1)=AFun(r,z);
            

            if i==1 || i==rNum || j==1 || j==zNum
                continue
            end
        
            M1(num, num)=Ce/dt + KeV/dr/dr + KeH/dz/dz;
            M1(num, num+zNum)=-KeV/2/dr/dr - KeV/4/(i-1)/dr/dr;
            M1(num, num-zNum)=-KeV/2/dr/dr + KeV/4/(i-1)/dr/dr;
            M1(num, num+1)=-KeH/2/dz/dz;
            M1(num, num-1)=-KeH/2/dz/dz;

            M0(num, num)=Ce/dt - KeV/dr/dr - KeH/dz/dz -g;
            M0(num, num+zNum)=KeV/2/dr/dr + KeV/4/(i-1)/dr/dr;
            M0(num, num-zNum)=KeV/2/dr/dr - KeV/4/(i-1)/dr/dr;
            M0(num, num+1)=KeH/2/dz/dz;
            M0(num, num-1)=KeH/2/dz/dz;
            M0(numA, num)=g;
            

            M1(numA, numA)=Ca/dt + KaV/dr/dr + KaH/dz/dz;
            M1(numA, numA+zNum)=-KaV/2/dr/dr - KaV/4/(i-1)/dr/dr;
            M1(numA, numA-zNum)=-KaV/2/dr/dr + KaV/4/(i-1)/dr/dr;
            M1(numA, numA+1)=-KaH/2/dz/dz;
            M1(numA, numA-1)=-KaH/2/dz/dz;

            M0(numA, numA)=Ca/dt - KaV/dr/dr - KaH/dz/dz - g;
            M0(numA, numA+zNum)=KaV/2/dr/dr + KaV/4/(i-1)/dr/dr;
            M0(numA, numA-zNum)=KaV/2/dr/dr - KaV/4/(i-1)/dr/dr;
            M0(numA, numA+1)=KaH/2/dz/dz;
            M0(numA, numA-1)=KaH/2/dz/dz;
            M0(numA, num)=g;

        end
    end
    T0=T1;
    
    T1=inv(M1)*(M0*T0 + A0);

    for i=rNum
        for j=1:zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            T1(num,1)=BoundaryCondition(i,j);
            T1(numA,1)=BoundaryCondition(i,j);
        end
    end

    for i=1:rNum
        for j=1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            T1(num,1)=BoundaryCondition(i,j);
            T1(numA,1)=BoundaryCondition(i,j);
        end
    end

    for i=1:rNum
        for j=zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            T1(num,1)=BoundaryCondition(i,j);
            T1(numA,1)=BoundaryCondition(i,j);
        end
    end
    t=t+dt;
    T1
    pause
end

end



