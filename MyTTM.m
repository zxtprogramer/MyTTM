function MyTTM()

Global();
StartCondition();

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global M C T1 T0;

t=tBegin;

while t<tEnd
    M=zeros(rNum*zNum*2);
    C=zeros(rNum*zNum*2,1);
    T0=T1;
    
%%%%%boundary condition%%%%%
    for i=1
        for j=2:zNum-1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;

            M(num, num)=Ce/dt + KeV/dr/dr + KeH/dz/dz;
            M(num, num+zNum)=-KeV/dr/dr;
            M(num, num+1)=-KeH/2/dz/dz;
            M(num, num-1)=-KeH/2/dz/dz;

            M(numA, numA)=Ca/dt + KaV/dr/dr + KaH/dz/dz;
            M(numA, numA+zNum)=-KaV/dr/dr;
            M(numA, numA+1)=-KaH/2/dz/dz;
            M(numA, numA-1)=-KaH/2/dz/dz;



            C(num,1)=(Ce/dt - KeV/dr/dr - KeH/dz/dz -g)*T0(num,1) + ...
                     (KeV/dr/dr)*T0(num+zNum,1) + ...
                     (KeH/2/dz/dz)*T0(num+1,1) + ...
                     (KeH/2/dz/dz)*T0(num-1,1) + ...
                     g*T0(numA,1);

            C(numA,1)=(Ca/dt - KaV/dr/dr - KaH/dz/dz - g)*T0(numA,1) + ...
                      (KaV/dr/dr)*T0(numA+zNum,1) + ...
                      (KaH/2/dz/dz)*T0(numA+1,1) + ...
                      (KaH/2/dz/dz)*T0(numA-1,1) + ...
                      g*T0(num,1);
        end
    end

    for i=rNum
        for j=2:zNum-1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=BoundaryCondition(r,z,t);
            C(numA)=BoundaryCondition(r,z,t);
        end
    end

    for i=1:rNum
        for j=1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=BoundaryCondition(r,z,t);
            C(numA)=BoundaryCondition(r,z,t);
        end
    end

    for i=1:rNum
        for j=zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=BoundaryCondition(r,z,t);
            C(numA)=BoundaryCondition(r,z,t);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=2:rNum-1
        for j=2:zNum-1

            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;
        
            M(num, num)=Ce/dt + KeV/dr/dr + KeH/dz/dz;
            M(num, num+zNum)=-KeV/2/dr/dr - KeV/4/(i-1)/dr/dr;
            M(num, num-zNum)=-KeV/2/dr/dr + KeV/4/(i-1)/dr/dr;
            M(num, num+1)=-KeH/2/dz/dz;
            M(num, num-1)=-KeH/2/dz/dz;

            M(numA, numA)=Ca/dt + KaV/dr/dr + KaH/dz/dz;
            M(numA, numA+zNum)=-KaV/2/dr/dr - KaV/4/(i-1)/dr/dr;
            M(numA, numA-zNum)=-KaV/2/dr/dr + KaV/4/(i-1)/dr/dr;
            M(numA, numA+1)=-KaH/2/dz/dz;
            M(numA, numA-1)=-KaH/2/dz/dz;


            C(num,1)=(Ce/dt - KeV/dr/dr - KeH/dz/dz -g)*T0(num,1) + ...
                   (KeV/2/dr/dr + KeV/4/(i-1)/dr/dr)*T0(num+zNum,1) + ...
                   (KeV/2/dr/dr - KeV/4/(i-1)/dr/dr)*T0(num-zNum,1) + ...
                   (KeH/2/dz/dz)*T0(num+1,1) + ...
                   (KeH/2/dz/dz)*T0(num-1,1) + ...
                   g*T0(numA,1) + ...
                   AFun(r,z,t);
            

            C(numA,1)=(Ca/dt - KaV/dr/dr - KaH/dz/dz)*T0(numA,1) + ...
                    (KaV/2/dr/dr + KaV/4/(i-1)/dr/dr)*T0(numA+zNum,1) + ...
                    (KaV/2/dr/dr - KaV/4/(i-1)/dr/dr)*T0(numA-zNum,1) + ...
                    (KaH/2/dz/dz)*T0(numA+1,1) + ...
                    (KaH/2/dz/dz)*T0(numA-1,1) + ...
                    g*T0(numA,1);

        end
    end
    
    T1=inv(M)*C;
    t=t+dt;


    ri=0:dr:rMax; zi=0:dz:zMax;
    [z,r]=meshgrid(zi,ri);
    Te=zeros(rNum,zNum); Ta=zeros(rNum,zNum);
    for i=1:rNum
        for j=1:zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            Te(i,j)=T0(num,1);
            Ta(i,j)=T0(numA,1);
        end
    end
    det(M)
    %pcolor(r,z,Te)
    figure(1)
    surf(r,z,Te);
    figure(2)
    surf(r,z,Ta);
    %surf(r,z,Ta);

    pause
end

end



