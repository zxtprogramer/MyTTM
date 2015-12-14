function BoundaryCondition1(r,z,t)
%const temperatue boundary
  
global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt t;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;

global TempEnv;

    
    for i=1
        for j=2:zNum-1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num, num)=Ce + KeV/dr/dr*dt + KeH/dz/dz*dt;
            M(num, num+zNum)=-KeV/dr/dr*dt;
            M(num, num+1)=-KeH/2/dz/dz*dt;
            M(num, num-1)=-KeH/2/dz/dz*dt;

            M(numA, numA)=Ca + KaV/dr/dr*dt + KaH/dz/dz*dt;
            M(numA, numA+zNum)=-KaV/dr/dr*dt;
            M(numA, numA+1)=-KaH/2/dz/dz*dt;
            M(numA, numA-1)=-KaH/2/dz/dz*dt;

            if(T0(num,1)>T0(numA,1))
                flag=1;
            else
                flag=0;
            end

            C(num,1)=(Ce - KeV/dr/dr*dt - KeH*dt/dz/dz)*T0(num,1) + ...
                     (KeV*dt/dr/dr)*T0(num+zNum,1) + ...
                     (KeH*dt/2/dz/dz)*T0(num+1,1) + ...
                     (KeH*dt/2/dz/dz)*T0(num-1,1) + ...
                     g*dt*(T0(numA,1)-T0(num,1))*flag + ...
                     AFun(r,z,t)*dt;

            C(numA,1)=(Ca - KaV*dt/dr/dr - KaH*dt/dz/dz)*T0(numA,1) + ...
                      (KaV*dt/dr/dr)*T0(numA+zNum,1) + ...
                      (KaH*dt/2/dz/dz)*T0(numA+1,1) + ...
                      (KaH*dt/2/dz/dz)*T0(numA-1,1) + ...
                      g*dt*(T0(num,1)-T0(numA,1))*flag;
        end
    end

    for i=rNum
        for j=2:zNum-1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=TempEnv;
            C(numA)=TempEnv;
        end
    end

    for i=1:rNum
        for j=1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=TempEnv;
            C(numA)=TempEnv;
        end
    end

    for i=1:rNum
        for j=zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;

            M(num,num)=1;
            M(numA,numA)=1;

            C(num)=TempEnv;
            C(numA)=TempEnv;
        end
    end


end
