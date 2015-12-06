function StartCondition()

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global M C T1 T0;

Temp0=300;
for i=1:rNum
    for j=1:zNum
        num=(i-1)*zNum + (j-1) + 1;
        numA=num + rNum*zNum;
        r=i*dr; z=j*dz;
        if z<2*zMax/3 && z>zMax/3
           T1(num,1)=Temp0+1000;
           T1(numA,1)=Temp0+1000;
           T0(num,1)=Temp0+1000;
           T0(numA,1)=Temp0+1000;
        else
           T1(num,1)=Temp0;
           T1(numA,1)=Temp0;
           T0(num,1)=Temp0;
           T0(numA,1)=Temp0;
        end

    end
end


end
