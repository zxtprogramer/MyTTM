function MyTTM()

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global EM0 AM0 M0 T0 EM1 AM1 M1 T1;
global Temp0;

t=tBegin;

while t<tEnd
    M1=zeros(rNum*zNum*2);
    M0=zeros(rNum*zNum*2);


    for i=2:rNum-1
        for j=2:zNum-1
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;


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

            T0=T1;
          
          
            
            
            
            
            
        end
    end
    t=t+dt;
end





end

