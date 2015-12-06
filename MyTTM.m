function MyTTM()

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global EM0 AM0 M0 T0 EM1 AM1 M1 T1;
global Temp0;

t=tBegin;

while t<tEnd
    M1=zeros(rNum*zNum*2);
    M0=zeros(rNum*zNum*2);

    for i=1:rNum
        for j=1:zNum
            num=(i-1)*zNum + (j-1);
            
            
        end
    end
    t=t+dt;
end





end

