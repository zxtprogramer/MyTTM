function Global()

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global M C T1 T0;

nm=1e-7; %cm
ps=1e-12; %s

Ce=1; Ca=0.42; %J/(cm3K)
KeV=2; KaV=0.01; %J/(cm s K)
KeH=2; KaH=0.01; %J/(cm s K)

g=1.25e13; %W/(cm3 K)

%%%%test%%%%
%g=0;
%Ce=Ca; KeV=KaV; KeH=KaH;
%%%%%%%%%%%%


rMin=0; rMax=10*nm; %cm
zMin=0; zMax=15*nm; %cm

Nr=21;Nz=21;

rNum=Nr+1; zNum=Nz+1;
dr=(rMax-rMin)/Nr;
dz=(zMax-zMin)/Nz;

tBegin=0; tEnd=1*ps; %s
Nt=1e3;
dt=(tEnd-tBegin)/Nt;

M=zeros(rNum*zNum*2);
T1=zeros(rNum*zNum*2,1);

T0=zeros(rNum*zNum*2,1);
C=zeros(rNum*zNum*2,1);

end
