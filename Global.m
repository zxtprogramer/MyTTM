function Global()

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global EM0 AM0 M0 T0 EM1 AM1 M1 T1;
global Temp0;

nm=1e-7; %cm
ps=1e-12; %s

Temp0=300; %K

Ce=1; Ca=0.42; %J/(cm3K)
KeV=2; KaV=0.01; %J/(cm s K)
KeH=2; KaH=0.01; %J/(cm s K)
g=1.25e13; %W/(cm3 K)

rMin=0; rMax=10*nm; %cm
zMin=0; zMax=10*nm; %cm
Nr=10; Nz=10;
rNum=Nr+1; zNum=Nz+1;
dr=(rMax-rMin)/Nr;
dz=(zMax-zMin)/Nz;


tBegin=0; tEnd=10*ps; %s
Nt=1e5;
dt=(tEnd-tBegin)/Nt;


EM0=zeros(rNum,zNum);
AM0=zeros(rNum,zNum);
M0=zeros(rNum*zNum*2);
T0=zeros(rNum*zNum,1);

EM1=zeros(rNum,zNum);
AM1=zeros(rNum,zNum);
M1=zeros(rNum*zNum*2);
T1=zeros(rNum*zNum,1);

end
