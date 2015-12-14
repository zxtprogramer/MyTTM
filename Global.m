function Global()

global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;




nm=1e-7; %cm
ps=1e-12; %s
eV=1.602176e-19; % J

Na=1.12875e23;  % /cm3
Ea=3*eV;
v0=4e14;  %Hz

%%%%%SiO2%%%%%%%%%%%%%%%%%
Ce=1; Ca=0.42; %J/(cm3K)
KeV=2; KaV=0.01; %J/(cm s K)
KeH=2; KaH=0.01; %J/(cm s K)

g=1.25e13; %W/(cm3 K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%HOPG%%%%%%%%%%%%%%%%%
Ce=3.73e-2; Ca=1.6; %J/(cm3K)
KeV=11; KaV=5e3; %J/(cm s K)
KeH=11; KaH=5e3; %J/(cm s K)

g=3e13; %W/(cm3 K)
g=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rMin=0; rMax=10*nm; %cm
zMin=0; zMax=10*nm; %cm

Nr=20;Nz=20;

rNum=Nr+1; zNum=Nz+1;
dr=(rMax-rMin)/Nr;
dz=(zMax-zMin)/Nz;

tBegin=0; tEnd=10*1e-15; %s
Nt=1e3;
dt=(tEnd-tBegin)/Nt;

M=zeros(rNum*zNum*2);
T1=zeros(rNum*zNum*2,1);

T0=zeros(rNum*zNum*2,1);
C=zeros(rNum*zNum*2,1);

end
