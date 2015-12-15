function Global()

global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt t;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;

global TempEnv; % K

TempEnv=300; %K

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
KeH=5.6e-2; KaH=1e-2; %J/(cm s K); z direction
KeV=5.6; KaV=5; %J/(cm s K); r direction

g=3.5e14; %W/(cm3 K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rMin=0; rMax=20*nm; %cm
zMin=0; zMax=10*nm; %cm

Nr=50;Nz=30;

rNum=Nr+1; zNum=Nz+1;
dr=(rMax-rMin)/Nr;
dz=(zMax-zMin)/Nz;

tBegin=0; tEnd=100*1e-15; %s
Nt=1e4;
dt=(tEnd-tBegin)/Nt;

M=sparse(rNum*zNum*2);
T1=sparse(rNum*zNum*2,1);

T0=sparse(rNum*zNum*2,1);
C=sparse(rNum*zNum*2,1);

end
