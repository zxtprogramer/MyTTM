global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;

Global();

rM=[0:dr:rMax];
zM=[0:dz:zMax];
tM=[0:dt:tEnd];


sumAFun=0;
num=0;
for t=tM
  for z=zM
    for r=rM
      sumAFun=sumAFun + 2*pi*r*dr*dz*dt*AFun(r,z,t)/eV;
    end
  end
  sumAFun
end

sumAFun
