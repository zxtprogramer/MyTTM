function a=AFun(r,z,t)
%Heat Source

global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;


dE=1e3*eV;

r0=1*nm;
z0=0.6*nm;
t0=1e-15;
b=dE/(sqrt(2*pi^3)*r0^2*z0*t0);



Sdis=exp(-r/r0)*exp(-z/z0);

Tdis=exp(-(t)^2/2/t0^2);
a=b*Sdis*Tdis;



  
end
