function a=AFun(r,z,t)
%Heat Source

global nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt;
global M C T1 T0;
global eV;


rd=3*nm;
Td=100e-15;
dE=1e10*eV;

S=dE/4/pi/rd^2*exp(-(r^2)/4/rd^2)*exp(-z/rd);
Tt=(1/Td)*exp(-t/Td);
a=S*Tt;

  
end
