function MyTTM()

Global();
StartCondition();

global eV nm ps Ce Ca KeV KeH KaV KaH g;
global rMin rMax zMin zMax Nr Nz rNum zNum dr dz tBegin tEnd Nt dt t;

global M C T1 T0;  %T1 is the Temp at next time, T0 is now Temp
global nT Tmax;  %nT total number of rearrangement atoms; Tmax the max Temp during process;
  
global Ea Na v0;  %active energy; atom density; phonon frequency;

t=tBegin;
nT=0;
Tmax=zeros(rNum,zNum);
Tmelt=zeros(rNum,zNum);
sumAFun=0;

while t<tEnd
    M=sparse(rNum*zNum*2);
    C=sparse(rNum*zNum*2,1);
    T0=T1;
    
%%%%%boundary condition%%%%%
    BoundaryCondition2();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=2:rNum-1
        for j=2:zNum-1

            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;
        
            M(num, num)=Ce + KeV*dt/dr/dr + KeH*dt/dz/dz;
            M(num, num+zNum)=-KeV/2/dr/dr*dt - KeV/4/(i-1)/dr/dr*dt;
            M(num, num-zNum)=-KeV/2/dr/dr*dt + KeV/4/(i-1)/dr/dr*dt;
            M(num, num+1)=-KeH*dt/2/dz/dz;
            M(num, num-1)=-KeH*dt/2/dz/dz;

            M(numA, numA)=Ca + KaV*dt/dr/dr + KaH*dt/dz/dz;
            M(numA, numA+zNum)=-KaV*dt/2/dr/dr - KaV*dt/4/(i-1)/dr/dr;
            M(numA, numA-zNum)=-KaV*dt/2/dr/dr + KaV*dt/4/(i-1)/dr/dr;
            M(numA, numA+1)=-KaH*dt/2/dz/dz;
            M(numA, numA-1)=-KaH*dt/2/dz/dz;


            if(T0(num,1)>T0(numA,1))
                flag=1;
            else
                flag=0;
            end

            C(num,1)=(Ce - KeV*dt/dr/dr - KeH*dt/dz/dz)*T0(num,1) + ...
                   (KeV/2/dr/dr*dt + KeV/4/(i-1)/dr/dr*dt)*T0(num+zNum,1) + ...
                   (KeV/2/dr/dr*dt - KeV/4/(i-1)/dr/dr*dt)*T0(num-zNum,1) + ...
                   (KeH/2/dz/dz*dt)*T0(num+1,1) + ...
                   (KeH/2/dz/dz*dt)*T0(num-1,1) + ...
                   g*dt*(T0(numA,1)-T0(num,1))*flag + ...
                   AFun(r,z,t)*dt;
            

            C(numA,1)=(Ca - KaV*dt/dr/dr - KaH*dt/dz/dz)*T0(numA,1) + ...
                    (KaV*dt/2/dr/dr + KaV*dt/4/(i-1)/dr/dr)*T0(numA+zNum,1) + ...
                    (KaV*dt/2/dr/dr - KaV*dt/4/(i-1)/dr/dr)*T0(numA-zNum,1) + ...
                    (KaH*dt/2/dz/dz)*T0(numA+1,1) + ...
                    (KaH*dt/2/dz/dz)*T0(numA-1,1) + ...
                    g*dt*(T0(num,1)-T0(numA,1))*flag;

        end
    end
    
%    T1=inv(M)*C;
    T1=M\C;
    
    t=t+dt;

%%%%%%%%%%%%%%%%%deal with data%%%%%%%%%%%%%%%%%%%%%%%
    ri=0:dr:rMax; zi=0:dz:zMax;
    [zM,rM]=meshgrid(zi,ri);
    Te=zeros(rNum,zNum); Ta=zeros(rNum,zNum);

    for i=1:rNum
        for j=1:zNum
            num=(i-1)*zNum + (j-1) + 1;
            numA=num + rNum*zNum;
            r=(i-1)*dr; z=(j-1)*dz;
            
            Te(i,j)=T0(num,1);
	    Ta(i,j)=T0(numA,1);
	    
	    
	    if Ta(i,j)>Tmax(i,j)
	      Tmax(i,j)=Ta(i,j);
	    end
	    
	    if Tmax(i,j)>4000
		Tmelt(i,j)=1;
	    else
		Tmelt(i,j)=0;
            end
	      

	    nT=nT + (2*pi*r*dr*dz*dt)*Na*v0*exp(-(Ea*Na)/(Ta(i,j)*Ca));
	    sumAFun=sumAFun + 2*pi*r*dr*dz*dt*AFun(r,z,t);
	    
        end
    end

    hold off
    figure(1)
    surfc(rM,zM,Te,'FaceAlpha',0.5);
    %surfc(rM,zM,Te);
    sumAFun/eV
    nT
    %pcolor(r,z,Te)
    hold on;

    %surfc(rM,zM,Ta);
    surfc(rM,zM,Tmax,'FaceAlpha',1);

    hold off;
    figure(2);
    surfc(rM,zM,Tmelt);
    %pcolor(r,z,Ta)

end

end



