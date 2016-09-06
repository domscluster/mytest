
%fileID = fopen('/home/milan/LAPLACE.TRASFORM.SCRIPT/temp/num2.txt','r');
%formatSpec = '%f %f';
%sizeC = [2 Inf];
%C = fscanf(fileID,formatSpec,sizeC)
%fclose(fileID);
%C=C'
%t=C(1:100000 , 1)
%y=C(1:100000 ,2)

A1=0.30623
tau1= 6.449449217
A2= 0.70
tau2 = 45.454545455
dt=0.01
tmax=1000;
t=0:dt:tmax-dt;
y=A1*exp(-t/tau1)+A2*exp(-t/tau2);
 
% Zero padding (to avoid periodicity of fft)
y=[y 0.*y];
 
% Make laplace{d/dt{...}} using fft, and cal. frequency and omega
f=(0:(length(y)-1))./(dt*length(y));
omega=2*pi.*f;
sy=fft(gradient(y));
% sy=i.*omega.*fft(y)*dt-y(1)
 
% remove last half 
y=y(1:end/2);
sy=sy(1:end/2);
f=f(1:end/2);
omega=omega(1:end/2);
 
% Analytical expression of Laplace transform of exponential decay
%expoDecay=i.*omega*A./(i.*omega+1/tau)-A;
 
% Plot imaginary part
figure(1)
hold on
	plot(omega,-imag(-sy),'.');
%	plot(omega,-imag(-expoDecay),'r');
hold off
 
% ... and real part
figure(2)
hold on
	plot(omega,real(-sy),'.');
%	plot(omega,real(-expoDecay),'r');
hold off

figure(3)
hold on
        plot(real(-sy),-imag(-sy),'.');
%       plot(omega,real(-expoDecay),'r');
hold off

ImagNumLapTrans.png
