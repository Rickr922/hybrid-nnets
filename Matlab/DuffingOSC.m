%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%           Simple and Duffing Oscillators
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear
close all

%%%%% Custom Parameters
SR = 1000;         %sample rate
k = 1 / SR;
durSec = 120.1;         %time of simulation (sec)

timeSamples = floor(durSec/k);
timeVec = (0:timeSamples-1)*k;
timeVecSec = (0:timeSamples-1)*SR/timeSamples;

% omega0 = 1;%sqrt(100); %radian freq (linear)
% % gamma = -0:100/timeSamples:500;        %gamma
% gamma = -49.506*ones(timeSamples,1);

omega0 = 1;
gamma = 1*ones(timeSamples,1);

steps = 5;
if steps == 1
    displacements = -1;
    velocities = -1;
elseif steps > 1
    displacements = -1:2/steps:1;          %initial displacement
    velocities = -1:2/steps:1;           %initial velocity
else
    disp('negative steps!')
    return
end

%+++++++++++++++++++++++++++++++++++
%%%% check stability condition
if k > 2/omega0
disp('stability condition violated!')
return
end

%+++++++++++++++++++++++++++++++++++
outDO = zeros(steps^2,timeSamples);
outSHO = zeros(steps^2,timeSamples);

pDO = zeros(steps^2,timeSamples);
pSHO = zeros(steps^2,timeSamples);

%-- main loop
for i = 1:steps
    for j = 1:steps
        %%%% Initialization
        % numerical initial conditions linearly-implicit scheme
        xPrev = displacements(i);
        vel0 = velocities(j);
        x = xPrev + k*vel0 - 0.5*k^2*(-omega0^2*xPrev-gamma(1)*xPrev^3); % implement second-order accurate initial condition
        % numerical initial conditions SHO
        vPrev = displacements(i);
        v = vPrev + k*vel0 - 0.5*k^2*(-omega0^2*vPrev);
        for n = 1 : timeSamples
            
            %Duffing
            xNext = x*(2-omega0^2*k^2)/(1+(gamma(n)*k^2*x^2/2)) - xPrev;
            outDO(j + 5*(i-1),n) = x;
            pDO(j + 5*(i-1),n) = (xNext - x)/k;
            xPrev = x; 
            x = xNext;
            
            %Linear
            vNext = v*(2-omega0^2*k^2) - vPrev;
            outSHO(j + 5*(i-1),n) = v;
            pSHO(j + 5*(i-1),n) = (vNext - v)/k;
            vPrev = v;
            v = vNext;
        end
    end
end

figure(1)
for i = 1:25
    plot(outDO(i,:));
    hold on
end
hold off
figure(2)
for i = 1:25
    plot(outSHO(i,:));
    hold on
end
hold off
figure(3)
for i = 1:25
    plot(outDO(i,:),pDO(i,:));
    hold on
end
hold off
figure(4)
for i = 1:25
    plot(outSHO(i,:),pSHO(i,:));
    hold on
end
hold off