clear;

format long

% -- Laser Properties --
c = 3.*10.^8; % Speed of light in m/s
wavelength = 800.*10.^-9;
k = 2.*pi./wavelength; % Wave vector
freq = 2*pi*(c./wavelength); % Frequency
T = 1./freq; % Period
sigma = 80*10^-15; % FWHM of pulse


% -- Time Settings --
numsteps = 90*10^4; % Total steps
tf = 800*10^-15; % Final time (in seconds)
h = tf./numsteps; % Step size
tvec = zeros(1, numsteps); % Initialize vector
Emax = 1; % Arbitrary
twave1 = 0.5*tf;%initial time of pulse 1 - makes graph look nice
twave2 = 2.*twave1;

for i=2:numsteps % Step through and set time vector
    tvec(i) = tvec(i-1)+h;
end

%-- Tau Settings --
numstepsTau = .4*10^3; %total number of steps
tfTau = tf; %final value of tau
hTau = tfTau./numstepsTau; %step size of tau
tauvec = zeros(1,numstepsTau); %initialize vector

Itot = zeros(1,numstepsTau);

for i=2:numstepsTau % Step through and set tau vector
    tauvec(i) = tauvec(i-1)+hTau;
end


% phase = zeros(1, numsteps); % Initialize phaseshift (position) vector
% for i=2:numsteps % Set the phase vector off of the times using function
%     phase(i) = time2pos(tvec(i));
% end

for i=2:numstepsTau
    Wave1 = E1(Emax,freq,tvec,twave1,sigma);
    Wave2 = E2(Emax,freq,tvec,tauvec(i),twave2,sigma);
    E = Wave1+Wave2;
    I = Intensity(E);
    Itot(i)= trapz(tvec,(I.^2)); %numerical integrate function
    %scatter(tvec(1,2:end), I(1,2:end))

    %E = Wave1 ;
    %scatter(phase(1,2:end), I(1,2:end))
    %scatter(tv(1,2:end), I(1,2:end))
    %scatter(tvec(1,2:end), E(1,2:end));

    % function x = time2pos(t) % Simulates the piezo micrometer, varying 1 micron per 1/2 second
    % %x = 0.0000000000001*sin(2*t);
    % x = 0.0000000000001*(t);
    % end
end
%scatter(tauvec,Itot,1);
plot((tauvec*c),Itot)

function I = Intensity(E)
I = (8*(E).^2);
end

function E = E2(Emax,freq,tvec,tau,twave2,sigma)
numsteps = length(tvec);
E = zeros(1,numsteps);
for i=1:numsteps
    %E(i) = Emax*cos(k*(0)-freq*tvec(i)+ phase(i));
    E(i) = Emax*cos(freq*(tvec(i)-tau-twave2))*exp((-(tvec(i)+tau-twave2).^2)./(2*sigma^2));
end
end

function E = E1(Emax,freq,tvec,twave1,sigma)
numsteps = length(tvec);
E = zeros(1,numsteps);
for i=1:numsteps
    %E(i) = Emax*cos(k*(0)-freq*tvec(i));
    E(i) = Emax*cos(freq*(tvec(i)-twave1))*exp(-((tvec(i)-twave1).^2)./(2*sigma^2));
end
end