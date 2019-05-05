

clear;

format long

% -- Laser Properties --
c = 3.*10.^8; % Speed of light in m/s
wavelength = 400.*10.^-9;
k = 2.*pi./wavelength; % Wave vector
freq = c./wavelength; % Frequency
omega = 2*pi*freq; % ANGULAR frequency
T = 1./freq; % Period
sigma = 54*10^-15; % FWHM of pulse 107*10^-25 for 800nm at 58


% -- Time Settings --
numsteps = 12*10^4; % Total steps
tf = (10)*sigma; %250*10^-15; % Final time (in seconds)
h = tf./numsteps; % Step size
tvec = zeros(1, numsteps); % Initialize vector
Emax = 1; % Arbitrary
twave1 = 0.5*tf;%initial time of pulse 1 - makes graph look nice
twave2 = 2.*twave1;

for i=2:numsteps % Step through and set time vector
    tvec(i) = tvec(i-1)+h;
end

%-- Tau Settings --
numstepsTau = 12*10^3; %total number of steps
tfTau = tf; %final value of tau
hTau = tfTau./numstepsTau; %step size of tau
tauvec = zeros(1,numstepsTau); %initialize vector

Itot = zeros(1,numstepsTau);

for i=2:numstepsTau % Step through and set tau vector
    tauvec(i) = tauvec(i-1)+hTau;
end


for i=2:numstepsTau
    Wave1 = E1(Emax,omega,tvec,twave1,sigma);
    Wave2 = E2(Emax,omega,tvec,tauvec(i),twave2,sigma);
    E = Wave1+Wave2;
    I = Intensity(E);
    Itot(i)= trapz(tvec,(I.^2)); %numerical integrate function
    %scatter(tvec(1,2:end), I(1,2:end))
end

plot((tauvec*c),Itot)

peaks = [];
fwhm = max(Itot)./2;
for i = 3:length(Itot)
    if Itot(i) > Itot(i-2) && Itot(i-1) > Itot(i)
        peaks = [peaks;Itot(i-1)];
    end
end
rangefwhm = peaks(end)-peaks(1)
numpeak = 0;
for i = 1:length(peaks)
    if peaks(i) > fwhm
        numpeak = numpeak+1;
    end
end
numpeak
function I = Intensity(E)
I = 8*(E).^2;
end

function E = E2(Emax,omega,tvec,tau,twave2,sigma)
numsteps = length(tvec);
E = zeros(1,numsteps);
for i=1:numsteps
    %E(i) = Emax*cos(k*(0)-freq*tvec(i)+ phase(i));
    E(i) = Emax*cos(omega*(tvec(i)-tau-twave2))*exp((-(tvec(i)+tau-twave2).^2)./(2*(sigma^2)));
end
end

function E = E1(Emax,freq,tvec,twave1,sigma)
numsteps = length(tvec);
E = zeros(1,numsteps);
for i=1:numsteps
    %E(i) = Emax*cos(k*(0)-freq*tvec(i));
    E(i) = Emax*cos(freq*(tvec(i)-twave1))*exp(-((tvec(i)-twave1).^2)./(2*(sigma^2)));
end
end




