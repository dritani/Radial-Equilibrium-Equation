function [] = PropPlot2D ()


% Number of Vertical Lines
NSTATN = 51;
% Number of Horizontal Lines
NSTRM = 11;
% Leading Edge Number
NLE = 21;
% Trailing Edge Number
NTE = 31;
% Radius of Hub
RHUB = 0.45;
% Radius of shroud
RSHROUD = 0.50;
% Define DELTAR by knowing the number of stations between hub and shroud
DELTAR = ( RSHROUD - RHUB ) / ( NSTRM - 1 );
DELTAZ = DELTAR;

% Dimension the necessary variables
RADIUS = zeros (NSTATN,NSTRM);
Z = zeros (NSTATN,NSTRM);

% Initialize the necessary variables
for i=1:NSTATN;
    for j=1:NSTRM;
        
        % Radius is straightforward; interpolation between RHUB and RSHROUD
        RADIUS(i,j) = RHUB + (j - 1) * ((RSHROUD - RHUB) / ( NSTRM - 1));
        Z(i,j) = DELTAZ * (i - 1);
        
    end
end


% Read the necessary data
CZ1=dlmread('cz_comp.txt');
CZ2=dlmread('cz_incomp.txt');
CZ3=dlmread('cz_analytical.txt');

CR1=dlmread('cr_comp.txt');
CR2=dlmread('cr_incomp.txt');
CR3=dlmread('cr_analytical.txt');

DENS1=dlmread('dens_comp.txt');
DENS2=dlmread('dens_incomp.txt');
DENS3=dlmread('dens_analytical.txt');

BETA1=dlmread('beta_comp.txt');
BETA2=dlmread('beta_incomp.txt');
BETA3=dlmread('beta_analytical.txt');


% Plot everything versus Radius
figure(1)
hold on
plot(RADIUS(NTE,1:NSTRM),CZ1(NTE,1:NSTRM),'-b');
plot(RADIUS(NTE,1:NSTRM),CZ2(NTE,1:NSTRM),'--g');
plot(RADIUS(NTE,1:NSTRM),CZ3(NTE,1:NSTRM),'or');
legend('Compressible','Incompressible','Analytical')
ylabel ('CZ')
xlabel ('Radius')
hold off

figure(2)
hold on
plot(RADIUS(NTE,1:NSTRM),CR1(NTE,1:NSTRM),'-b');
plot(RADIUS(NTE,1:NSTRM),CR2(NTE,1:NSTRM),'--g');
plot(RADIUS(NTE,1:NSTRM),CR3(NTE,1:NSTRM),'or');
legend('Compressible','Incompressible','Analytical')
ylabel ('CR')
xlabel ('Radius')
hold off

figure(3)
hold on
plot(RADIUS(NTE,1:NSTRM),DENS1(NTE,1:NSTRM),'-b');
plot(RADIUS(NTE,1:NSTRM),DENS2(NTE,1:NSTRM),'--g');
plot(RADIUS(NTE,1:NSTRM),DENS3(NTE,1:NSTRM),'or');
legend('Compressible','Incompressible','Analytical')
ylabel ('Density')
xlabel ('Radius')
axis([0.45 0.50 1 2.5])
hold off

figure(4)
hold on
plot(RADIUS(NTE,1:NSTRM),BETA1(NTE,1:NSTRM),'-b');
plot(RADIUS(NTE,1:NSTRM),BETA2(NTE,1:NSTRM),'--g');
plot(RADIUS(NTE,1:NSTRM),BETA3(NTE,1:NSTRM),'or');
legend('Compressible','Incompressible','Analytical')
ylabel ('Beta')
xlabel ('Radius')
hold off


end