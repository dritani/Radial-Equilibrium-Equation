function [] = BladePlot3D ()


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


% Opens the text files containing the BETA values for each run
% BETA1 = Compressible with losses
% BETA2 = Incompressible without losses
% BETA3 = Analytical solution
BETA1=dlmread('beta_comp.txt');
BETA2=dlmread('beta_incomp.txt');
BETA3=dlmread('beta_analytical.txt');

% Define color matrices for the 3 cases
C1=zeros(NSTATN,NSTRM,3);
C2=zeros(NSTATN,NSTRM,3);
C3=zeros(NSTATN,NSTRM,3);

% Initialize color matrices. Case 1 is red, Case 2 is green, Case 3 is blue
% The shading is done based on proximity to maximum or minimum value.
% Shading is only used to make blade shape more apparent.
for i=1:NSTATN;
    for j=1:NSTRM;
        C1(i,j,1)=abs(BETA1(i,j)/abs(max(BETA1(:))))^2;
        C2(i,j,2)=abs(BETA2(i,j)/abs(max(BETA2(:))))^2;
        C3(i,j,3)=abs(BETA3(i,j)/abs(min(BETA3(:))))^2;
    end
end


% Plots the 3D Blade shape using R vs Z vs Beta
hold on
figure(1)
surf(Z(NLE:NTE,1:NSTRM),RADIUS(NLE:NTE,1:NSTRM),BETA1(NLE:NTE,1:NSTRM),C1(NLE:NTE,1:NSTRM,:))
surf(Z(NLE:NTE,1:NSTRM),RADIUS(NLE:NTE,1:NSTRM),BETA2(NLE:NTE,1:NSTRM),C2(NLE:NTE,1:NSTRM,:))
surf(Z(NLE:NTE,1:NSTRM),RADIUS(NLE:NTE,1:NSTRM),BETA3(NLE:NTE,1:NSTRM),C3(NLE:NTE,1:NSTRM,:))
xlabel('Z')
ylabel('R')
zlabel('Beta')
legend('Compressible','Incompressible','Analytical')
grid on

end