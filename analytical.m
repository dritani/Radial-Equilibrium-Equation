function [] = analytical()

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
% Initial Stagnation density
DENSITY0 = 1.5;
% RPM
N=6000;
% Rotational speed
OMEGA = 2 * pi * N / 60;
% Define DELTAR by knowing the number of stations between hub and shroud
DELTAR = (RSHROUD-RHUB)/(NSTRM-1);

% Dimension necessary variables
RADIUS = zeros (NSTATN,NSTRM);
BETA = zeros (NSTATN,NSTRM);
RCU = zeros (NSTATN,NSTRM);
CZ = zeros (NSTATN,NSTRM);
CR = zeros (NSTATN,NSTRM);
DENSITY = zeros (NSTATN,NSTRM);

for i=1:NSTATN;
    for j=1:NSTRM;
        
        RADIUS(i,j) = RHUB + (j - 1) * DELTAR;
        
        % CZ; constant for Quasi-1D, incompressible, no intra-blade stations
        CZ(i,j) = 135.812;
        
        % CR = 0 for Quasi-1D, incompressible, no intrabladestations
        CR(i,j) = 0;
        
        % Density is constant
        DENSITY(i,j) = DENSITY0; 
        
        % Swirl distribution is given
        if (i<=NLE)
           RCU(i,j)= 39.3 ;
        elseif (i>NTE)
           RCU(i,j) = 117.8;
        else
           RCU(i,j) = 39.3 + ((i - NLE) / (NTE - NLE)) * (117.8 - 39.3);
        end
        
        BETA(i,j)= - rad2deg ( atan (( OMEGA * RADIUS(i,j) - RCU(i,j) / RADIUS(i,j) ) / CZ(i,j) ));
        

    end 
end

% Print the results to a text file
dlmwrite('cz_analytical.txt',CZ);
dlmwrite('cr_analytical.txt',CR);
dlmwrite('dens_analytical.txt',DENSITY);
dlmwrite('beta_analytical.txt',BETA);


end