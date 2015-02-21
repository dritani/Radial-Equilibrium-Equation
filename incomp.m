function [] = incomp ()


% Number of Vertical Lines
NSTATN = 51;
% Number of Horizontal Lines
NSTRM = 11;
% Leading Edge Number
NLE = 21;
% Trailing Edge Number
NTE = 31;
% Mass flow rate
XMASS = 30.4;
% Inlet Stagnation density
DENSITY0 = 1.5;
% Radius of Hub
RHUB = 0.45;
% Radius of shroud
RSHROUD = 0.50;
% RPM
N=6000;
% Rotational speed
OMEGA = 2 * pi * N / 60;
% Gas constants
RGAS = 287;
GAMAM = 3.5;
CP = GAMAM * RGAS;
% Inlet Stagnation conditions
T0 = 288;
P0 = DENSITY0 * T0 * RGAS;
% Define DELTAR by knowing the number of stations between hub and shroud
DELTAR = ( RSHROUD - RHUB ) / ( NSTRM - 1 );
DELTAZ = DELTAR;


% Dimension the necessary variables
PSI = zeros (NSTATN,NSTRM);
RHS = zeros (NSTATN,NSTRM);
CZ = zeros (NSTATN,NSTRM);
RCU = zeros (NSTATN,NSTRM);
PTOTAL = zeros (NSTATN,NSTRM);
RADIUS = zeros (NSTATN,NSTRM);
Z = zeros (NSTATN,NSTRM);
CR = zeros (NSTATN,NSTRM);
DENSITY = zeros (NSTATN,NSTRM);
ENTROPY = zeros (NSTATN,NSTRM);
HTOTAL = zeros (NSTATN,NSTRM);
A = zeros (NSTATN,NSTRM);
B = zeros (NSTATN,NSTRM);
BETA = zeros (NSTATN,NSTRM);


% Initialize the necessary variables
for i = 1:NSTATN;
    for j = 1:NSTRM;
        
        % Radius is straightforward; interpolation between RHUB and RSHROUD
        RADIUS(i,j) = RHUB + (j - 1) * DELTAR;
        Z(i,j) = DELTAZ * (i - 1);
        
        % We are given an intitial guess for the stream function
        PSI(i,j) = (RADIUS(i,j)^2 - RHUB^2) / (RSHROUD^2 - RHUB^2);
        DENSITY(i,j) = DENSITY0;
        PTOTAL(i,j) = P0;
        HTOTAL(i,j) = T0*CP;
        ENTROPY(i,j) = 0;
        RHS(i,j) = 0;
        
        % Before and after the rotor, swirl is constant (design)
        % Inside the rotor, interpolate
        if (i<=NLE)
           RCU(i,j)= 39.3 ;
        elseif (i>NTE)
           RCU(i,j) = 117.8;
        else
           RCU(i,j) = 39.3 + ((i - NLE) / (NTE - NLE)) * (117.8 - 39.3);
        end
        
    end
end


% The main loop. Maximum number of iterations is 100, but convergence
% is usually reached earlier than that
%
% Since this loop is so large, there will be no indentation for any
% statements inside it (merely a matter of coding legibility)
for p = 1:100;

% Step 2: Find Cz, Cr. At the boundaries, use forward or backward
% difference formula as explained by the TA
for i = 1:NSTATN;
    for j = 1:NSTRM;
        
        if (i==1)
            CR(i,j)=-(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAZ))...
                *(-3*PSI(i,j)+4*PSI(i+1,j)-PSI(i+2,j));
        elseif (i==NSTATN)
            CR(i,j)=-(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAZ))...
                *(3*PSI(i,j)-4*PSI(i-1,j)+PSI(i-2,j));
        else
            CR(i,j)=-(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAZ))...
                *(PSI(i+1,j)-PSI(i-1,j));            
        end

        
        
        if (j==1)
            CZ(i,j)=(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAR))...
                *(-3*PSI(i,j)+4*PSI(i,j+1)-PSI(i,j+2));
        elseif (j==NSTRM)
            CZ(i,j)=(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAR))...
                *(3*PSI(i,j)-4*PSI(i,j-1)+PSI(i,j-2));
        else
            CZ(i,j)=(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*(1/(2*DELTAR))...
                *(PSI(i,j+1)-PSI(i,j-1));            
        end
        
    end
end


% Define the norm of the density matrix. This will serve as a convergence
% criterion when compared with the updated density norm later on
OLDDENS = norm(DENSITY);

% Update density at inlet
% Ignore this for incompressible
%for j = 1:NSTRM;
%    CSQ = (RCU(1,j) / RADIUS(1,j))^2 + CZ(1,j)^2 + CR(1,j)^2;
%    HSTATIC = HTOTAL(1,j)-CSQ/2;
%    PSTATIC = PTOTAL(1,j) * (HSTATIC/HTOTAL(1,j))^GAMAM;
%    DENSITY(1,j) = PSTATIC /(RGAS*HSTATIC/CP);
%end


% Trace thermodynamic variables
% Code translated directly from FORTRAN77 code found in the notes
for i = 2:NSTATN;
    
    LEFT=i-1;
    
    if (i>NLE && i<=NTE)
        LEFT=NLE;
    end
    
    NSTART=1;
    PSIDN=PSI(LEFT,NSTART);
    PSIUP=PSI(LEFT,NSTART+1);
    
    for j = 1:NSTRM;
        
        DESIRED=PSI(i,j);
        
        while (1);
            if (DESIRED <= PSIUP && DESIRED >= PSIDN)
                break;
            end
        
            NSTART=NSTART+1;
        
            if (NSTART == NSTRM)
                fprintf('Cannot trace streamline\n');
                return;
            end
        
            PSIDN=PSI(LEFT,NSTART);
            PSIUP=PSI(LEFT,NSTART+1);
        end
        
        DELTA = (DESIRED - PSIDN) / (PSIUP - PSIDN);
        ROTATE = 0;
        
        if (i>NLE && i<=NTE) 
            ROTATE=OMEGA;
        end
        
        RCU1=DELTA*(RCU(LEFT,NSTART+1)-RCU(LEFT,NSTART))+RCU(LEFT,NSTART);
        HTOTAL1=DELTA*(HTOTAL(LEFT,NSTART+1)-HTOTAL(LEFT,NSTART))+HTOTAL(LEFT,NSTART);
        PTOTAL1=DELTA*(PTOTAL(LEFT,NSTART+1)-PTOTAL(LEFT,NSTART))+PTOTAL(LEFT,NSTART);
        RAD1=DELTA*(RADIUS(LEFT,NSTART+1)-RADIUS(LEFT,NSTART))+RADIUS(LEFT,NSTART);
        
        
        DENS1=DELTA*(DENSITY(LEFT,NSTART+1)-DENSITY(LEFT,NSTART))+DENSITY(LEFT,NSTART);
        CZ1=DELTA*(CZ(LEFT,NSTART+1)-CZ(LEFT,NSTART))+CZ(LEFT,NSTART);
        CR1=DELTA*(CR(LEFT,NSTART+1)-CR(LEFT,NSTART))+CR(LEFT,NSTART);
        ENTROP1=DELTA*(ENTROPY(LEFT,NSTART+1)-ENTROPY(LEFT,NSTART))+ENTROPY(LEFT,NSTART);
        
        C1SQ = CZ1^2+CR1^2+(RCU1/RAD1)^2;
        HSTAT1=HTOTAL1-C1SQ/2;
        PSTAT1=PTOTAL1*((HSTAT1/HTOTAL1)^GAMAM);
        
        ROTALP1 = HTOTAL1-ROTATE*RCU1;
        HOR2=ROTALP1+((ROTATE*RADIUS(i,j))^2)/2;
        HOR1=ROTALP1+((ROTATE*RAD1)^2)/2;
        POR1=PTOTAL1*(HOR1/HTOTAL1)^GAMAM;
        POR2IDL=POR1*((HOR2/HOR1)^GAMAM);
        
        % Accounting for losses
        % Ignore this for incompressible
        %if (i>NLE && i<=NTE)
        %   OMEGLOS=0.03*(i-NLE)/(NTE-NLE);
        %   PLOSS=OMEGLOS*(POR1-PSTAT1);
        %else
            OMEGLOS=0;
            PLOSS=0;
        %end
        POR2=POR2IDL-PLOSS;
        
        RCU(i,j)=RCU1;
        % Redefine the swirl matrix inside the rotor
        if (i>NLE && i<=NTE)
            RCU(i,j)=39.3+((i-NLE)/(NTE-NLE))*(117.8-39.3);
        end

        HTOTAL(i,j)=HTOTAL1+ROTATE*(RCU(i,j)-RCU1);
        PTOTAL(i,j)=RCU1+(HTOTAL(i,j)-HTOTAL1)/ROTATE;

        CU=RCU(i,j)/RADIUS(i,j);
        VU=CU-ROTATE*RADIUS(i,j);
        V2SQ=VU^2+CZ(i,j)^2+CR(i,j)^2;
        C2SQ=CU^2+CZ(i,j)^2+CR(i,j)^2;
        
        HSTATIC=HOR2-V2SQ/2;
        HTOTAL(i,j)=HSTATIC+C2SQ/2;
        PTOTAL(i,j)=POR2*(HTOTAL(i,j)/HOR2)^GAMAM;
        PSTATIC=PTOTAL(i,j)*(HSTATIC/HTOTAL(i,j));
        
        % Ignore this for incompressible
        % DENSITY(i,j)=PSTATIC/(RGAS*HSTATIC/CP);
        
        ENTROPY(i,j)=CP*log(HTOTAL(i,j)/HTOTAL1)-RGAS*log(PTOTAL(i,j)/PTOTAL1)+ENTROP1;
        
    end
end


% Recalculate the norm of the density matrix now that it has been updated
NEWDENS = norm(DENSITY);

% Calculate the norm of the RHS matrix. This is the 2nd convergence
% criterion
OLDRHS = norm(RHS);

% Calculate RHS of the REE equation by using the data computed above
for i = 1:NSTATN;        
    for j = 2:NSTRM-1;
        CU=RCU(i,j)/RADIUS(i,j);
        TSTATIC=(HTOTAL(i,j)-CZ(i,j)^2+CR(i,j)^2+CU^2)/(2*CP);
        RHS(i,j) = (( 2*pi ) / ( XMASS * CZ(i,j) * (-1) * 2 * DELTAR ))*...
            ( CU / RADIUS(i,j ) * ( RCU(i,j+1) - RCU(i,j-1) ) -...
            ( HTOTAL(i,j+1) - HTOTAL(i,j-1) ) + TSTATIC *...
            ( ENTROPY(i,j+1) - ENTROPY(i,j-1) ));
    end
end

% Recalculate the norm of the RHS
NEWRHS = norm(RHS);

% k is a variable that determines maximum number of iterations
k = 0;

% Next we have to update the stream function at every point using the RHS
% found above.
%
% Again, no indentation will be used inside this loop to improve
% code legibility
while(1);
    
k = k + 1;

% Max 20 iterations
if (k==20)
    break;
end

% All points except the inlet, since the conditions there are fixed
for i = 2:NSTATN;
    for j = 2:NSTRM-1;
        
        % The following variables represent the density and the radius
        % a half-step away from the current point. They are needed for
        % the PSI computation.
        % RJMINUS, for example, represents the RADIUS 1/2 units below
        % the current point, while DIPLUS represents the density
        % 1/2 units to the right of the current point
        if (i<NSTATN)
            RIPLUS = (RADIUS(i,j) + RADIUS(i+1,j)) / 2;
            DIPLUS = (DENSITY(i,j) + DENSITY(i+1,j)) / 2;
        end
        
        RIMINUS = (RADIUS(i,j) + RADIUS(i-1,j)) / 2;
        DIMINUS = (DENSITY(i,j) + DENSITY(i-1,j)) / 2;
        RJPLUS = (RADIUS(i,j) + RADIUS(i,j+1)) / 2;
        DJPLUS = (DENSITY(i,j) + DENSITY(i,j+1)) / 2;
        RJMINUS = (RADIUS(i,j) + RADIUS(i,j-1)) / 2;
        DJMINUS = (DENSITY(i,j) + DENSITY(i,j-1)) / 2;
        
        if (i<NSTATN)          
            A(i,j)=1/(1/(DIPLUS*RIPLUS)+1/(DIMINUS*RIMINUS)+1/(DJPLUS*RJPLUS)+...
                1/(DJMINUS*RJMINUS));
            B(i,j)=(PSI(i+1,j)/(DIPLUS*RIPLUS)+PSI(i-1,j)/(DIMINUS*RIMINUS)+...
                PSI(i,j+1)/(DJPLUS*RJPLUS)+PSI(i,j-1)/(DJMINUS*RJMINUS));
            
        % Here we account for the boundary conditions. The Neumann boundary
        % condition is needed for the exit of the duct.
        % The Dirichlett boundary conditions are automatically satisfied,
        % since the boundary itself does not change (PSI_hub=0,
        % PSI_shroud=1, always). Thus for a point right above the hub,
        % PSI(i,j-1) would be zero every time, meaning that the Dirichlett
        % condition is accounted for
        else
            A(i,j)=1/(2/(DIMINUS*RIMINUS)+1/(DJPLUS*RJPLUS)+1/(DJMINUS*RJMINUS));            
            B(i,j)=(2*PSI(i-1,j)/(DENSITY(i,j)*RADIUS(i,j))+...
                PSI(i,j+1)/(DJPLUS*RJPLUS)+PSI(i,j-1)/(DJMINUS*RJMINUS));        
        end
        
        % Update PSI
        PSIOLD=PSI(i,j);
        PSI(i,j)=A(i,j)*(B(i,j)+RHS(i,j)*(DELTAZ^2));
        
        % Optional point relaxation technique
        % DELTAPSI=PSI(i,j)-PSIOLD;
        % PSI(i,j)=PSI(i,j)+DELTAPSI*0.7;
        
    end
end


    % The convergence criterion for PSI. Breaks the infinite while(1)
    % loop it's contained in. Alternatively, if this does not converge then
    % the loop will terminate after 20 iterations from k defined right
    % below the while(1) loop
    if (abs(PSI(i,j)-PSIOLD)<10^(-5))
        break;
    end
    
% End of while(1) loop which updates PSI
end

% So PSI has been updated. Now it's time to see if it has fully
% converged. This occurs when RHS reaches a certain tolerance. The
% simultaneous convergence of density is imposed here. These values for
% tolerance are requested in page 14 of lecture 8
if (abs(NEWRHS-OLDRHS)<0.01 && abs(NEWDENS-OLDDENS)<0.001)
    break;
end
    
% End of big loop
end


% Calculate the angles
for i=1:NSTATN;
    for j=1:NSTRM;
        CU=RCU(i,j)/RADIUS(i,j);
        ROTATE=0;
        if (i>=NLE && i<=NTE)
           ROTATE=OMEGA;
        end
        VU=CU-ROTATE*RADIUS(i,j);
        % Compute the angle Beta and convert it to degrees
        BETA(i,j)= rad2deg(atan(VU/CZ(i,j)));
    end
end

% Print the results to a text file
dlmwrite('cz_incomp.txt',CZ);
dlmwrite('cr_incomp.txt',CR);
dlmwrite('dens_incomp.txt',DENSITY);
dlmwrite('beta_incomp.txt',BETA);


end