%% Liquid Rocket Engine Chamber Geometry Model | Truman Abbe | 7/9/2026

%% Model Inputs:
% F      max thrust (N)
% P_c    chamber pressure (Pa)
% C_F    thrust coefficient ()
% D_c    chamber diameter (m)
% L_c    chamber length (m)
% Th     nozzle angle (deg)
% D_e    exit diameter (m)
% P_amb  ambient pressure (Pa)
% ga     specific heat ratio ()
%
%% Model Outputs:
% D_t    throat diameter (m)
% eps_c  contraction ratio ()
% L_n    nozzle length (m)
% L_ch   characteristic length (m)
% LDrat  chamber length to diameter ratio ()
% eps    expansion ratio ()
% M_e    exit Mach number ()
% P_e    exit pressure (Pa) 
% eps_ep exit to ambient pressure ratio ()
%
%% Design Considerations:
% The characterist length L_ch is recommended to be between 0.75 m and 1.25
% m. The L/D ratio LDrat is recommended to be between 4 and 6. The contraction 
% ratio eps_c is recommended to be between 2.5 to 10. Exit pressure ratio eps_ep 
% is recommended to be between 1 and 1.5.  
%
%% Candidate Geometries:
% Candidate 4:
% The contraction ratio of 13.6 is too high, and there could be rapid
% erroding of the nozzle convergent section.
% F     = 30.000000 N
% P_c   = 1000.000000 kPa
% C_F   = 1.300000 
% D_c   = 20.000000 mm
% L_c   = 100.000000 mm
% Th    = 15.000000 deg
% D_t   = 5.420558 mm
% eps_c = 13.613568 
% L_n   = 54.411218 mm
% L_ch  = 1.361357 m
%
% Candidate 5:
% The target max thrust is increased to 50 N to provide a contraction ratio
% of 8.16, within the recommended range. The chamber length is decreased to
% provide L/D ratio of 3.
% F = 50.000000 N
% P_c = 1000.000000 kPa
% C_F = 1.300000 
% D_c = 20.000000 mm
% L_c = 60.000000 mm
% Th = 15.000000 deg
% D_t = 6.997911 mm
% eps_c = 8.168141 
% L_n = 48.524459 mm
% L_ch = 0.490088 m
% LDrat = 3.000000
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Candidate 6:
% The L/D ratio is adjusted back to 5. The model is reworked to account for
% the exit pressure P_e and the exit to ambient pressure ratio eps_ep. An
% exit diameter D_e smaller than the chamber diameter D_c is specified. The 
% exit to ambient pressure ratio is 1.25, which is underexpanded. 
% C6 Inputs:
% F = 50.000000 N        | max thrust
% P_c = 1000.000000 kPa  | chamber pressure
% C_F = 1.300000         | thrust coefficient
% D_c = 20.000000 mm     | chamber diameter
% L_c = 100.000000 mm    | chamber length
% Th = 15.000000 deg     | nozzle angle
% D_e = 10.000000 mm     | exit diameter
% P_amb = 85.000000 kPa  | ambient pressure for SLC
% ga = 1.270000          | specific heats ratio
% C6 Outputs:
% D_t = 6.997911 mm      | throat diameter
% eps_c = 8.168141       | contraction ratio
% L_n = 29.864204 mm     | nozzle length
% L_ch = 0.816814 m      | characteristic length
% LDrat = 5.000000       | L/D ratio
% eps = 2.042035         | expansion ratio
% M_e = 2.124141         | exit Mach number
% P_e = 106.726746 kPa   | exit pressure 
% eps_ep = 1.255609      | exit to ambient pressure ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

F   = 50;      % (N)
P_c = 1000000; % (Pa)
C_F = 1.3;     % ()
D_c = 0.02;    % (m)
L_c = 0.1;     % (m)
Th  = 15;      % (deg)
D_e = 0.01;    % (m)
P_amb = 85000; % (Pa)
ga = 1.27;     % ()

H = chamber_geometry(F, P_c, C_F, D_c, L_c, Th, D_e, P_amb, ga);

fprintf('New Run: \n\n');

fprintf('Inputs:\n');
fprintf('F = %0.6f N\n', F);
fprintf('P_c = %0.6f kPa\n', P_c / 1000);
fprintf('C_F = %0.6f \n', C_F);
fprintf('D_c = %0.6f mm\n', D_c * 1000);
fprintf('L_c = %0.6f mm\n', L_c * 1000);
fprintf('Th = %0.6f deg\n', Th);
fprintf('D_e = %0.6f mm\n', D_e * 1000);
fprintf('P_amb = %0.6f kPa\n', P_amb / 1000);
fprintf('ga = %0.6f \n\n', ga);

fprintf('Outputs:\n');
fprintf('D_t = %0.6f mm\n',  H(1) * 1000);
fprintf('eps_c = %0.6f \n',  H(2));
fprintf('L_n = %0.6f mm\n',  H(3) * 1000);
fprintf('L_ch = %0.6f m\n',  H(4));
fprintf('LDrat = %0.6f \n',  H(5));
fprintf('eps = %0.6f \n',    H(6));
fprintf('M_e = %0.6f \n',    H(7));
fprintf('P_e = %0.6f kPa\n', H(8) / 1000);
fprintf('eps_ep = %0.6f \n', H(9));


function H = chamber_geometry(F, P_c, C_F, D_c, L_c, Th, D_e, P_amb, ga)

    D_t = sqrt( (4 * F) / (pi * P_c * C_F) );

    A_t = ( pi * D_t^2 ) / 4;
    
    eps_c = ( D_c^2 * pi ) / ( 4 * A_t );

    L_n = ( D_c - D_t ) / ( 2 * tan(deg2rad(Th)) ) + ( D_e - D_t ) / (2 * tan(deg2rad(Th)));

    L_ch = L_c * eps_c;

    LDrat = L_c / D_c;

    eps = ( D_e / D_t )^2;

    M_e = find_exit_Mach(eps, ga);

    P_e = P_c * (1 + ((ga-1)/2) * M_e^2)^(-1 * (ga/(ga-1)));

    eps_ep = P_e / P_amb;

    H = [D_t,   ... % H(1)
         eps_c, ... % H(2)
         L_n,   ... % H(3)
         L_ch,  ... % H(4)
         LDrat, ... % H(5)
         eps,   ... % H(6)
         M_e,   ... % H(7)
         P_e,   ... % H(8)
         eps_ep];   % H(9)

end

function M_e = find_exit_Mach(eps, ga)
    % esp = A_e / A_t
    % ga = specific heat ratio (e.g., 1.27)
    
    % Define the area-Mach equation as a function
    area_func = @(M) (1/M) * ((2/(ga+1)) * (1 + ((ga-1)/2) * M^2))^((ga+1)/(2*(ga-1)));
    
    % Solve for M where area_func(M) = esp
    % Initial guess: 2.0
    M_e = fzero(@(M) area_func(M) - eps, 2.0);
end