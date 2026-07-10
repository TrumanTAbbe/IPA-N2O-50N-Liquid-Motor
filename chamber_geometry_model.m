%% Liquid Rocket Engine Geometry Model | Truman Abbe | 7/9/2026

% Inputs:
% F     max thrust (N)
% P_c   chamber pressure (Pa)
% C_F   thrust coefficient ()
% D_c   chamber diameter (m)
% L_c   chamber length (m)
% Th    nozzle angle (deg)
%
% Outputs:
% D_t   throat diameter (m)
% eps_c contraction ratio ()
% L_n   nozzle length (m)
% L_ch  characteristic length (m)

F   = 30;      % (N)
P_c = 1000000; % (Pa)
C_F = 1.3;     % ()
D_c = 0.02;    % (m)
L_c = 0.1;     % (m)
Th  = 15;      % (deg)

H = chamber_geometry(F, P_c, C_F, D_c, L_c, Th);

fprintf('new run \n\n\n\n\n');

fprintf('D_t = %0.6f mm\n', H(1) * 1000);
fprintf('eps_c = %0.6f \n', H(2));
fprintf('L_n = %0.6f mm\n', H(3) * 1000);
fprintf('L_ch = %0.6f m\n\n', H(4));

fprintf('F = %0.6f N\n', F);
fprintf('P_c = %0.6f kPa\n', P_c / 1000);
fprintf('C_F = %0.6f \n', C_F);
fprintf('D_c = %0.6f mm\n', D_c * 1000);
fprintf('L_c = %0.6f mm\n', L_c * 1000);
fprintf('Th = %0.6f deg\n', Th);

function H = chamber_geometry(F, P_c, C_F, D_c, L_c, Th)

    D_t = sqrt( (4 * F) / (pi * P_c * C_F) );

    A_t = ( pi * D_t^2 ) / 4;
    
    eps_c = ( D_c^2 * pi ) / ( 4 * A_t );

    L_n = ( D_c - D_t ) / ( tan(deg2rad(Th)) );

    L_ch = L_c * eps_c;

    H = [D_t, eps_c, L_n, L_ch];

end

