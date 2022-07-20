% PROGRAMA 5.1
clear all, close all, clc


% Pipe properties
L = 10.; %length 
D = 0.0254; %diameter
angle = 0.;  % angle in degrees
theta = angle * pi / 180.; % angle in radians
rug = 0.001; %rugosity

% Other data
To = 300.; % temperature
g = 9.81; % acceleration of gravity 
Ro = 287.; % constant of the gas (air)

% Iteration Calculations
dz0 = D; % approximated delta z, 40 times diameter
N = floor(L / dz0 + 1); % number of intervales (integer)
dz = L / ( N-1) ; % e   xact delta z
tol = 1e-5; % error tolerance of calculations

% First Geometrical calculations
S = pi * D; % perimeter of the pipe
A = pi / 4 * (D^2); % area of the pipe
eD = rug/D ;



% Flow properties 
% % % JL = .5; % superficial velocity of liquid
DenL = 1000.;  % density of liquid
VisG = 1e-5; % viscosity of gas
VisL = 1e-3;  % viscosity of liquid
surTen = 0.7; % surface tension of liquid in contact with gas


JJ = 4 ; 

% Perda de carga do escoamento de Água escoando sozinha
ReJ = DenL * JJ * D / VisL ;
fL = ( -1.8 * log(  (eD/3.7)^1.11 + 6.9/ReJ ) )^(-2) ; % Haaland Eq.
DPL = fL * L/D * DenL*JJ^2/2 ;
Pref = DPL ;



ss = 1.5:.25:3.5% JL
for t = 1:length(ss) ;
JL = JJ - ss(t) ;
JG(N) = ss(t) ;
t, jgas = JG(N)



% Outlet properties
% % % JG(N) = .5; % Superficial velocity of gas
P(N) = 1e5; % Pressure at outlet

% Second Calculations  
% % % J(N)  = JL + JG(N); % Superficial velocity of the mixture 
J(N) = JJ ;
DenG(N) = P(N) / (Ro * To); % Density of the ideal gas 
dRho = DenL-DenG(N); % difference of density 

GG = DenG(N) * JG (N); % Mass flux of gas 
GL = DenL * JL; % Mass flux of liquid

global JL DenL VisG VisL surTen D S A L theta rug g dRho

% Routine flowPattern
% Algorithm 4.2
 pattern=ModelGBar87(JG(N), DenG(N));
% fPat = ModelGBar87(jG, jL, DenG, DenL, VisL, VisG, sTen, D, AngR, g)    

z = L;
% PL = P(L);
Pwater(N) = P(N);
Pair(N) = P(N);
% % 
% % % Routines for VoidFraction and FrictionForce for each flow pattern
% % if pattern == 1 % Dispersed Flow
% %    %% algorithm 5.3 
% %    [alfa(N), TW(N)] = alfaTauDispersed(JG(N), JL, J(N), DenG(N), DenL, ...
% %        VisG, VisL, surTen, D, S, theta) ;
% % elseif pattern == 2 % Separated Flow
% %     %% algorithm 5.4
% %     [alfa(N), TW(N)] = alfaTauSeparated(JG(N), J(N), DenG(N), DenL, ...
% %         VisG, VisL, surTen, D, S, theta) ;
% % else % Intermittent Flow
% %     %% algorithm 5.5
    [alfa(N), TW(N)] = AlphaTauIntermittent_v2(JG(N), JL, J(N), DenG(N), DenL,...
        VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,P(N),P(N))  ;
% % end

% actual phase velocities (eq. 3.3 & 3.5)  - Line11 
UG(N) = JG(N) / alfa(N) ; % actual velocity of gas
UL(N) = JL / ( 1-alfa(N) ) ; % actual velocity of liquid

% Psi function - Line 12
Psi(N) = P(N) + GG * UG(N) + GL * UL(N) ;


%Density of the mixture & dPsidz - Line 13 
Den(N) = alfa(N) * DenG(N) + ( 1-alfa(N) ) * DenL ;

dPsidz(N) = -TW(N) - Den(N) * g * sin(theta) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations for i from N-1 to 1 
i= N ; 
xx(i) = (i-1)*dz ; % position at i=N


x = 1e5-1000:10:1e5+Pref;

% for k = 1:length(x) ;
    
    
    
for i = N-1:-1:1
    xx(i) = (i-1)*dz ;
    Psi(i) = Psi(i+1) - dPsidz(i+1)*dz ;
       
    [P(i),DenG(i),JG(i),J(i),UG(i),UL(i)] = ...
       bisectionP(x, Ro, To, GG, JL, DenL,...
        VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,P(N),Psi(i)) ;
    
    [del, TW(i)] = AlphaTauIntermittent_v2(JG(i), JL, J(i), DenG(i), DenL,...
        VisG, VisL, surTen, D, S, A,  theta,rug, g, dRho,P(N),P(i));
    clear del
    
    %Density of the mixture & dPsidz - Line 13 
    Den(i) = alfa(i) * DenG(i) + ( 1-alfa(i) ) * DenL ;

    dPsidz(i) = -TW(i) - Den(i) * g * sin(theta) ;
    
    
    

    % Perda de carga do escoamento de Água escoando sozinha
    ReJL = DenL * J(i) * D / VisL ;
    fL = ( -1.8 * log(  (eD/3.7)^1.11 + 6.9/ReJL ) )^(-2) ; % Haaland Eq.
    dPL = fL * dz/D * DenL*J(i)^2/2 ;
    Pwater(i) = Pwater(i+1) + dPL ;

    % Perda de carga do escoamento de Ar escoando sozinho
    ReJG = DenG(i) * J(i) * D / VisG;
    fG = ( -1.8 * log(  (eD/3.7)^1.11 + 6.9/ReJG ) )^(-2); % Haaland Eq.
    dPG= fG * dz/D * DenG(i)*J(i)^2/2 ;
    Pair(i) = Pair(i+1) + dPG ;
    
end
hold on
plot(xx,P-1e5,'x-')



dPm = ( P(1) - P(N) ) / L ;



tDP(t) = dPm;
tP(:,t)=P-1e5;

end  
plot(xx,Pwater(i)-1e5 ,'y')
grid on
xlabel('x (m)')
ylabel('P_{man} (Pa)')
title(['Pressão manométrica para J =',num2str(JJ),'m/s'])
legend('J_G = 1,50 m/s','J_G = 1,75 m/s','J_G = 2,00 m/s','J_G = 2,25 m/s','J_G = 2,50 m/s','J_G = 2,75 m/s','J_G = 3,00 m/s','J_G = 3,25 m/s','J_G = 3,50 m/s','J_G = 0 m/s (Só Líquido)')

patch( [ 0 0 L ], [Pair(1)-1e5 Pwater(1)-1e5 0 ] ,'y','EdgeColor','y','EdgeAlpha',0.4 ), alpha(.2)
set(gcf,'Position',[0 500 600 400])