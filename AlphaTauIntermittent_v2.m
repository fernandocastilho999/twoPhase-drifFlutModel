function [alpha, TW] = AlphaTauIntermittent_v2(JG, JL, J, DenG, DenL,...
        VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,P) 


DenC = DenG ; 
DenF = DenL ; 
VisC = VisG ; 
VisF = VisL ; 

% Frequency, Gregory & Scott (1969)
f = fun_freqSchulkes(DenL, VisL, D,  theta, g,JL, J); %Schulkes(2011)

RLS = 1 ;



% % % % 
% calculation of UT
Fr = J / ( dRho*g*D/DenL ) ;
Eo = dRho*g*D^2/surTen ;
if Fr<3.5 
    CoT = 1 + 0.2 * (sin(theta))^2  ;
    CinfT =  (0.542 - 1.76 / (Eo^0.56)) * cos(theta) + ...
        ( 0.345*sin(theta) ) / ( 1+3805*(Eo^-3.06) )^0.58 ;
else
    CoT = 1.2  ;
    CinfT = ( 0.345*sin(theta) ) / ( 1+3805*(Eo^-3.06) )^0.58 ;
end
VinfT = CinfT * ( dRho*g*D/DenL )^(1/2) ;
% % % % UT =  1.2 + 0.35 * sqrt(g*D) * sin(theta)
% % % % 
% % % % % calculation of UT
% % % % regimen = 1 ; % 1:agitado , 2:distorcido
% % % % if regimen==1
% % % %     n = 0;
% % % % else
% % % %     n=7/4;
% % % % end
% % % % CoB = 1 + 0.2 * (sin(theta))^2 ;
% % % % CinfB = 1.54 * (RS^n) * (Eo^(-1/4)) * sin(theta) ;
% % % % VinfB = CinfB * ( dRho*g*D/DenL )^(1/2) ;
% % % % UB = CoB * J + VinfB 
% % % % 
% % % % % Calculation of Alpha
Co = CoT ;
VGJ = VinfT ; 
alpha =  JG / ( Co*J+VGJ ) ;
% % % % 
Co = 1.12 ;
ene = 0 ;

contlinhas = 0 ;
linhas = 1;

[res1,res2,LF, LS, RGB, TWC, TWF] = horizontal_v2(Co, RLS, DenL, VisL, DenG, VisG, ...
            surTen, D, g, theta, JL, JG, J, f, ene);
 
        
        
% Liquid fraction in the slug
% Gregory et al (1978)
RS = ( 1 + (J/8.66)^1.39 ) ^(-1)  ;

alpha = (1-RS)*LS + (RGB)*LF ;

%%%
% [LF,LS,TWC,TWF] = ...     
% 	horizontal...
%     (JG, JL, J , DenC , DenF , VisC , VisF, ...
%     D, A , theta, rug, g, surTen, RS, dRho, UB, f) ;




DenS = (1-RS)*DenG + RS * DenL; 
VisS = (1-RS)*VisG + RS * VisL;

ReS = DenS * J * D / VisS  ;
DS = D; 
CfS = ( -3.6 * log( (rug / (3.7*DS)).^1.11 + 6.9/ReS )) .^(-2) ;
tauWS =  CfS * DenS * J * abs(J) / 2 ;
TWS =  tauWS * S / A ;

%beta <1; 
beta = (LF * PL/P) / ( LF *PL/P + LS ) ;
TW = beta * (TWC + TWF ) + (1-beta)*TWS ;

% TW= RLS;



%%%%%%%%
% TW=0;
