
% function fPat = ModelGBar87(D, g, DenL, DenG, VisL, VisG, jL, jG, sTen, AngR)

function pattern = ModelGBar87(JG, DenG)
% function pattern=ModelGBar87(JG(N), JL, J(N) , DenG(N) , DenL, VisL, VisG, ...
%     sTen, D, L, S, A, Ang, rug, dRho, g);           

    global JL DenL VisG VisL surTen D S A L theta rug g dRho


    jG = JG;
    jL = JL;
    sTen=surTen;
    AngR = theta;

    fPat = 0;
    pattern =0 ;

    % Reynolds numbers
    ReL = DenL * jL * D / VisL;
    ReG = DenG * jG * D / VisG;
    % Determination of CL - CG and n-m exponents of PressureDrop equation
         if (ReL<4100)
            CL=16. ;
            n = 1. ;
         else

            CL = 0.046 ;
            n = 0.2 ;
         end

         if (ReG<4100)
             CG=16. ;
             m = 1. ;
         else
             CG = 0.046 ;
             m = 0.2 ;
         end

    % Determination of
    %%%%%%%%

        fM = CL * ( (jL+jG) * D / (VisL/DenL) )^-n ;


        dCD = 2. * ( 0.4*sTen / ( (DenL-DenG)* g ) )^0.5 ;

        dCB = 3./8. * (DenL / (DenL-DenG)) * (fM * (jL+jG)^2.) / (g * max(1e-3,cos(AngR))) ;

        dC = min(dCB,dCD) ;

        Crit1a = ( 0.725 + 4.15 * (jG/(jL+jG))^0.5 ) * ( sTen/DenL )^(3./5.) * ( 2*fM/D * (jL+jG)^3 )^(-2./5.) ; 

        %Print*,'Crit1',Crit1

   %%%%%%%%
    % Parameter of Lockhart-Martinelli
    X2 = ( (4.*CL/D) * ( DenL*jL*D / VisL )^-n * (DenL*(jL^2.)/2.) ) / ...
        ( (4.*CG/D) * ( DenG*jG*D / VisG )^-m * (DenG*(jG^2.)/2.) ) ;

    %print*, X2 
    X = sqrt(X2) ;

    Y = ( (DenL-DenG)*g*sin(AngR) ) / ( (4.*CG/D) * ( DenG*jG*D / VisG )^-m * (DenG*(jG^2)/2) ) ;


    % Obtain hL
    h_L=Calc_h_L2(X2, Y, m , n, jL, jG) ;
    %PRINT*, 'jL=',jL,'jG=',jG,'X=',X,'Y=',Y,'h_L=',h_L


    % Calculate F
    F = sqrt( (DenG)/(DenL-DenG) ) * (jG) / sqrt(D * g * cos(AngR)) ;

    %Criterion*

    h_x = 2. * h_L - 1. ;
    A_L = 0.25 * ( pi - acos(h_x) + h_x * sqrt(1 - h_x^2.) ) ;
    A_G = 0.25 * ( acos(h_x) - h_x * sqrt(1 - h_x^2.) ) ;
    A_P = A_L + A_G ;
    u_L = A_P / A_L ;
    u_G = A_P / A_G ;

    S_I = sqrt( 1 - h_x^2 ) ;

    S_L = pi - acos(h_x) ;
    D_L = 4 * A_L / S_L ;

    Crit2a = F^2. * ( 1./((1-h_L)^.2) * (u_G^2. * S_I)/(A_G)) ;


    Z = ( (4.*CL/D) * ( DenL*jL*D / VisL )^(-n) * (DenL*(jL^2)/2.) ) /...
            ( DenL * g * cos(AngR) ) ;
    Crit3a = 2.* (A_L/A_P)^2. * (1.-h_L) ;


    s = 0.01 ; 
    K = F * sqrt(ReL) ;
    Crit4a = 2. / ( (sqrt(u_L))*(u_G)*(sqrt(s)) ); 

    W = jL / sqrt(g*D) ;
    Crit4b = 1.5 * sqrt(h_L) * ( A_L / A_P ) ;


    % ANULAR OR NOT%
    alfaL = Calc_alfaL(X2,Y) ;
    %print*,jL,jG,'alfaL = ', alfaL
    Crit5a = ( (2. - 3./2.*alfaL) / ( (alfaL^3.) * (1.-3./2.*alfaL) ) ) * X2 ;

    Rsm = 0.48 ;
    Crit5b = alfaL/Rsm ;

    % THERE IS BUBBLE?
    Crit6a = 19. * ( ( ( (DenL-DenG)*sTen ) / ( (DenL^2)*g ) )^(1./2.) ) ;

    Clift = 0.8 ;
    gam = 1.3 ;
    Uo = 1.53 * ( ( ( g*(DenL-DenG)*sTen ) / (DenL^2.) )^(1./4.) ) ;
    Crit6b = (3./4.) * cos(pi/4.) * ( (Uo^2.)/g ) * ( Clift*(gam^2.)/D ) ; %#ok<NASGU>
    Crit6c = cos(AngR) / (sin(AngR))^2. ; %#ok<NASGU>

    % BUBBLE OR INTERMITTENT
    alfaC = 0.25 ;

    Crit7a = ((1.-alfaC)/alfaC)*jG - 1.53 * (1.-alfaC) * ...
        ( ( ( g*(DenL-DenG)*sTen ) / (DenL^2.) )^(1./4.) ) * sin(AngR) ;


    % SLUG OR ELONGATED%
    alfaS = 0.058 * ( ( dC * ( ( (2.*fM/D) * ((jL+jG)^3.) )^(2./5.) ) * ( (DenL/sTen)^(3./5.) ) - 0.725 )^2. ) ;
    Rs = 1. - alfaS ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE MODELS%%
    if ( (dC>=Crit1a) && (jL>=(jG*0.48/0.52)) ) %%%%%FECHADO1%%%%

        fPat = 20 ;
        pattern = 1 ;

    elseif ( (Crit2a<1.) && ( Z<Crit3a ) && ...
            (( (AngR<0.) && (W<Crit4b) ) | ...
            ( (AngR>=0.) && (K<Crit4a) )) ) 

        fPat = 21 ;
        pattern = 2 ;

    elseif ( (Crit2a<1.) && ( Z<Crit3a ) ) 

        fPat = 22 ;
        pattern = 2 ;

    elseif ( (Y<Crit5a) && (Crit5b<0.5) ) 

        fPat = 23 ;
        pattern = 2 ;

    %(AngR>=(80.*pi/180.))
    elseif ( ( (D<=Crit6a) ) && ...
            ( ( (jL>Crit7a) && (AngR>=(60.*pi/180.)) ) ) ) 

        fPat = 24 ;
        pattern = 1 ;

    elseif ( (Rs<1) && (AngR==0.) && ( Rs>0.48 ) ) 

        fPat = 27 ;
        pattern = 3 ;

    %elseif ( Rs<=0.48 ) 

     %   ModelGBar87 = 25

    elseif ( (Rs<0.48)&&(AngR>0.) && (AngR>=(70.*pi/180.))  )  

        fPat = 26 ;
        pattern = 3 ;

       
    else
        fPat = 25 ;
        pattern = 3 ;
    end
    
end


function h_L=Calc_h_L2(X2, Y, m , n, jL, jG)
%     IMPLICIT NONE
%     DOUBLE PRECISION :: jL, jG
%     DOUBLE PRECISION ::  f, f1, f2
%     DOUBLE PRECISION :: m, n
%     DOUBLE PRECISION :: h_L, h_L1, h_L2
%     DOUBLE PRECISION:: X2, Y
%  %   INTEGER :: nx
%     DOUBLE PRECISION :: fMSec2
%     DOUBLE PRECISION :: tol, var

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_L1 = 1e-5  ;
    h_L2 = 1. - 1.e-5 ;

    tol = 1e-7 ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        var = 1. ;
        while (var>=tol) 
            f1 = fMSec2(h_L1, X2, Y, m, n) ;
            f2 = fMSec2(h_L2, X2, Y, m, n) ; %#ok<NASGU>

            h_L = (h_L1 + h_L2) / 2. ;
            f = fMSec2(h_L, X2, Y, m, n) ;

            if (f*f1 < 0.)  ;
                h_L2 = h_L ;
            else
                h_L1 = h_L ;
            end

            var = (abs(h_L1-h_L2)) ;

        end

%     Calc_h_L2 = h_L ;
    

end


function f=fMSec2(h_L, X2, Y, m, n)
%     IMPLICIT NONE
%     DOUBLE PRECISION :: h_L, h_x
%     DOUBLE PRECISION :: m, n, X2, Y
%     DOUBLE PRECISION :: A_L, A_G, S_L, S_G, S_I, A_P, u_L, u_G, D_L, D_G
%     DOUBLE PRECISION, PARAMETER :: pi = 3.145159265


    h_x = 2.*h_L-1. ;

 
    A_L = 0.25 * ( pi - acos(h_x) + h_x * sqrt(1. - h_x^2.) )  ;
    A_G = 0.25 * ( acos(h_x) - h_x * sqrt(1. - h_x^2.) ) ;
    S_L = pi - acos(h_x) ;
    S_G = acos(h_x) ;
    S_I = sqrt( 1. - h_x^2. ) ;
    A_P = A_L + A_G ;
    u_L = A_P / A_L ;
    u_G = A_P / A_G ;
    D_L = 4. * A_L / S_L ;
    D_G = 4. * A_G / (S_G + S_I) ;

    f = (X2^1.) * ( ((u_L*D_L)^(-n)) * (u_L^2.) * (S_L / A_L) ) -...
        ( ((u_G*D_G)^(-m)) * (u_G^2.) * (S_G/A_G + S_I/A_L + S_I/A_G) ) + 4.*Y ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alfaL = Calc_alfaL(X2,Y)

% IMPLICIT NONE
% DOUBLE PRECISION :: fMBis2
% DOUBLE PRECISION :: X2, Y, alfaL1, alfaL2, alfaL, f1, f2, f
% DOUBLE PRECISION :: tol, var

alfaL1 = 1e-4 ;
alfaL2 = 1. - 1e-4 ;

tol = 1e-7 ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        var = 1. ;
        while (var>=tol)
            f1 = fMBis2(alfaL1,X2,Y) ;
            f2 = fMBis2(alfaL2,X2,Y) ;

            alfaL = (alfaL1 + alfaL2) / 2. ;
            f = fMBis2(alfaL,X2,Y) ;

            if (f*f1 < 0.)  
                alfaL2 = alfaL ;
            else
                alfaL1 = alfaL ;
            end

            var = (abs(alfaL1-alfaL2)) ;

        end

%     Calc_alfaL = alfaL ;
    

end


function f = fMBis2(alfaL,X2,Y)
% 
% IMPLICIT NONE
% DOUBLE PRECISION :: alfaL, X2, Y


    f = ( (1.+75.*alfaL) / ( ((1.-alfaL)^(5./2.))*alfaL ) ) - ( (1./(alfaL^3.))*X2 ) - Y ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
