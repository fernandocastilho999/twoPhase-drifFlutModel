function [res1,res2,LB, LS, RGB, TWC, TWF] = horizontal_v2(Co, RLS, rhoL, miL, rhoG, miG, ...
    tensup, D, g, inclinacaor, jL, jG, J, f, ene)

%global Co RLS rhoL miL rhoG miG tensup D g inclinacaor jL jG J f ene

% DADOS DE ENTRADA %
%    inclinacaor = inclinacao*(PI)/180 

%Passo
dh = D/1e3 ;
primeira = 0 ;

%Velocidades Superficiais



alturaFilme = D*0.9  ;  %0.9 altura inicial
Vdrift = 1.54 *(tensup*g*(rhoL-rhoG)/rhoL^2  )^0.25 *RLS^ene*sin(inclinacaor) ;
UGS = J + Vdrift ;
ULS = J - UGS*(1-RLS) ;
%Co = CoM
VD = 0.54 *sqrt(g*D)*cos(inclinacaor)+0.35 *sqrt(g*D)*sin(inclinacaor) ; %Bendiksen 1984 
UT = Co*J + VD ;
%    CInf = CInfM
%    VD = CInf * sqrt(g*D)
%    UT = Co*J + VD
altB = alturaFilme ;
somaULB     = 0  ;
voltotal     = 0 ; 
invIncDB    = 0 ;
somaRLB        = 0 ;
compB       = 0 ;
loop        = 1   ;
incdBolha = 0 ;
k = 0 ;
ii=0;

condicao = 1 ; %Condição para entrada do primeiro loop
while (condicao == 1)
   

    if (RLS >= 1 ) 
        UGS = 0 ;
    else
        Vdrift = 1.54 *(tensup*g*(rhoL-rhoG)/rhoL^2  )^0.25 *RLS^ene*sin(inclinacaor);
        UGS = J + Vdrift ;
    end
    angInterno = 2 *acos(1 -2 *altB/D) ;
    RLB = (angInterno - sin(angInterno))/2 /pi ;

    %Cálculo interativo dos parâmetros associados ao desenho da bolha
    ULS =  J - UGS*(1-RLS) ;
    if (loop == 1) 
        ULBo = UT-((UT-ULS)*RLS/RLB) ;
        ALo = D^2 *(angInterno - sin(angInterno))/8  ;
    else
        ULBo = ULB;
        ALo = AL;
    end
    ULB = UT-((UT-ULS)*RLS/RLB) ; %Velocidade do líquido na bolha
    UGB =  UT - ((UT-UGS)*(1 -RLS)/(1 -RLB)) ; %Velocidade do gás na bolha

    %Liquido
    SL = angInterno*D/2 ; %Perimetro molhado do líquido
    AL = D^2 *(angInterno - sin(angInterno))/8  ; %Área do líquido
    DhL = 4 *AL/SL ; %Diametro hidráulico de liquido
    ReL = rhoL*ULB*DhL/miL ; %Reynolds do liquido
    fL = fun_fBlasius(ReL) ;
    %fl = fator_atrito_hall(Rel,dhl)       
    talL = rhoL*fL*abs(ULB)*ULB/2 ;

    %Gas
    SG = (2 *pi-angInterno)*D/2 ;  %Perimetro do Gas
    AG = D^2 *(2 *pi- angInterno + sin(angInterno))/8 ; %Area do gás
    SI = D*sin(angInterno/2 ) ; %Pertimetro interfacial
    DhG = 4 *AG/(SG+SI) ; %Diametro Hidráulico de gas
    ReG = rhoG*UGB*DhG/miG ; %Reynolds do gas
    fG = fun_fBlasius(ReG)  ;
    %fg = fator_atrito_hall(ReG,dhg)        
    talG = rhoG*fG*abs(UGB)*UGB/2 ;

    %Interfacial
    fI = 0.014 ; %Fator de atrito interfacial - Sugerido por: Cohen e Hanratty(1968); Shoham e Tailtel(1984) - Citado no artigo de Tailtel e Barnea(1990)
    talI = rhoG*fI*abs(UGB-ULB)*(UGB-ULB)/2 ;
    dR = 4 *sin(angInterno/2 )/(pi*D) ; %Derivada dRLB/dh
    termoFilmeL = talL*SL/AL ;
    termoFilmeG = talG*SG/AG ;
    termoFilmeI = talI*SI/((AL*AG/(AL+AG))) ;
    termoGravit = (rhoL - rhoG ) * g * sin(inclinacaor) ;
    inerciaL = rhoL * (UT-ULB)^2  * dR /  RLB   ;
    inerciaG = rhoG  * (UT-UGB)^2  * dR/ (1-RLB) ;
    forHidro = (rhoL - rhoG ) * g * cos(inclinacaor) ;
    numerador = termoFilmeL - termoFilmeG - termoFilmeI + termoGravit ;
    denominador = forHidro - inerciaG - inerciaL ;
    incdBolha   =   numerador / denominador ; %dh/dz
    %write(*,*),"Inclinacao",rls,rlsg,rlsb
    %pause


    %PARA CHEGAR NO NARIZ DA BOLHA   
    if(incdBolha > 0)  
        altB = altB - dh  ;
        somaULB     = 0  ;
        voltotal     = 0  ;
        invIncDB    = 0  ;
        somaRLB        = 0 ; 
        compB       = 0  ;
        loop        = 1   ;
        incdBolha = 0  ;
        %clear perfil

    %DH AJUSTÁVEL 
    else

        if(incdBolha < 0 ) 
            ii = ii +1 ;
            verificador =  abs(atand(numerador/denominador)) ;
            %write(11,*),compB,altB/D  %Plota o perfil
            res1(ii) = compB ;
            res2(ii) = altB/D  ;
            dCompB      =   invIncDB * dh  ; %Ponto anterior
            invIncDB    =   1 /incdBolha  ; %Novo ponto - Apenas usado no próximo loop

            %Integração Dos Parametros
            prodVelFilm =   ((ULB+ULBo)/2 )*((AL+ALo)/2 )*dCompB  ;
            somaULB     =   somaULB + prodVelFilm  ;
            volsecao    =   ((AL+ALo)/2 *dCompB) ;
            voltotal = voltotal + volsecao ;
            dCompRLB    =  dCompB * RLB ;
            somaRLB        =   somaRLB + abs(dCompRLB)                  ;
            LB = abs(compB) ;
            LS = UT/f-LB ;
            if (voltotal == 0 )   
                ULBm = 0  ;
                RLBm = 0  ;

            else
                ULBm = abs(somaULB/voltotal) ;
                RLBm = abs(somaRLB / LB) ;
            end

            %Cálculo do erro
            betam = 1 -(jL-ULBm*RLBm)/UT/RLS/(1 -RLBm/RLS) ;
            beta = LB/(LB+LS) ;
            ebeta = (betam-beta)*100 /beta  ; %ebeta = (beta-betam)*100/beta


            %Critério de Parada
            if (abs(ebeta) > 0.1 )  
                condicao = 1 ;
            else
                condicao = 0 ;
            end


            if (ebeta < 0 )    %Se erro < 0 não existe solução real
                if(k == 1)   %Verifica se precisa voltar 1 ponto
                    %Ponto novo = Ponto antigo
                    compB = compB - dCompB ;
                    altB = altB + dh ;
                    somaULB = somaULB - prodVelFilm  ;
                    voltotal = voltotal - volsecao ;
                    somaRLB    = somaRLB - abs(dCompRLB) ;
                    k = 0 ;
                end
                dh = dh/10 ; %Refina o passo
            else
                loop = loop + 1 ;
                compB = compB + invIncDB* dh ;
                altB  = altB - dh ;
                k = 1 ;
            end


         else
            if(k == 1)  
                compB = compB - invIncDB* dh ;
                altB = altB + dh ;
                somaULB = somaULB - prodVelFilm ;
                voltotal = voltotal - volsecao ;
                somaRLB    = somaRLB - abs(dCompRLB) ;
                loop = loop - 1 ;
                k = 0 ;
            end
            dh = dh/10  ;
           % pause


        end
    end


end

RGB = 1 - RLBm     ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DenL = rhoL ; VisL = miL ;
DenG = rhoG ; VisG = miG ;
JL = jL ; JG = jG ;
surTen = tensup ;
theta = inclinacaor ; 
DenC = DenG ; VisC = VisG ;
DenF = DenL ; VisF = VisL ;
A = pi/4 * D^2 ;

dRho = DenL - DenG ;
rug = 0 ;

    delta = altB/D ;
    lambda = 2 * acos(1 - (2*delta) ) ;

    SC = D * (pi-lambda/2) ; 
    SF = D * lambda/2 ;
    SI = D * sin(lambda/2) ;
    DC = D * ( 1 + ( sin(lambda/2)-sin(lambda)/2 ) / ...
        ( 2 * ( 2*pi-lambda*sin(lambda) ) )  )^-1 ;
    DF = D * ( 1-sin(lambda)/lambda ) ;
%     RF = ( lambda - sin(lambda) ) / (2*pi) 
    RF = 1 - RGB ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ED = funED(4, JG, JL , DenG , DenL, VisG, VisL, g, D,  surTen, dRho ) ;
    
UC =  ( JG + ED*JL ) / (1-RF) ;
UF =  ( JL - ED*JL ) / RF ;

J = JG + JL ;
Fr = J / ( dRho*g*D/DenL )^(1/2) ;
Eo = dRho*g*(D^2)/surTen  ;

if Fr<3.5
    CoT = 1+0.2*(sin(theta))^2 ;
    CinfT = (0.542-1.76/(Eo^.56))*cos(theta) + ...
        (0.345 / (1+3805*(Eo^3.06))^0.58 ) * sin(theta);
else
    CoT = 1.2 ;
    CinfT = (0.345 / (1+3805*(Eo^3.06))^0.58 ) * sin(theta);
end
VinfT = CinfT * (dRho*g*D/DenF)^(1/2) ;

UT = CoT * J + VinfT ;

ReC = DenC * UC * DC / VisC ;
ReF = DenF * UF * DF / VisF ;

CfC = ( -3.6 * log( rug/(3.7*DC)^1.11 +  6.9/ReC ) )^-2 ;
CfF = ( -3.6 * log( rug/(3.7*DF)^1.11 +  6.9/ReF ) )^-2 ;

tauWC = CfC * DenC * UC * abs(UC) / 2 ;
tauWF = CfF * DenF * UF * abs(UF) / 2 ;

TWC = tauWC * SC / A ;
TWF = tauWF * SF / A ;

end


function ED = funED(model, JG, JL , DenC , DenF, VisC, VisF, g, D,  surTen, dRho )

    
    % Wallis (Tabela4.3)
    ReL = DenF * JL * D / VisF;
    We1 = ( DenC * JG^2 * D / surTen ) * ( dRho / DenC )^(1/3);
    
    if model==1
        % Wallis (1969)
        phi = 1e4 * JG * VisC / surTen * (DenC/DenF)^0.5;
        ED = 1 - exp( -0.125 * (phi-1.5) ) ;
    elseif model==2
        % Oliemans (1986)
        Coef = (10^-2.52)*(DenC^0.18)*(DenF^1.08)*(VisC^0.28)*...
            (VisF^0.27)*(JG^1.44)*(JL^0.7)*(surTen^-1.8)*(g^0.46)*(D^1.72) ;
        ED = Coef/(Coef+1);
    elseif model == 3
        % Ishii & Mishima (1989)
        ED = tanh ( 7.25e-7 * ReL^(1/4) * We1*(5/4) ) ;
    else
        % Sawant et al. (2008)
        ReLmin = 250*log(ReL) - 1265 ;
        Em = 1 - ReLmin/ReL ;
        a = 2.31e-4 * ReL^-0.35 ;
        ED = Em * tanh (a * We1^1.25);
    end
    
    % Sawant et al. (2009)
%     Nmu = VisF / ( ( DenF^2 * surTen^3 ) / ( dRho*g ) )^(1/4) ;
%     ReLmin1 = 13*(Nmu^(-0.5)) + 0.3* (ReL-13*(Nmu^(-1/2)))^0.95  ;
%     We2 = ( DenC * JG^2 * D / surTen ) * ( dRho / DenC )^(1/4) ;
%     ED5 = ( 1-ReLmin1/ReL ) * tanh( 2.31e-4 * ReL^-0.35 * We2^(5/4) ) ;

end