function [Pk,DenGk,JGk,Jk,UGk,ULk]= bisectionP(x, Ro, To, GG, JL, DenL,...
    VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,Psik)
    k = 0;
    err = 100 ;
    
    a = min(x) ;
    b = max(x) ;
    
    
   
    tol = 1e-4;
    z=0;
    
    while abs(err) > tol 
        z = z+1;
        p = (a+b)/2. ;
        
        [fa,DenGk,JGk,Jk,UGk,ULk] = funP(a, Ro, To, GG, JL, DenL,...
            VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,Psik) ;

        [fb,DenGk,JGk,Jk,UGk,ULk] = funP(b, Ro, To, GG, JL, DenL,...
            VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,Psik) ;

        [fp,DenGk,JGk,Jk,UGk,ULk] = funP(p, Ro, To, GG, JL, DenL,...
            VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,Psik) ;

        %Bisection Method
        if fa*fp<0 
           b = p ;
        else
           a = p ;          
        end
        
        err = abs(b-a);
    end
    Pk = p ;
end

function [f,DenGk,JGk,Jk,UGk,ULk]=funP(xk, Ro, To, GG, JL, DenL,...
            VisG, VisL, surTen, D, S, A, theta,rug, g, dRho,PL,Psik) 

   xk-1e5;
    
    DenGk = xk / (Ro*To) ;
    
    JGk = GG / DenGk ; 
    Jk = JL + JGk;
    
    
    [alfak, TWk] = AlphaTauIntermittent_v2(JGk, JL, Jk, DenGk, DenL,...
        VisG, VisL, surTen, D, S, A,  theta,rug, g, dRho,PL,xk) ;

    
    UGk = JGk/alfak ;
    ULk = JL/(1-alfak) ;
    
    f = xk + alfak*DenGk* UGk^2 + (1-alfak)*DenL*ULk^2 - Psik;

        
%     

end