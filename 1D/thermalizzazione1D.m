function [spin_pari, spin_dispari] = thermalizzazione1D(therm,spin_pari,spin_dispari,J,B,b,N)
 %processo di termalizzazione

for i = 1:therm 
   
    
    %Metropolis per gli spin pari
    vicini_disp = circshift(spin_dispari,-1);
    r = rand(1,N/2);
    dE_pari = 2*J.*spin_pari.*(spin_dispari + vicini_disp +B); %variazione energia
    prob_pari = exp(-dE_pari*b);  %la probabilità è e^-deltaE*beta
    accept = r < prob_pari;        %mossa di metropolis per tutti i pari
    %accept = r < exp(-2*b.*(spin_dispari + vicini_disp + B).*spin_pari);
    spin_pari(accept) = -spin_pari(accept);
    %exp(-2*b(cnt) .* (spin_dispari + vicini_disp + B).*spin_pari);
    
    %Metropolis per gli spin dispari
    vicini_pari = circshift(spin_pari, 1);
    
    r = rand(1,N/2);
    dE_dispari = 2*J.*spin_dispari.*(spin_pari + vicini_pari +B); %variazione energia
    prob_dispari = exp(-dE_dispari*b);
    accept = r < prob_dispari;
    %accept = r < exp(-2*b.*(spin_pari + vicini_pari + B).*spin_dispari);
    spin_dispari(accept) = -spin_dispari(accept);
    
end