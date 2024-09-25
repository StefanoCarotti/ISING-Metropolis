function  [spin_pari, spin_dispari, Magmedia] = evoluz1D(sweeps, spin_pari, spin_dispari,Magmedia,J,B,b,skip,N) 
%evoluzione del sistema con calcolo magnetizzazione media in base al passo
%di campionamento (skip)

for s = 1:sweeps
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
    if rem(s,skip) == 0
        Magmedia(s/skip) = (sum(spin_pari)+sum(spin_dispari))/(N);
    end
    
    
end