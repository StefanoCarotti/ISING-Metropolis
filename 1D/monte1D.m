function [Magn, errMag, beta] =  monte1D(N,J,B,sweeps,therm,skip)
%montecarlo per ising1D
beta = linspace(0.2,5,50);
Nbeta = length(beta);

% Pre-allocazione degli oggetti richiesti per velocizzare il tempo di computazione
Magmedia = zeros([1,floor(sweeps/skip)]); % Definizione del vettore con lunghezza opportuna
Magn     = zeros(size(beta)); % inizializzazione del vettore
errMag   = zeros(size(beta)); % inizializzazione del vettore


spin_pari    = sign(0.5 -rand(1,N/2)); %creo N/2 spin pari
spin_dispari = sign(0.5 -rand(1,N/2)); %creo N/2 spin dispari

hw    = waitbar(0,'Musichetta di attesa...'); 

for cnt = 1:Nbeta
    b = beta(cnt);    
    [spin_pari, spin_dispari] = thermalizzazione1D(therm,spin_pari,spin_dispari,J,B,b,N);
    %Inizio evoluzione sistema
    
    [spin_pari, spin_dispari, Magmedia] = evoluz1D(sweeps, spin_pari, spin_dispari,Magmedia,J,B,b,skip,N); 
   
    Magn(cnt) = mean(Magmedia);         %magnetizzazione media per spin
    errMag(cnt) = std(Magmedia);        %errore sulla magnetizzazione media
    
    waitbar(cnt/Nbeta);
    %dE = 2*J*spin(row,col)*sum(vicini); 
end
close(hw);