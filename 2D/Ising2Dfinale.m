
% Montecarlo per modello di Ising in 2 dimensioni
% dichiaro le variabili
Ntot = 512;
N = Ntot/4;
J = 1;
sweeps = 2e4; %ripetizioni
therm = 1e3;
skip = 50;
bcrit = 0.44069;
% definisco l'intervallo di valori di beta per cui studiamo il sistema
beta = [0.2:0.03:0.8];  %linspace(0.2,0.8,10);
Nbeta = length(beta);


x = linspace(bcrit,0.8,128);
Ons = (1-sinh(2*x).^(-4)).^(1/8); %valore teorico di L. ONSAGER

% creo il reticolo come uno struct di 4 matrici differenti
% in una dimensione dividevamo in pari e dispari per trovare i primi vicini
% in 2 dimensione ci serve dividere in 4 matrici differenti:
% matrice dispari-dispari (1,1) (1,3) (3,1) (3,3) etc
%matrice dispari-pari (1,2) (1,4) (3,2) etc
%matrice pari-dispari (2,1) (2,3) (4,1) etc
%matrice pari-pari (2,2) (2,4) (4,2) (4,4) etc

spin.pp = sign(0.5 - rand(N,N));  %pari-pari
spin.dp = sign(0.5 - rand(N,N));  %dispari-pari
spin.dd = sign(0.5 - rand(N,N));  %dispari-dispari
spin.pd = sign(0.5 - rand(N,N));  %pari-dispari


%prealloco spazio per vettori (ottimizzazione)
M_media = zeros(1,sweeps/skip);
Mg = zeros(1,Nbeta);
errMg = zeros(1, Nbeta);

E_media = zeros(1,sweeps/skip);
Energia = zeros(1,Nbeta);
errEn = zeros(1,Nbeta);

hw = waitbar(0, "musichetta di attesa...");
for cnt = 1:Nbeta
    b = beta(cnt);
       
    for t = 1:therm    %fase di termalizzazione
        
       spin = termalizzazione2D(spin, J, b, N);
        
    end
    
    for s = 1:sweeps     %evoluzione vera e propria del sistema
       [spin, E_totale] = evoluz2D(spin, J, b, N);
       
       if rem(s,skip) == 0
            
            E_media(s/skip) = (-1/2).*sum(E_totale(:))/(4*N^2);
            
            reticolo = cell2mat(struct2cell(spin));            
            M_media(s/skip) = sum(reticolo(:));
            
        end
        
    end
    
    Mg(cnt) = abs(mean(M_media))/(4*N^2);
    errMg(cnt) = std(Mg)/sqrt(sweeps/skip-1);
    
    Energia(cnt) = mean(E_media);
    errEn(cnt) = std(Energia)/sqrt(sweeps/skip-1);
    
    waitbar(cnt/Nbeta);
end
close(hw)

%Parte dei grafici

fig1=figure; % Creiamo la prima figura
    set(fig1,'Windowstyle','docked') % Mettiamola 'docked'
    title(['Ising 2D con N = 512 ' ,...
           ' e B = 0'],...
          'FontSize',16,'FontWeight','normal','FontAngle','it') 
    grid on 
    grid minor % grigla fitta
    hold on    
    errorbar(beta,Mg,errMg,...
             'Color','green','LineStyle','-','Marker','.',...
             'MarkerSize',8) % Plot della magnetizzazione media (B)
    plot(x,Ons,'r-','LineWidth',2) % Plot della magnetizzazione scondo la formula esatta
    xline (bcrit, 'b-', 'LineWidth', 1);
    xlabel('\beta')
    ylabel('Magnetizzazione spontanea')          
    legend({'Simulazione','Previsione', '\beta critica'},'Location','NorthWest')
    hold off
    
    fig2=figure;
    set(fig2,'Windowstyle','docked') % Mettiamola 'docked'
    title('Andamento energia interna',...
          'FontSize',16,'FontWeight','normal','FontAngle','it') 
    grid on 
    grid minor % grigla fitta
    hold on    
    errorbar(beta,Energia,errEn,...
             'Color','green','LineStyle','-','Marker','.',...
             'MarkerSize',8) % Plot dell'energia interna
    xlabel('\beta')
    ylabel('Energia media per spin')          
    legend({'Simulazione'},'Location','NorthEast')