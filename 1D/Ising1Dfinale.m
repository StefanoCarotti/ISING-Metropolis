
function [setup] = Ising1Dfinale(setup)
% Montecarlo che simula modello di Ising in una dimensione
tic;

entries={'Dimensione del reticolo N' ...
            'Costante di accoppiamento J'...
            'Campo mangnetico B'...
            'Cicli di evoluzione'...
            'Cicli di termalizzazione'...
            'Passo di campionamento'... % "skip"
            
        }; 

   if nargin<1
      default={
          '512'...
          '1'  ...
          '0.001'...
          '1.0e+5'...
          '1.0e+3'...
          '200'...
          
              };
   else
      default=setup;                      
   end

   addopt.Resize='on';
addopt.WindowStyle='normal';
addopt.Interpreter='tex';
header = 'Algoritmo di Metropolis per il Modello di Ising 1D ';
lines = 1;
setup = inputdlg(entries, header, lines, default, addopt);


N = str2double(setup{1});
J = str2double(setup{2});
B = str2double(setup{3});
sweeps = str2double(setup{4});
therm = str2double(setup{5});
skip = str2double(setup{6});

[Magn, errMag, beta] =  monte1D(N,J,B,sweeps,therm,skip);



%momento grafici
ms = sinh(B.*beta)./sqrt(sinh(B.*beta).^2+exp(-4.*beta)); % magnetizzazione teorica

graf1=figure; % prima figura
    set(graf1,'Windowstyle','docked') % Mettiamola 'docked'
    title(['Ising1D Metropolis - Magnetizzazione media - N = 512 ',...
           '- B = 0'],...
          'FontSize',16,'FontWeight','normal','FontAngle','it') 
    grid on 
    grid minor % grigla fitta
    hold on    
plot(beta,ms,'r-','LineWidth',2) % Plot della magnetizzazione scondo la formula esatta
    xlabel('\beta')
    ylabel('Magnetizzazione') 
    ylim([-0.1 1.2])
   
    
    errorbar(beta,Magn,errMag./sqrt(sweeps/skip),...
             'Color','blue','LineStyle','-','Marker','p',...
             'MarkerSize',4)
%     plot(beta,Magn,...
%              'Color','blue','LineStyle','-','Marker','p',...
%              'MarkerSize',8)
    if B < 0
       legend({'Simulazione','Previsione'},'Location','NorthEast') %
       % Mettiamo la legenda in modo che non intralci il diagramma
    else
       legend({'Simulazione','Previsione'},'Location','NorthWest') 
    end
tempo = toc;