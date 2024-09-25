function  dE = trovaviciniDispariDispari(dispari_dispari, pari_dispari, dispari_pari, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx

%dd_vicini = pari_dispari + pari_dispari(d2p,:) + dispari_pari + dispari_pari(:,d2p);

dd_vicini = pari_dispari + circshift(pari_dispari, [0 -1]) + dispari_pari + circshift(dispari_pari, [-1 0]);
          
      dE = J.*dispari_dispari.*(dd_vicini); %variazione energia se inverto spin