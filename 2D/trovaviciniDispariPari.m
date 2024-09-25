function  dE = trovaviciniDispariPari(dispari_pari, pari_pari, dispari_dispari, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx

%dp_vicini = pari_pari + pari_pari(d2p,:) + dispari_dispari + dispari_dispari(:,p2d);

dp_vicini = pari_pari + circshift(pari_pari, [0 -1]) + dispari_dispari + circshift(dispari_dispari, [1 0]);
          
      dE = J.*dispari_pari.*(dp_vicini); %variazione energia se inverto spin