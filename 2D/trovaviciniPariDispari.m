function  dE = trovaviciniPariDispari(pari_dispari, dispari_dispari, pari_pari, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx

%pd_vicini = dispari_dispari + dispari_dispari(p2d,:) + pari_pari + pari_pari(:,d2p);
pd_vicini = dispari_dispari + circshift(dispari_dispari, [0 1]) + pari_pari + circshift(pari_pari, [-1 0]);
          
      dE = J.*pari_dispari.*(pd_vicini); %variazione energia se inverto spin