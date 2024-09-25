function  dE = trovaviciniPariPari(pari_pari, dispari_pari, pari_dispari, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

pp_vicini = dispari_pari + circshift(dispari_pari, [0 1]) + pari_dispari + circshift(pari_dispari, [1 0]);


dE = J.*pari_pari.*(pp_vicini); %variazione energia se inverto spin


