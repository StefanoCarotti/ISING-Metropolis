function  dE = trovaviciniPDP(pdp, ppp, ddp, pdd, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia


pdp_vicini = ppp + circshift(ppp, [0 -1 0]) + ddp + circshift(ddp, [1 0 0]) + pdd + circshift(pdd, [0 0 1]);


dE = J.*pdp.*(pdp_vicini); %variazione energia se inverto spin