function  dE = trovaviciniPDD(pdd, ppd, ddd, pdp, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



pdd_vicini = ppd + circshift(ppd, [0 -1 0]) + ddd + circshift(ddd, [1 0 0]) + pdp + circshift(pdp, [0 0 -1]);
%  NN = S.eeo + S.ooo + S.eeo(:,o2e,:) + S.ooo(e2o,:,:) + S.eoe + S.eoe(:,:,o2e);

dE = J.*pdd.*(pdd_vicini); %variazione energia se inverto spin