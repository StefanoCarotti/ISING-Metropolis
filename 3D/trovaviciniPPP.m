function  dE = trovaviciniPPP(ppp, pdp, dpp, ppd, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



ppp_vicini = pdp + circshift(pdp, [0 1 0]) + dpp + circshift(dpp, [1 0 0]) + ppd + circshift(ppd, [0 0 1]);
%NN = S.eoe + S.oee + S.eoe(:,e2o,:) + S.oee(e2o,:,:) + S.eeo + S.eeo(:,:,e2o);

dE = J.*ppp.*(ppp_vicini); %variazione energia se inverto spin


