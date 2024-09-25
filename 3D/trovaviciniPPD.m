function  dE = trovaviciniPPD(ppd, pdd, dpd, ppp, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



ppd_vicini = pdd + circshift(pdd, [0 1 0]) + dpd + circshift(dpd, [1 0 0]) + ppp + circshift(ppp, [0 0 -1]);
% NN = S.eoo + S.oeo + S.eoo(:,e2o,:) + S.oeo(e2o,:,:) + S.eee + S.eee(:,:,o2e);

dE = J.*ppd.*(ppd_vicini); %variazione energia se inverto spin


