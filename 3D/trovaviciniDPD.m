function  dE = trovaviciniDPD(dpd, ppd, ddd, dpp, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



dpd_vicini = ppd + circshift(ppd, [-1 0 0]) + ddd + circshift(ddd, [0 1 0]) + dpp + circshift(dpp, [0 0 -1]);
%  NN = S.eeo + S.ooo + S.eeo(o2e,:,:) + S.ooo(:,e2o,:) + S.oee + S.oee(:,:,o2e);

dE = J.*dpd.*(dpd_vicini); %variazione energia se inverto spin