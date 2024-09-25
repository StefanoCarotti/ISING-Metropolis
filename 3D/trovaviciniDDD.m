function  dE = trovaviciniDDD(ddd, pdd, dpd, ddp, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



ddd_vicini = pdd + circshift(pdd, [-1 0 0]) + dpd + circshift(dpd, [0 -1 0]) + ddp + circshift(ddp, [0 0 -1]);
%  NN = S.eoo + S.oeo + S.eoo(o2e,:,:) + S.oeo(:,o2e,:) + S.ooe + S.ooe(:,:,o2e);

dE = J.*ddd.*(ddd_vicini); %variazione energia se inverto spin