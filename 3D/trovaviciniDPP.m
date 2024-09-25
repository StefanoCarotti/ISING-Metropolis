function  dE = trovaviciniDPP(dpp, ddp, ppp, dpd, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



dpp_vicini = ddp + circshift(ddp, [0 1 0]) + ppp + circshift(ppp, [-1 0 0]) + dpd + circshift(dpd, [0 0 1]);
%  NN = S.ooe + S.eee + S.ooe(:,e2o,:) + S.eee(o2e,:,:) + S.oeo + S.oeo(:,:,e2o);

dE = J.*dpp.*(dpp_vicini); %variazione energia se inverto spin