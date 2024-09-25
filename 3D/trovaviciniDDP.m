function  dE = trovaviciniDDP(ddp, pdp, dpp, ddd, J,N)
%trova i vicini degli elementi pari-pari e calcola la variazione di energia

d2p  = [2:N,1];                      % 2->3, 4->5, permutazione verso sx
p2d  = [N,1:N-1];                    % 1->2N, 3->2, 5->4, permutazione verso dx



ddp_vicini = pdp + circshift(pdp, [-1 0 0]) + dpp + circshift(dpp, [0 -1 0]) + ddd + circshift(ddd, [0 0 1]);
% NN = S.eoe + S.oee + S.eoe(o2e,:,:) + S.oee(:,o2e,:) + S.ooo + S.ooo(:,:,e2o);

dE = J.*ddp.*(ddp_vicini); %variazione energia se inverto spin