function [spin, E_totale] = evoluz3D(spin, J, b, N)
%studiamo la variazione di Energia
        E_totale = 0.0;
        
        %troviamo i vicini di pari-pari-pari
        dE = trovaviciniPPP(spin.ppp, spin.pdp, spin.dpp, spin.ppd, J,N);
        
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.ppp(accept) = -spin.ppp(accept);
        
        
       %troviamo i vicini di pari-pari-dispari
        dE = trovaviciniPPD(spin.ppd, spin.pdd, spin.dpd, spin.ppp, J,N);
        
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.ppd(accept) = -spin.ppd(accept);
        
        %troviamo i vicini di pari-dispari-pari
        dE = trovaviciniPDP(spin.pdp, spin.ppp, spin.ddp, spin.pdd, J,N);
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.pdp(accept) = -spin.pdp(accept);
        
       %troviamo i vicini di pari-dispari-dispari
        dE =trovaviciniPDD(spin.pdd, spin.ppd, spin.ddd, spin.pdp, J,N);
        
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.pdd(accept) = -spin.pdd(accept);
        
       %troviamo i vicini di dispari-pari-pari
        dE =trovaviciniDPP(spin.dpp, spin.ddp, spin.ppp, spin.dpd, J,N);
        
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.dpp(accept) = -spin.dpp(accept);
        
        %troviamo i vicini di dispari-pari-dispari
        dE =trovaviciniDPD(spin.dpd, spin.ppd, spin.ddd, spin.dpp, J,N);
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.dpd(accept) = -spin.dpd(accept);
        
        %troviamo i vicini di dispari-dispari-pari
        dE =trovaviciniDDP(spin.ddp, spin.pdp, spin.dpp, spin.ddd, J,N);
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.ddp(accept) = -spin.ddp(accept);
        
        %troviamo i vicini di dispari-dispari-dispari
        dE =trovaviciniDDD(spin.ddd, spin.pdd, spin.dpd, spin.ddp, J,N);
        
        E_totale = E_totale + dE;
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.ddd(accept) = -spin.ddd(accept);
        
        
        
        
        
        
        
     