function spin = termalizzazione2D(spin, J, b, N)
%troviamo i vicini di pari-pari
        dE = trovaviciniPariPari(spin.pp, spin.dp, spin.pd, J,N);
        
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        R = rand(N,N);
        accept = R < prob;
        spin.pp(accept) = -spin.pp(accept);
        
        %troviamo i vicini di dispari-pari
        dE = trovaviciniDispariPari(spin.dp, spin.pp, spin.dd, J,N);
        
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.dp(accept) = -spin.dp(accept);
        
        %troviamo i vicini di dispari-dispari
        dE = trovaviciniDispariDispari(spin.dd, spin.pd, spin.dp, J,N);
        
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.dd(accept) = -spin.dd(accept);
        
        %troviamo i vicini di pari-dispari
        dE =trovaviciniPariDispari(spin.pd, spin.dd, spin.pp, J,N);
        
        R = rand(N,N);
        prob = exp(-2*dE*b);      %Probabilità di invertire spin
        
        accept = R < prob;
        spin.pd(accept) = -spin.pd(accept);