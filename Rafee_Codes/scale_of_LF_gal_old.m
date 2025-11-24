function  [WLF] = scale_of_LF_gal_old(Beta,z,Flg)

h=0.7;
% faintLim=19.8;%appmagn;
faintlimflux=19.8; 
faintLim=absMag(faintlimflux,z);
if Flg

Mstar=-21.35+5*log10(h);
alpha=-1.3;
Nstar=1.49E-2*h^3;
lschechter=@(M) Nstar.* (10.^((Mstar - M)/2.5)).^(alpha + 1) .* log(10)/2.5 .* exp(-(10.^((Mstar - M)/2.5)));
        for hello=1 :length(faintLim)
            M=faintLim(hello);
            flux(hello)= integral(lschechter,-25,M); 
        end
    LF=flux;
    WLF=(LF(:).^(Beta))./sum((LF).^(Beta));
else
        
        Mstar=-20.44+5*log10(h);
        alpha=-1.05;
        Nstar=1.49E-2*h^3;
        lschechter=@(M) Nstar.* (10.^((Mstar - M)/2.5)).^(alpha + 1) .* log(10)/2.5 .* exp(-(10.^((Mstar - M)/2.5)));
        for hello=1 :length(faintLim)
            M=faintLim(hello);
            flux(hello)= integral(lschechter,-25,M); 
        end
    LF=flux;
    WLF=(LF(:).^(Beta))./sum((LF).^(-Beta));
    
end
end

