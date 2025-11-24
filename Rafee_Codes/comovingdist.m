function Comovdis = comovingdist( z,unit )

 [Omega_m,Omega_l,Omega_k,~,H0,c,pc2ly,~]=Collect;
nargs = 2;
for k = nargin:nargs-1
    switch k
        case 0
            z =z;
        case 1
            unit='mpc';
         otherwise
    end
    
end
s=1;
 
        if strcmp(unit, 'mpc')
        s=1;
        elseif strcmp(unit, 'kpc')
        s=1000;
        elseif strcmp(unit, 'pc')
                s=1e+06;
        elseif strcmp(unit, 'ly') 
                s=pc2ly * 1e+06;
        end
                    
     
    for hello=1 :length(z)
    Ez=@ (x) 1./(sqrt(Omega_m*(1+x).^3+Omega_k*(1+x).^2+Omega_l)); %For specific z
%     Fun= 1/Ez;
    x=z(hello);
    Ezz= integral(Ez,0,x); %min , max
    Comdis(hello)=s*(c/H0) * Ezz;  %%%find the comoving distance
    end
    
    Comovdis=Comdis;
   end
   



% comoving.dist <- function (z, unit = "mpc") {
%     if (!(unit %in% c("pc", "kpc", "mpc", "ly", "m"))) {
%         message("Unknown physical unit. Please use one of either pc, kpc, mpc, or ly")
%         stop()
%     }
%     ifelse(unit == "mpc", 1,
%         ifelse(unit == "kpc", 1000,
%             ifelse(unit == "pc", 1e+06,
%                 ifelse(unit == "ly", const$pc2ly * 1e+06, const$mpc2m)
%             )
%         )
%     ) 
%     * const$c.kms/const$H0 * integrate(function(x) 1/E.z(x), lower = 0, upper = z)$value
% }

