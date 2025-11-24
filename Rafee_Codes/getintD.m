function [ res ] = getintD(z_min,z_max) 
    
    cosmopar;
    n=1000;
    res=[];
    buf=z_max;
    for k=1 : length(buf)
        z_max=buf(k);
        dz=(z_max-z_min)/n;
        zz=[z_min:dz:z_max]; 
        E=sqrt(omg_m*(1+zz).^3+omg_k*(1+zz).^2+omg_lambda);
        f=1./E;
        res=[res sum(f.*dz)];
    end
end

