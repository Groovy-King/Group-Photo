function [absMagv]=absMag(appmag, z, d, dunit,e,k )

nargs = 6;
for k = nargin:nargs-1
    switch k
        case 0
            appmag =appmag;
        case 1
            z=0;
        case 2
            d=0;
        case 3
            unit='pc';
        case 4
            e=0;
        case 5
            k=0;
         otherwise
    end
end
pc2ly=3.261636;

        absMagv=appmag - 5 * log10(lumdistt( z,'mpc' )/(1e-05)) - k - e;

         
end

