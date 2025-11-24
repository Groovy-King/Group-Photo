function [ lumdist ] = lumdistt( z,unit )
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
    covdis=comovingdist( z,unit )';
    
    lumdist=(1+z).*covdis;
end

% lum.dist <- function (z, unit = "mpc") {
%     (1 + z) * comoving.dist(z, unit)
% }
