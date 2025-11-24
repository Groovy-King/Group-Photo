function [ energy ] = E( z )

 cosmopar;
energy=sqrt(omg_m.*(1+z).^3+omg_k.*(1+z).^2+omg_lambda);

end

