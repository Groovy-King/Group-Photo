function  [H_z,D_a]=getHzDa(z)

     cosmopar; %invoke the cosmological parameters
     E=sqrt(omg_m*(1+z).^3+omg_k*(1+z).^2+omg_lambda);
     H_z=H0*E;
     D_c=(c/H0).*getintD(0,z) ;% Comoving distance, c speed of light
     D_a=(D_c'./(1+z)); %Angular Diameter Distance

end

