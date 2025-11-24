function [ len_prp, r500] =NFW_cylinder_find_quantile_old(q,MAS,H_z,Da,zcenter)

                                %Initialization
                                    c200=(6.71*(MAS/(2*10^14))^(-.091) *(1+zcenter)^(-0.44))/2;
                                    c500=0.6675*c200;
                                    cosmopar; %invoke the cosmological parameters
                                    r_min=0.001;    %The minimum radius of the cluster
                                    rho_crit=(3*(H_z^2))/(8*pi*G); %The critical density of the universe at a redshift(z)
                                    r500=(MAS/((4/3)*pi*500*rho_crit))^(1/3);
                                    r_max=1.5*r500;   %The maximum radius of the cluster
                                    r=linspace(r_min,r_max,1000);
                                    phi=r/Da;%phys2angl(r,zcenter,'radians');%r/Da;

                                    r_s=r500/c500;
                                    
                                %NFW
                                    phi_s=r_s/Da;
                                    x=r./r_s;
                                    xphi=phi/phi_s; 
                                    
                                     r200=r500/0.67;
%                                     c200=r200/r_s;
                                    delta_c=(200/3)*(c200^3)/(log(1+c200)-(c200/(1+c200)));
                                    rho_s= rho_crit*delta_c;
                                    
%                                     Sigmaphi=getsigma(x',r_s,rho_s);
%                                      dr=r(3)-r(2);
%                                     Lphi=Sigmaphi/sum(2*pi*dr.*r'.*Sigmaphi);

                                    Sigmaphi=getsigma(xphi',phi_s,rho_s);
                                    dphi=phi(3)-phi(2);
                                    Lphi=Sigmaphi/sum(2*pi*dphi.*phi'.*Sigmaphi);
                                    
                                %Normalization    
                                    Prphi=Lphi/sum(Lphi);
%                                plot(phi,Prphi);
                                                                                                       
                                    %%%%%
                                     for k=1:length(Prphi)
                                             comp=sum(Prphi(1:k));
                                             if q<=comp
                                                  len_prp=phi(k);
%                                                     len_prp1=phys2angl(r(k),zcenter,'radians');
                                                 break;
                                             end
                                        
                                     end
          
end

