function   prob = Noise_modified_Disc_adjusted2_old(XYZd,GP,mass)%(GP,XYZd,mass)  %phi_s
            
                cosmopar
                [H_z,Da]=getHzDa(XYZd(3));
                XYZd=[XYZd(1)*(pi/180) XYZd(2)*(pi/180) XYZd(3)];
                XYZd=XYZd';
                len_i=sqrt(XYZd(1)^2+XYZd(2)^2+XYZd(3)^2);
                u_i=XYZd./len_i;
                GP=[GP(:,1)*(pi/180) GP(:,2)*(pi/180) GP(:,3)];
                XYZg=GP';
                len_g=sqrt(GP(:,1).^2+GP(:,2).^2+GP(:,3).^2);
                u_g=XYZg(1:3,:);
                u_g(1,:)=u_g(1,:)./len_g(:)';
                u_g(2,:)=u_g(2,:)./len_g(:)';
                u_g(3,:)=u_g(3,:)./len_g(:)';
                
                phi=acos(dot(u_g,repmat(u_i,1,length(u_g(1,:)))));
                X_t=repmat(XYZd,1,length(XYZg(1,:)))-XYZg;
                 a_z=zeros(length(X_t(1,:)),1);
                
                 
                 z_len=abs(XYZg(3,:)-repmat(XYZd(3),1,length(XYZg(3,:))));
          
                for k=1:length(X_t(1,:))
                    a_z(k)=abs(u_g(1:3,k)'*X_t(1:3,k));
                end
          
          
                c= 299792.458 ; %Km/S Speed of light



                                sss=0.26*log(10)/3;
                                sig_v2=((1.55e8/10^0.0403)*(mass*E(XYZd(3))/10^14)^0.94)^(1/3)*exp(sss/2);
                                sig_z=(1+XYZd(3))*(sig_v2/c);
                                idx_zlocal=find(z_len<=(2*sig_z));

                                den=1./(sig_z.*sqrt(2.*pi));
                                buf=(den.*exp(-((a_z(idx_zlocal)).^2)./((2*sig_z.^2))));
                                
                                prob1=zeros(1,length(a_z));
                                prob1(idx_zlocal)=buf;
                                prob1=prob1';
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                    rho_crit=(3*(H_z.^2))/(8*pi*G); 
                                    r500=(mass./((4/3)*pi*500.*rho_crit)).^(1/3);
                                    
                                    r=real(phi.*Da');
                                    cutoff=1.5*r500;
                                    idx_rlocal=find(r<=cutoff & r~=0);
                                    r=r(idx_rlocal);
                                    if ~isempty(r) && length(r)>1 
                                            c200=(6.71*(mass/(2*10^14))^(-.091) *(1+XYZd(3))^(-0.44))/2;
                                            c500=0.6675*c200;
                                            r_s=r500./c500;
                                            phi_s=r_s;
                                            xphi=r'./phi_s; 
                                            r200=r500./0.67;
                                            delta_c=(200./3).*(c200.^3)./(log(1+c200)-(c200./(1+c200)));
                                            rho_s= rho_crit.*delta_c;
                                            Sigmaphi=getsigma(xphi,phi_s,rho_s);
                                            prob2=(1:length(GP(:,1)))*0;
                                            prob2(idx_rlocal)=Sigmaphi./sum(2*pi*r'.*Sigmaphi);
                                    else
                                        prob2=(1:length(GP(:,1)))*0;
                                    end
                                  
                prob=prob1.*prob2';
        
end

