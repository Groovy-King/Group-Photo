% matlabpool('open','AttachedFiles',{'E.m'})

disp('Starting Run!')
for vol=1:9
    for slice=1:3
            cosmopar;
            load(strcat('Input/RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));
            Xdata=Data(:,2:4);%[ra,dec,z];
            [Xdata,ia,ic]=unique(Xdata,'rows','stable');
            
            Pos=Xdata;

            Set=cell(1,1);
            
            massspace=logspace(10,15,100);
            idxidx=round(linspace(41,90,10));

            q=0.999999; %the tolarance for the radial dispersion
            
            for masloop=1:length(idxidx)
                    mass=massspace(idxidx(masloop));
                       for k=1:length(Pos)
                            mu=Pos(k,:);
                            [current_data,ia]=setdiff(Xdata,Pos(k,:),'rows','stable');
                            rho=mu(3);
                            sss=0.26*log(10)/3;
                            sig_v2=((1.55e8/10^0.0403)*(mass*E(rho)/10^14)^0.94)^(1/3)*exp(sss/2);
                            sig_z=((1+rho)*sig_v2/c);
                            [ idx_flag ] = isittrue_test(mu,current_data,sig_z,mass,q);
                            if length(idx_flag)>3
                                buf=[mu ;current_data(idx_flag,:)];
                             Set{k}=buf;
                            else
                             Set{k}=[];
                            end
                       end
                       save(strcat('Output_Sets/Set_mass_vol',num2str(vol),'_slice_',num2str(slice),'_',num2str(idxidx(masloop))),'Set');
                       Set=cell(1,1);
            end

                clearvars -except slice vol
    end
end
% matlabpool close;
