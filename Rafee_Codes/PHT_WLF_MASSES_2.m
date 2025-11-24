% P(G|M) based on P(M|G)
%The correct one


% matlabpool ;
for vol=1:9

    for slice=1:3
            
            dif_mass=cell(1,1);
            dif_post=cell(1,1);
            load(strcat('RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));
%             run('Find_the_whole_sets_11.m')
            
            %Take the absolute mag.
            appsolutemag=Data(:,9);
            
            %Take the gals luminusity 
            Lum=Data(:,8);
            
            %find the unique gals
            Xdata=Data(:,2:4);
            [Xdata,ia,ic]=unique(Xdata,'rows','stable');
            X_orig=Xdata;
            
            %filter out the gal luminusities and the absolute mag
%             Lum=Lum(ia);
            lum_check=appsolutemag(ia);
            appsolutemag=appsolutemag(ia);
            Lum=Lum(ia);
            
            lum_norm=Lum./sum(Lum);
            [num,cent]=hist(lum_norm);
            numbertotake=round(100*sum(num)/100);
            [B,IX]=sort(Lum,'descend');
            Xdata=Xdata(IX(1:numbertotake),:);
            appsolutemag= appsolutemag(IX(1:numbertotake));
             
             
            % filter out the significan position based on the neighbourhood
            % richness (i.e. if there is very poor neighbourhood I will
            % neglect the particular galaxy to be potential galaxy group centre)
            % the output are the expected potential galaxy group centres and
            % their indices in the original data
            
            %Find the weight for the expected potential galaxy groups

            
            Beta=-0.4;
            
%             [Xgrid,New_idx]=check_gals_lum( X_orig,Xdata,10^12,IX(1:numbertotake),lum_check);
%             save('newnew_lum','Xgrid','New_idx');
%                 load('newnew_lum.mat');

%              [Xgrid,New_idx]=check_gals( X_orig,Xdata,10^12,IX(1:numbertotake));
%             save('newnew','Xgrid','New_idx');
%             load('newnew.mat')

             
%              [ss,iz,ix]=unique(New_idx);
%              New_idx=New_idx(iz);
%              Xgrid=Xgrid(iz,:);


            load(strcat('Expectedpr3_virAllT_mass_no_5_',num2str(vol),'_',num2str(slice),'.mat'));

             %Group WLF=Flg=1
            Flg=1;
%             appsolutemag_group=appsolutemag(New_idx);
            [WLF1]=scale_of_LF_gal_old(Beta,Xdata(:,3),Flg);
%             Flg=1;
%             [lia,~]=ismember((1:length(appsolutemag))',New_idx);
%             idd=find(lia==0);
%             appsolutemag_bg=appsolutemag(idd);
%             [WLF2]=scale_of_LF_gal_old(Beta,appsolutemag_bg,Flg);
%             WLF=appsolutemag;
%             WLF(idd)=WLF2;
%             WLF(New_idx)=WLF1;
            WLF=WLF1;%./sum(WLF);
            
            mass_no_final=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             [P_MlG,idx_new]=Find_PMlG_after_filt(vol,slice,New_idx); 
% % %             lia=ismember(New_idx,idx_new(:,1));
% % %             chosen_idx=find(lia);
% % %             Xgrid=Xgrid(chosen_idx,:);
% % %             New_idx=idx_new(:,1);
% % %             bufff=P_MlG;
% % %             for GD=1 : length(P_MlG(:,1))
% % %                 for ms=1: length(P_MlG(1,:))
% % %                     if sum(P_MlG(GD,:))>0
% % %                         bufff(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(GD,:));
% % %                     else
% % %                         bufff(GD,ms)=0;
% % %                     end
% % %                 end
% % %             end
% % %             
% % %             P_MlG=bufff;
% % %             ms_bands=logspace(10,15,100);
% % %             idxidx=round(linspace(41,90,10));
% % %             ms_bands=ms_bands(idxidx);
% % %             mass_no_final=ms_bands*P_MlG';
% % %             st=strcat('ExpectedT1_mass_no_',num2str(5),'_',num2str(vol),'_',num2str(slice));
% % %             save(st, 'mass_no_final');
% % %                                 
% % %                                 
% % %             for GD=1 : length(P_MlG(:,1))
% % %                 for ms=1: length(P_MlG(1,:))
% % %                     if sum(P_MlG(:,ms))>0
% % %                         P_GlM(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(:,ms));
% % %                     else
% % %                         P_GlM(GD,ms)=0;
% % %                     end
% % %                 end
% % %             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% mass_virial
%                         [idx_new,P_MlG,MASS_EST]=Find_PMlG_after_filter_viral_th(vol,slice,New_idx,Xgrid);
%             
%                         lia=ismember(New_idx,idx_new(:,1));
%                         chosen_idx=find(lia);
%                         Xgrid=Xgrid(chosen_idx,:);
%                         New_idx=idx_new(:,1);
%                         bufff=P_MlG;
%                         for GD=1 : length(P_MlG(:,1))
%                             for ms=1: length(P_MlG(1,:))
%                                 if sum(P_MlG(GD,:))>0
%                                     bufff(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(GD,:));
%                                 else
%                                     bufff(GD,ms)=0;
%                                 end
%                             end
%                         end
%             
%                         P_MlG=bufff;
%                         for GD=1 : length(P_MlG(:,1))
%                             for ms=1: length(P_MlG(1,:))
%                                 if sum(P_MlG(:,ms))>0
%                                     P_GlM(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(:,ms));
%                                 else
%                                     P_GlM(GD,ms)=0;
%                                 end
%                             end
%                         end
%                         
%                         for gpr=1:length(MASS_EST(1,:))
%                                mass_no_final(gpr)=sum(P_MlG(gpr,:)'.*MASS_EST(:,gpr));
%                         end
%                         st=strcat('ExpectedT_mass_no_',num2str(5),'_',num2str(vol),'_',num2str(slice));
%                         save(st, 'mass_no_final');
%                         

        

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            post=zeros(length(Xgrid),1);
            mass_val=logspace(10,15,100);
            idxidx=round(linspace(41,90,10));
            mass_val=mass_val(idxidx);
            H_total=zeros(length(Xgrid),1);

            for mass=1:length(idxidx)     
%             
                        Likelihood=[];
                        Posterior1=[];
                            
                            
                
                            P_g_i_m_500_sum=zeros(length(X_orig),1);
                            prior_G=P_GlM(:,mass)';
                

                            for k=1:length(Xgrid)
                                %mass_virial
                                P_g_i_m_500=Noise_modified_Disc_adjusted2_old(Xgrid(k,:),X_orig,MASS_EST(mass,k)).*prior_G(k); % the prior to this k
                                
%                                 P_g_i_m_500=Noise_modified_Disc_adjusted2_old(Xgrid(k,:),X_orig,mass_val(mass)).*prior_G(k); % the prior to this k
                                P_g_i_m_500_sum=P_g_i_m_500_sum+P_g_i_m_500;
                            end

                            %% Find non zeros



                            for k=1:length(Xgrid)            
                                        %mass_virial
                                        
                                        buf_temp_prep=((Noise_modified_Disc_adjusted2_old(Xgrid(k,:),X_orig,MASS_EST(mass,k)).*prior_G(k))./P_g_i_m_500_sum);
                                        
%                                         buf_temp_prep=((Noise_modified_Disc_adjusted2_old(Xgrid(k,:),X_orig,mass_val(mass)).*prior_G(k))./P_g_i_m_500_sum);
                                        fnan=isnan(buf_temp_prep);
                                        buf_temp_prep(fnan==1)=0;
                                        post(k,:)=sum(buf_temp_prep.*WLF) ; % the prior to this k
                                        %X_orig is all gals


                            end

                      
                                P_M_l_g=P_MlG(:,mass);
                                H_total=H_total+(post.*P_M_l_g);
            end
            
            gals_pos=X_orig;
            Expected_group_idx=New_idx;
            save(strcat('Prob_paperAllT_',num2str(vol),'_',num2str(slice)),'gals_pos','Expected_group_idx','-v7.3');
            GP=Xgrid;
            Hough=H_total;
            save(strcat('PHT_paperAllT_',num2str(vol),'_Slice_',num2str(slice),'.mat'),'GP','Hough');
            clearvars -except slice vol;
            H_total=[];
            Xgrid=[];
            post=[];
            P_g_i_m_500_sum=[];
            
            

    end
end

ans = 1