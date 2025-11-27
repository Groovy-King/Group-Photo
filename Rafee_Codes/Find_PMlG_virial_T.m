    
for vol=1:9

    for slice=1:3
            
%             dif_mass=cell(1,1);
%             dif_post=cell(1,1);
%             load(strcat('RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));
% 
%             appsolutemag=Data(:,9);
% 
%             Lum=Data(:,8);
% 
%             Xdata=Data(:,2:4);
%             [Xdata,ia,ic]=unique(Xdata,'rows','stable');
%             X_orig=Xdata;
%             lum_check=appsolutemag(ia);
%             appsolutemag=appsolutemag(ia);
%             Lum=Lum(ia);
%              lum_norm=Lum./sum(Lum);
%             [num,cent]=hist(lum_norm);
%             numbertotake=round(100*sum(num)/100);
%             [B,IX]=sort(Lum,'descend');
%             Xdata=Xdata(IX(1:numbertotake),:);
%             appsolutemag= appsolutemag(IX(1:numbertotake));
%                [Xgrid,New_idx]=check_gals( X_orig,Xdata,10^12,IX(1:numbertotake));
%              [ss,iz,ix]=unique(New_idx);
%              New_idx=New_idx(iz);
%              Xgrid=Xgrid(iz,:);
%               mass_no_final=0;
%               mass_no_final_mn=0;
%               mass_no_final_md=0;

                load(strcat('Expectedpr3_vir_mass_no_5_',num2str(vol),'_',num2str(slice),'.mat'));

                [idx_new,P_MlG,MASS_EST]=Find_PMlG_after_filter_viral_AllT(vol,slice,New_idx,Xgrid);%Find_PMlG_after_filter_viral_th(vol,slice,New_idx,Xgrid);
            
                        lia=ismember(New_idx,idx_new(:,1));
                        chosen_idx=find(lia);
                        Xgrid=Xgrid(chosen_idx,:);
                        New_idx=idx_new(:,1);
                        bufff=P_MlG;
                        for GD=1 : length(P_MlG(:,1))
                            for ms=1: length(P_MlG(1,:))
                                if sum(P_MlG(GD,:))>0
                                    bufff(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(GD,:));
                                else
                                    bufff(GD,ms)=0;
                                end
                            end
                        end
            
                        P_MlG=bufff;
                        for GD=1 : length(P_MlG(:,1))
                            for ms=1: length(P_MlG(1,:))
                                if sum(P_MlG(:,ms))>0
                                    P_GlM(GD,ms)=P_MlG(GD,ms)/sum(P_MlG(:,ms));
                                else
                                    P_GlM(GD,ms)=0;
                                end
                            end
                        end
                        
                        for gpr=1:length(MASS_EST(1,:))
                               mass_no_final_mn(gpr)=sum(P_MlG(gpr,:)'.*MASS_EST(:,gpr));
                               idx_max=find(max(P_MlG(gpr,:))==P_MlG(gpr,:));
                               mass_no_final_md(gpr)=MASS_EST(idx_max(1),gpr);
                        end
                        
                        
                        
                        
                       st=strcat('Expectedpr3_virAllT_mass_no_',num2str(5),'_',num2str(vol),'_',num2str(slice));
                        save(st, 'mass_no_final_md','P_MlG','P_GlM','mass_no_final_mn','Xgrid','New_idx','MASS_EST');
                        clearvars -except slice vol;
    end
end