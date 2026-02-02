function [idx_new,p,MASS_EST] = Find_PMlG_after_filter_viral_All(vol,slice,idx,New_idx)


                        %%%%%%%%%%%%%%%%%
%                          All_data=[];
%                                 for all_cone_vol=1: 9
%                                     for all_cone_slice=1:3
%                                         load(strcat('RA_Dec_Z_Vol_',num2str(all_cone_vol),'_Slice',num2str(all_cone_slice),'.mat'));
%                                         Xdataall=Data(:,2:4);
%                                         [Xdatall,ia,ic]=unique(Xdataall,'rows','stable');
%                                         All_data=[All_data ; Xdatall repmat(all_cone_vol,1,length(Xdatall(:,1)))' repmat(all_cone_slice,1,length(Xdatall(:,1)))'];
%                                         
%                                     end
%                                 end

                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               idx_new=[];
                                load(strcat('Input/RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));
                                X_orig=Data;
                                Xdata=Data(:,2:4);
                                [Xdata,ia,ic]=unique(Xdata,'rows','stable');
                                X_orig=X_orig(ia,:);
                                
                                Lum=X_orig(:,8);
                                [~,IX]=sort(Lum,'descend');
                                Xdata=Xdata(IX,:);
%                                 load('New_idx_file.mat')
%                                 [New_idx]=check_gals(X_orig(:,2:4),Xdata,10^12,IX);
%                                 save('New_idx_file','New_idx');
                                indexx=round(linspace(1,31,10));
                                idxidx=round(linspace(41,90,10));            

                                 tolar=1;
                                 
                                
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 %%%%%%%%%%%%%%%%%
                                 
                                 p=cell(1,1);     
                                
                                for pos=1:length(idx)
                                        
                                             idx_rf=find(Xdata(:,3)>(Xdata(idx(pos),3)-(tolar/1000)) & Xdata(:,3)<(Xdata(idx(pos),3)+(tolar/1000)));
                                             
                                             
                                             
                                             
%                                              idx_rf=find(All_data(:,3)>(Xdata(idx(pos),3)-(tolar/1000)) & All_data(:,3)<(Xdata(idx(pos),3)+(tolar/1000)));
                                             
                                            %%Filter idx_rf if it is group
                                            %%or likely to be will remain
                                            %%otherwise remove it from the
                                            %%list
                                            
                                            
                                            
                                            
                                            if ~isempty( idx_rf)
%                                                                     ididid=find(ismember(Xdata(idx_rf,:),New_idx,'rows'));
%                                                                     idx_rf=idx_rf(ididid);
                                                                    new_list=[];
                                                                    
                                                                    while ~isempty(idx_rf)

                                                                        %check how many galaxies
                                                                        %around and check the
                                                                        %gaussainity of them if
                                                                        %they are okay then keep
                                                                        %the position
                                                                        XYZd=Xdata(idx_rf(1),:);
                                                                        GP=Xdata;
                                                                        c= 299792.458 ; %Km/S Speed of light
                                                                        sss=0.26*log(10)/3;
                                                                        mass=10^(12);
                                                                        sig_v2=((1.55e8/10^0.0403)*(mass*E(XYZd(3))/10^14)^0.94)^(1/3)*exp(sss/2);
                                                                        sig_z=(1+XYZd(3))*(sig_v2/c);
                                                                        [H_z,Da]=getHzDa(XYZd(3));

                                                                        
                                                                        index1=[GP(:,1)*(pi/180) ,GP(:,2)*(pi/180), GP(:,3)];
                                                                        mu1=[XYZd(1)*(pi/180) XYZd(2)*(pi/180) XYZd(3)]; 
                                                                        cor=index1;
                                                                        vector_from_center_to_point=repmat(mu1,length(cor(:,1)),1)-cor; 
                                                                        Y = sqrt(sum(bsxfun(@minus, repmat(mu1,length(cor(:,1)),1),cor).^2,2));
                                                                        len_u=sqrt(mu1(1)^2+mu1(2)^2+mu1(3)^2);
                                                                        unit_vec_alongtheline_ofsighrt=mu1./len_u;   
                                                                        
                                                                        vec_alongtheline_ofsighrt=0;
                                                                        normal_vector=0;
                                                                        radian=0;
                                                                        
                                                                        for dot_vec=1:length(cor(:,1))
                                                                            vec_alongtheline_ofsighrt(dot_vec)=dot(vector_from_center_to_point(dot_vec,:),unit_vec_alongtheline_ofsighrt');
                                                                            normal_vector(dot_vec)=norm(vector_from_center_to_point(dot_vec,:)-(unit_vec_alongtheline_ofsighrt*vector_from_center_to_point(dot_vec,:)'*unit_vec_alongtheline_ofsighrt));
                                                                            radian(dot_vec)=sqrt(Y(dot_vec).^2-vec_alongtheline_ofsighrt(dot_vec).^2);
                                                                        end
                                                                        [expected_radius]=NFW_cylinder_find_quantile_old(0.99999,mass,H_z,Da,XYZd(3));
                                                                        f=normal_vector<=expected_radius;
                                                                        s=vec_alongtheline_ofsighrt>=(-2*sig_z);
                                                                        t=vec_alongtheline_ofsighrt<=(+2*sig_z);
                                                                        Flag=find(f==1 & s==1 & t==1);
                                                                        if length(Flag)>=3
                                                                            new_list=[new_list idx_rf(1)];
                                                                            idx_rf=idx_rf(2:end);
                                                                        else
                                                                            idx_rf=idx_rf(2:end);
                                                                        end

                                                                    end
                                                                    idx_rf=new_list;
                                            end
                                             
                                          
                                             
                                            if ~isempty( idx_rf)
                                                                Q=cell(1,1);
                                                                numelem=cell(1,1);
                                                                xcentr=cell(1,1);
                                                                
                                                                for curr_ms=1:length(indexx)

                                                                                mass=idxidx(curr_ms);
									%	fprintf('Loading Output Sets %d\n', num2str(mass));
                                                                                load(strcat('Output_Sets/Set_mass_vol',num2str(vol),'_slice_',num2str(slice),'_',num2str(mass),'.mat'));
                                                                                Set_b=Set(idx_rf);
                                                                                non=cellfun(@isempty,Set_b);
                                                                                idxnonz=find(non==0);
                                                                                Set_b=Set_b(idxnonz);
                                                                                
                                                                                   %%%%%%%%%%%%%%%%%
                                                                                   Flagidx=[];
                                                                            for expenum=1:length(Set_b)
                                                                        
                                                                            %%check the gaussainity
                                                                                     XX=Set_b{expenum};
                                                                                    v=XX(:,3)*c;
                                                                                    [H, pValue, W] = swtest(v,0.1);
                                                                                    if pValue>0.1 && H==0%pValue>0.1
                                                                                         Flagidx=[Flagidx ; expenum];
                                                                                    end
                                                            
                                                 
                                                                            end
                                             
                                                                            %%%%%%%%%%%
                                                                                Set_b=Set_b(Flagidx);
                                                                                                                                                                
                                                                                Max_C=max(cellfun(@length,Set_b));
                                                                                Min_C=min(cellfun(@length,Set_b));

                                                                                %Smoothing before normalazation //doesn't woork
                                                                                if ~isempty(Max_C)
                                                                                            All_c_lengths=(cellfun(@length,Set_b));
                                                                                            [numel,ed,~]=histcounts(All_c_lengths);
                                                                                            edgcenter=0;
                                                                                                for kkk=2:length(ed)
                                                                                                    edgcenter(kkk-1)=(ed(kkk-1)+ed(kkk))/2;
                                                                                                end
                                                                                                xcen=edgcenter;
                                                                                            

                                                                                            idxhist=find(numel~=0);
                                                                                            numelem{curr_ms,:}=numel(idxhist);
                                                                                            xcentr{curr_ms,:}=xcen(idxhist);
                                                                                           h_new=zeros(1,length(xcentr{curr_ms,:}));
                                                                                           

                                                                                        h_new=smooth(numel(idxhist),'rloess');
                                                                                        if h_new==0
                                                                                            Q{curr_ms,:}=0;
                                                                                        else
                                                                                           Q{curr_ms,:}=h_new(:)./sum(h_new);
                                                                                        end
                                                                                        
                                                                                else
                                                                                    Q{curr_ms,:}=[];
                                                                                    xcentr{curr_ms,:}=[];    
                                                                                end

                                                                end

                                                                buf_idx=cell(1,1);
                                                                coun=1;

                                                                    
                                                                    %find the closest val index
                                                                    
                                                                    
                                                                    for u=1:length(Q)
                                                                     for v=1:length(Q{u}(:))
                                                                         
                                                                         val=xcentr{u}(v);
                                                                         idx_tosave=0;
                                                                         for u_cop=1:length(Q)
                                                                            temp_buf=abs(repmat(val,length(Q{u_cop}(:)),1)-xcentr{u_cop}(:));
                                                                            if ~isempty(temp_buf)
                                                                                 ind_cand=find(min(temp_buf)==temp_buf);
                                                                                 idx_tosave(u_cop)=ind_cand(1);
                                                                            end
                                                                         end
                                                                         buf_idx{coun,:}=idx_tosave;
                                                                         coun=coun+1;
                                                                     end
                                                                    end
                                                                    
                                                                    coun=1;
                                                                    for u=1:length(Q)
                                                                     for v=1:length(Q{u}(:))
                                                                         curr_idx=buf_idx{coun,:};
                                                                         qq_buf=0;
                                                                         for qq=1:10
                                                                             if ~isempty(Q{qq})
                                                                                qq_buf =qq_buf +Q{qq}(curr_idx(qq));
                                                                             end
                                                                         end
                                                                        P_M_l_z_c(u,v)=Q{u}(v)/qq_buf;
                                                                        coun=coun+1;
                                                                     end
                                                                    end
                                                                         

                                                                
                                            

                                            %% Find the corresponding index of the closest count from (xcentr) then get the corresponding value from the (P_M_l_z_c)
                                            
                                                        for msms=1:10
                                                            mass=idxidx(msms);
                                                            sdt=load(strcat('Output_Sets/Set_mass_vol',num2str(vol),'_slice_',num2str(slice),'_',num2str(mass),'.mat'));
                                                                if ~isempty(sdt.Set{idx(pos)})
                                                                    len=length(sdt.Set{idx(pos)});
                                                                    xcc=length(xcentr);
                                                                    if xcc>1
                                                                            temp_buf=abs(repmat(len,length(xcentr{msms}(:)),1)-xcentr{msms}(:));
                                                                            if ~isempty(temp_buf)
                                                                                ind_cand=find(min(temp_buf)==temp_buf);
                                                                                p{msms,pos}=P_M_l_z_c(msms,ind_cand(1));
                                                                            else
                                                                                 p{msms,pos}=0;
                                                                            end
                                                                    else
                                                                        p{msms,pos}=0;
                                                                    end
                                                                else
                                                                    p{msms,pos}=0;
                                                                end

                                                        end
                                            
                                                        chooseone=find(cell2mat(p(:,pos))==max(cell2mat(p(:,pos))));
                                                        massspace(pos)=chooseone(1);

                                                        clear P_M_l_z_c;
                                                        idx_new=[idx_new ;idx(pos) pos];
                                            end
                                            %pos
                                end
                                p=cell2mat(p)';
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                    ms=logspace(10,15,100);
                                    idxidx=round(linspace(41,90,10)); 
                                    ms=ms(idxidx);

                                    load(strcat('Input/RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));

                                    Xdata=Data(:,2:4);
                                    [Xdata,ia,ic]=unique(Xdata,'rows','stable');
                                    Data_all=[Xdata(:,1:2)*pi/180 Xdata(:,3)];
                                    Lum_all=Data(:,8);
                                    PMlG=cell(1,1);   
                                    prval=[];
                                    expected_grp=cell(1,1);
                                    
                                              for ms_gg=1:10

                                                  sd=load(strcat('Output_Sets/Set_mass_vol',num2str(vol),'_slice_',num2str(slice),'_',num2str(idxidx(ms_gg)),'.mat'));
                                                  Net=sd.Set(idx_new(:,1));
                                                  
                                                  
                                                  
                                                        Group=Net;

                                                            G=4.3E-9;
                                                            c= 299792.458;
                                                            v=cell(1,1);
                                                            for group_in=1:length(Group)

                                                                      Group_pht=Group{group_in};
                                                                      if ~isempty(Group_pht)
                                                                                              group_rad=[Group_pht(:,1)*pi/180 Group_pht(:,2)*pi/180 Group_pht(:,3)];
                                                                                              Loa=find(ismember(Data_all,group_rad,'rows'));
                                                                                              Lum=Lum_all(Loa);

                                                                                                    Centre_ID=find(max(Lum)==Lum);

                                                                                                    radian=[];
                                                                                                    normal_vector=0;vec_alongtheline_ofsighrt=0;
                                                                                                    mu=group_rad(Centre_ID,:);
                                                                                                    cor= [group_rad(1:Centre_ID-1,:); group_rad(Centre_ID+1:end,:)];
                                                                                                    vector_from_center_to_point=repmat(mu,length(cor(:,1)),1)-cor; 
                                                                                                    Y = sqrt(sum(bsxfun(@minus, repmat(mu,length(cor(:,1)),1),cor).^2,2));
                                                                                                     len_u=sqrt(mu(1)^2+mu(2)^2+mu(3)^2);
                                                                                                    unit_vec_alongtheline_ofsighrt=mu./len_u;        
                                                                                                    for dot_vec=1:length(cor(:,1))
                                                                                                        vec_alongtheline_ofsighrt(dot_vec)=dot(vector_from_center_to_point(dot_vec,:),unit_vec_alongtheline_ofsighrt');
                                                                                                        normal_vector(dot_vec)=norm(vector_from_center_to_point(dot_vec,:)-(unit_vec_alongtheline_ofsighrt*vector_from_center_to_point(dot_vec,:)'*unit_vec_alongtheline_ofsighrt));
                                                                                                        radian(dot_vec)=sqrt(Y(dot_vec).^2-vec_alongtheline_ofsighrt(dot_vec).^2);
                                                                                                    end

                                                                                                     r=mean(radian);%sum(Weight'.*radian);
                                                                                                     [~,Da]=getHzDa(mean(Group_pht(:,3)));
                                                                                                     r=(r*Da);
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                                                        curr_grp=Group_pht;
                                    %                                                                     curr_grp=[curr_grp(:,2:3) curr_grp(:,4)];
                                                                                                        curr_grp=curr_grp';
                                                                                                        mm_grp=mean(curr_grp');
                                                                                                        [~,Da]=getHzDa(mm_grp(3));
                                                                                                        XYZd=[mm_grp(1) mm_grp(2) mm_grp(3)];
                                                                                                        XYZd=XYZd';
                                                                                                        len_i=sqrt(XYZd(1)^2+XYZd(2)^2+XYZd(3)^2);
                                                                                                        u_i=XYZd./len_i;
                                                                                                        buffer1=curr_grp';
                                                                                                        GP=[buffer1(:,1) buffer1(:,2) buffer1(:,3)];
                                                                                                        XYZg=GP';
                                                                                                        len_g=sqrt(GP(:,1).^2+GP(:,2).^2+GP(:,3).^2);
                                                                                                        u_g=XYZg(1:3,:);
                                                                                                        u_g(1,:)=u_g(1,:)./len_g(:)';
                                                                                                        u_g(2,:)=u_g(2,:)./len_g(:)';
                                                                                                        u_g(3,:)=u_g(3,:)./len_g(:)';
                                                                                                        phi=acos(dot(u_g,repmat(u_i,1,length(u_g(1,:)))));
                                                                                                        R1=phi*Da;
                                                                                                        R1=mean(R1);

                                                                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                % find velocity dispercion for TG

                                                                                                 itr=group_in;
                                                                                                 v{itr}=Group_pht(:,3)*c;
                                    % %                                                               if length(v{itr}) >=4

                                                                                                          final_vel=[];
                                                                                                          ss=sort(v{itr});
                                                                                                        for j=1 : length(v{itr})-1
                                                                                                               ggg(j)=ss(j+1)-ss(j);
                                                                                                               w(j)=j*(length(v{itr})-j);
                                                                                                        end
                                                                                                        sigma_gap(itr)=(sqrt(pi)/(length(v{itr})*(length(v{itr})-1)))*sum(ggg.*w);
                                                                                                        sig_fin_f=sqrt((length(v{itr})/(length(v{itr})-1))*sigma_gap(itr)^2);
                                                                                                        sig_fin_f=sigma_gap(itr);
                                                                                                        ggg=[];w=[];
                                                                                                        FF=3;
                                                                                                        MASS_EST(ms_gg,itr)=FF*(sig_fin_f^2)*r/G; %% or directly sigma_gap instead of sig_fin_f

                                                                                                        compare=(log10(ms')-log10(repmat(MASS_EST(ms_gg,itr),length(ms),1))).^2;
                                                                                                        idxddx=find(min(compare)==compare);
                                                                                                        prval=[prval ;ms_gg group_in idxddx MASS_EST(ms_gg,itr)];
                                                                      else
                                                                          PMlG{ms_gg,group_in}=0;
                                                                          MASS_EST(ms_gg,group_in)=0;
                                                                          prval=[prval ;ms_gg group_in 0 MASS_EST(ms_gg,group_in)];
                                                                      end            

                                                            end
                                                            expected_grp{ms_gg}=prval;
                                                            clear prval;
                                                            clearvars -except p  idx_new idxidx PMlG Data_all Lum_all ms_gg vol slice ms New_idx Xdata expected_grp MASS_EST
                                                            prval=[];

                                              end
          
                                
                                
                                
%                                             for gpr=1:length(MASS_EST(1,:))
%                                                 
%                                                 grp_final(gpr)=sum(p(gpr,:)'.*MASS_EST(:,gpr));
%                                                 
%                                             end
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
end
