%this code is the last version based on the local density mass

for vol=1:9
    for slice=1:3

         load(strcat('RA_Dec_Z_Vol_',num2str(vol),'_Slice',num2str(slice),'.mat'));
         Xdata=Data(:,2:4); 
         Xdata=unique(Xdata,'rows','stable');

        %The positions
        Pos=Xdata;

        indexx=round(linspace(1,31,10));
        idxidx=round(linspace(41,90,10));
        %loop for masses
     
        % contains all the expeted masses of the locations
                load(strcat('Expected_mass_no_5_',num2str(vol),'_',num2str(slice),'.mat'));

% contains the prec. vs .recall vals via changing the threshold value 
                load(strcat('PHTALK_paper_corr_LDm_',num2str(vol),'_Slice_',num2str(slice),'_',num2str(5),'.mat'));
                
% just to select the interested location via their own indices
                 load(strcat('Prob_paper_',num2str(vol),'_',num2str(slice),'.mat'));
% each predected location with its corresponding PHT value
                load(strcat('PHT_paper_',num2str(vol),'_Slice_',num2str(slice),'.mat'));
                
                idx=1;

%% Prior 1  P

 
 %for each set we will investigate how many neighbours to each memeber's
 %set and find the maximum no. to use it later
 
         %Hough_vals=cell(1,1);
         New_set=cell(1,1);
         coun=1;
%          mass_no_final=mass_no_final(Expected_group_idx);
         index_zero=find(mass_no_final~=0);
         mass_no_final=mass_no_final(index_zero);
         Expected_group_idx=Expected_group_idx(index_zero);
          cosmopar;
         Group_Index=cell(1,1);
         for k=1:length(Expected_group_idx)
                    mass=mass_no_final(k);%idxidx(mass_no_final(k));
                    curr=Pos(Expected_group_idx(k),:);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     load(strcat('Set_mass_vol',num2str(vol),'_slice_',num2str(slice),'_',num2str(mass),'.mat'));
%                     buf=Set{Expected_group_idx(k)};
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    sss=0.26*log(10)/3;
                    sig_v2=((1.55e8/10^0.0403)*(mass*E(curr(3))/10^14)^0.94)^(1/3)*exp(sss/2);
                    sig_z=((curr(3)+1)*sig_v2/c);
                    [ idx_flag ] = isittrue_test(curr,Xdata,sig_z,mass,0.9999999);
                    buf=Xdata(idx_flag,:);
                    
                    if ~isempty(buf)
%                        
                        len_final=length(buf(:,1));
                        %find the index for each memeber in the original mat to use it
                        %to find the number of its memeber in the set matrix
                            [~,Locb] = ismember(buf,Pos,'rows');
                            county=[];
                        
                            if ~isempty(Locb) &&  length(Locb)>4
                                
                                               %assign a probability value from
                                               %the model
                                               centre=Pos(Expected_group_idx(k),:);

                                                    P_final1=Hough(k);
                                                        
                                                idx_th=find(P_final1>=(threshold(idx)));
                                                if idx_th
                                                    mass_fn(k)=mass_no_final(k);
                                                    New_set{k}=[ Pos(Locb,:)];
                                                    %Hough_vals{k}=[Hough(Locb)];
                                                   Group_Index{k}=Expected_group_idx(k);%k;%Expected_group_idx(k);
                                                end
                                
                            end
                    end
                        
         end
         %save the prior probability of the galaxies
         
         New_set=New_set(~cellfun('isempty',New_set));
         %Hough_vals=Hough_vals(~cellfun('isempty',Hough_vals));
         idid=find(~cellfun('isempty',Group_Index));
         Group_Index=Group_Index(~cellfun('isempty',Group_Index));
         mass_fn=mass_fn(idid);
         mass_no_final=mass_fn;
         
         %if there are more than half of them similar then I will combine
         %them into one group
        mass_cell_est=num2cell(mass_no_final);
         for loop1=1:length(New_set)
             buff_set1=New_set{loop1};
             buff_idx1=Group_Index{loop1};
             %buff_hough1=Hough_vals{loop1};
             
             if (~isempty(buff_set1))
                     len1=length(buff_set1(:,1));
                     for loop2=loop1+1:length(New_set)
                        buff_set2=New_set{loop2};
                        buff_idx2=Group_Index{loop2};
                        %buff_hough2=Hough_vals{loop2};
                        
                        if (~isempty(buff_set2))
                            flag=ismember(buff_set1,buff_set2,'rows');
                            len2=length(buff_set2(:,1));
                            if length(find(flag==1))>=(.5*len2)
                                buff_merge=[buff_set1;buff_set2];
                                buff_merge_idx=[buff_idx1 ; buff_idx2];
                                 %buff_merge_hough=[buff_hough1; buff_hough2 ];
                                
                                
                                [New_set{loop1},ia,~]=unique(buff_merge,'rows','stable');
                                New_set{loop2}=[];
                                
                                %Hough_vals{loop1}=buff_merge_hough(ia);
                                %Hough_vals{loop2}=[];
                                
                                mass_cell_est{loop2}=[];
                                Group_Index{loop1}=buff_merge_idx(1);
                                Group_Index{loop2}=[];
                                buff_set1=New_set{loop1};
                                buff_idx1=Group_Index{loop1};
                                len1=length(buff_set1(:,1));
                            end
                        end
                    end
             end
         end         
          New_set=New_set(~cellfun('isempty',New_set));
          %Hough_vals=Hough_vals(~cellfun('isempty',Hough_vals));
          
          Group_Index=Group_Index(~cellfun('isempty',Group_Index));
         mass_cell_est=mass_cell_est(~cellfun('isempty',mass_cell_est));
         mass_no_final=cell2mat(mass_cell_est);
            save(strcat('G_V_paper_','Vol_',num2str(vol),'_S_',num2str(slice)),'New_set','Group_Index','mass_no_final');%'Hough_vals','New_set','Group_Index','mass_no_final');

    clearvars -except slice vol
    end
end
% matlabpool close;
