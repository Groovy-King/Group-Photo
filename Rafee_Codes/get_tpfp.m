function [ tp,fp,T,GT_not_detected, expected_radius, r500 ] = get_tpfp(GT_mass,indexorg,index,~,q,GTALL,GT_ALLmass)
    tp=0;
    fp=0;
    GT_not_detected=[];
    emp=0;
    tpsubforfp=0;
    gold=length(index);
    cor=index;
    cosmopar;
    T=[];F=[];
            for k=1: length(indexorg(:,1))
                        x=indexorg(k,1);
                        y=indexorg(k,2);
                        z=indexorg(k,3);
                        mu=[x,y,z];
                        %in radian
%                         rho=sqrt(mu(1).^2+mu(2).^2+mu(3).^2);
%                             GT_mass(k)=GT_mass(k)/1.5;%mass;10^14;%
%                             sig_v=((1.55e8/10^0.0403)*((GT_mass(k))*E(z)/10^14)^0.95)^(1/3)*10^(1.15192*0.0081);
%                             sig_z=(z+1)*sig_v/c;%/1.5;
%                             
                                sss=0.26*log(10)/3;
                                sig_v2=((1.55e8/10^0.0403)*(GT_mass(k)*E(z)/10^14)^0.94)^(1/3)*exp(sss/2);
                                sig_z=(1+z)*(sig_v2/c);
%                                 sig_z=(sig_v2/c);
                            
%                             Min_length_cy=-(2*sig_z);
%                             Max_length_cy=(2*sig_z);
%                             len_cy=Max_length_cy-Min_length_cy;
                            %in cartesian
                            [idx_flag, expected_radius, r500]=isittrue_test(mu,index,sig_z,(GT_mass(k)),q);
                        %save the nearest peak of the Grand truth
                                if ~isempty(idx_flag)
                                            T=[T index(idx_flag,:)'];
                                            %%% if i delete this bellow line I will consider that I detect maybe two close GT even if I've detect just one of them 
%                                             [index,~] = setdiff(index,index(idx_flag,:),'rows'); 
                                            tp=tp+1;
                                            if isempty(index)
                                                break;
                                            end
                                else
                                    GT_not_detected=[GT_not_detected; mu];
                                end
            end
            temporg=indexorg;
%             for peak1=1:length(index(:,1))
            if ~isempty(index) && ~isempty(T)
                [index,~] = setdiff(index,T','rows'); 
            else
                if isempty(index)
                    emp=1;
                end
            end
            %%%%%%%%%%%%%Just to filter out the real groups but have less
            %%%%%%%%%%%%% members till =2
            if ~isempty(GTALL) && ~isempty(index)
                TP_filter_out=[];
                indexorg=GTALL;
                for k=1: length(indexorg(:,1))
                        x=indexorg(k,1);
                        y=indexorg(k,2);
                        z=indexorg(k,3);
                        mu=[x,y,z];
%                         GT_ALLmass(k)=GT_ALLmass(k)/1.5; %mass;10^14;%
%                         rho=sqrt(mu(1).^2+mu(2).^2+mu(3).^2);
%                             sig_v=((1.55e8/10^0.0403)*((GT_ALLmass(k))*E(z)/10^14)^0.95)^(1/3)*10^(1.15192*0.0081);

                                sss=0.26*log(10)/3;
                                sig_v2=((1.55e8/10^0.0403)*(GT_ALLmass(k)*E(z)/10^14)^0.94)^(1/3)*exp(sss/2);
                                
                                sig_z=(1+z)*sig_v2/c;%/1.5;
%                             sig_z=(sig_v2/c);
%                             Min_length_cy=-(2*sig_z);
%                             Max_length_cy=(2*sig_z);
%                             len_cy=Max_length_cy-Min_length_cy;
                        idx_flag=isittrue_test(mu,index,sig_z,GT_ALLmass(k),q);
                        %save the nearest peak of the Grand truth
                                if idx_flag~=0
                                            TP_filter_out=[TP_filter_out index(idx_flag,:)'];
                                            if isempty(index)
                                                break;
                                            end
                                end
                end
                if ~isempty(TP_filter_out)
                    [index,~] = setdiff(index,TP_filter_out','rows'); 
                end
                 if isempty(index)
                    emp=1;
                end
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%
peak1=0;
if (emp==1)
    F=[];
    fp=0;
else
% % %                 while length(index(:,1))>=1
% % %                         peak1=peak1+1;
% % %                         x=index(peak1,1);
% % %                         y=index(peak1,2);
% % %                         z=index(peak1,3);
% % %                         mu=[x,y,z];
% % %                         idx_flag=isittrue(mu,index,len_cy,mass,q);
% % %                                 if idx_flag~=0
% % %                                     F=[F index(idx_flag,:)'];
% % %                                       %%% if i delete this bellow line I will consider that I detect maybe two close GT even if I've detect just one of them 
% % % %                                             [temporg,~] = setdiff(temporg(idx_flag,:),temporg,'rows'); 
% % %                                             fp=fp+1;
% % %                                             [index,~] = setdiff(index,index(idx_flag,:),'rows'); 
% % %                                             peak1=0;
% % %                                             if isempty(index)
% % %                                                 break;
% % %                                             end
% % %                                 end
% % %                 end
        fp=length(index);
            
end
end

