function [ New_p,idx ] = check_gals( X_orig,p,mass,idx_vals)


cosmopar;
ind=1;
kk=1;
flg=0;

           X_orig=[X_orig(:,1)*(pi/180)  X_orig(:,2)*(pi/180) X_orig(:,3)];
           p=[p(:,1)*(pi/180)  p(:,2)*(pi/180) p(:,3)];
           p_orig=p;
           
q=0.99999;
New_p=[];
idx=[];
    while ~isempty(p)
            mu=p(kk,:);
            
                                sss=0.26*log(10)/3;
                                sig_v2=((1.55e8/10^0.0403)*(mass*E(mu(3))/10^14)^0.94)^(1/3)*exp(sss/2);
                                sig_z=(1+mu(3))*(sig_v2/c);
                                Max_width=2*sig_z;
  
            [H_z,Da]=getHzDa(mu(3));

            [Wid1]=NFW_cylinder_find_quantile_old(q,mass,H_z,Da,mu(3));

            expected_radius=Wid1;
            
            
                    cor=X_orig;%p(kk:end,:);
                    vector_from_center_to_point=repmat(mu,length(cor(:,1)),1)-cor;
                    len_u=sqrt(mu(1)^2+mu(2)^2+mu(3)^2);
                    unit_vec_alongtheline_ofsighrt=mu./len_u;        
                    for dot_vec=1:length(cor(:,1))
                        vec_alongtheline_ofsighrt(dot_vec)=dot(vector_from_center_to_point(dot_vec,:),unit_vec_alongtheline_ofsighrt');
                        normal_vector(dot_vec)=norm(vector_from_center_to_point(dot_vec,:)-(unit_vec_alongtheline_ofsighrt*vector_from_center_to_point(dot_vec,:)'*unit_vec_alongtheline_ofsighrt));
                    end
                    f=normal_vector<=expected_radius;
                    s=vec_alongtheline_ofsighrt>=(-Max_width);
                    t=vec_alongtheline_ofsighrt<=(+Max_width);
                    idx_univ=find(f==1 & s==1 & t==1);
                    if ~isempty(idx_univ)

                        
                        if length(idx_univ)<=3
                           
                                    [~,idx2]=setdiff(p(:,1:3),p(1,1:3),'rows','stable');
                                    p=p(idx2,:);
                                    idx_vals=idx_vals(idx2,:);
                        else
                                    
                                            New_p=[New_p; p(1,1:3)];
                                            idx=[idx ; idx_vals(1)];

                            [~,idx2]=setdiff(p(:,1:3),p(1,1:3),'rows','stable');
                            p=p(idx2,:);
                            idx_vals=idx_vals(idx2,:);
                        end

                    else
                        [~,idx2]=setdiff(p(:,1:3),p(1,1:3),'rows','stable');
                        p=p(idx2,:);
                        idx_vals=idx_vals(idx2,:);
                    end

                        clear vec_alongtheline_ofsighrt normal_vector f s t;
    end
  New_p=[New_p(:,1)*(180/pi)  New_p(:,2)*(180/pi) New_p(:,3)];
end


