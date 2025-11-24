function [ idx_flag, expected_radius, r500 ] = isittrue_test(mu,index,sig_z,mass,q)

            cosmopar;
            ind=1;
            Flag=0;
            [H_z,Da]=getHzDa(mu(3));
            [Wid1, r500]=NFW_cylinder_find_quantile_old(q,mass,H_z,Da,mu(3));
            expected_radius=Wid1;
                    index1=[index(:,1)*(pi/180) ,index(:,2)*(pi/180), index(:,3)];
                    mu1=[mu(1)*(pi/180) mu(2)*(pi/180) mu(3)]; 
                    cor=index1;
                    vector_from_center_to_point=repmat(mu1,length(cor(:,1)),1)-cor; 
                    Y = sqrt(sum(bsxfun(@minus, repmat(mu1,length(cor(:,1)),1),cor).^2,2));
                    len_u=sqrt(mu1(1)^2+mu1(2)^2+mu1(3)^2);
                    unit_vec_alongtheline_ofsighrt=mu1./len_u;        
                    for dot_vec=1:length(cor(:,1))
                        vec_alongtheline_ofsighrt(dot_vec)=dot(vector_from_center_to_point(dot_vec,:),unit_vec_alongtheline_ofsighrt');
                        normal_vector(dot_vec)=norm(vector_from_center_to_point(dot_vec,:)-(unit_vec_alongtheline_ofsighrt*vector_from_center_to_point(dot_vec,:)'*unit_vec_alongtheline_ofsighrt));
                        radian(dot_vec)=sqrt(Y(dot_vec).^2-vec_alongtheline_ofsighrt(dot_vec).^2);
                    end
                    f=normal_vector<=expected_radius;
                    s=vec_alongtheline_ofsighrt>=(-2*sig_z);
                    t=vec_alongtheline_ofsighrt<=(+2*sig_z);
                    Flag=find(f==1 & s==1 & t==1);
            if (Flag~=0)
               idx_flag=Flag;
            else
                idx_flag=[];
            end
end