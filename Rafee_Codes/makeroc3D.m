function [tp,fp,fn,precision, recal, thresh,GT_not_det, R, R500 ] = makeroc3D(GT_mass,indexorg,val,index,mass,mass1,GTALL,GT_ALLmass) 

w=1;
TPR=0;
FPR=0;
FNR=0;
thresh=0;
flfl=0;
mm=0;
fn=0;
mm=(max(val)-min(val))/100;
thresh(1)=min(val);
coun=0;

q=0.99999;
        for th=1:100
                   tp(th)=0;
                   fn(th)=0;
                   fp(th)=0;
                   idx_over= find(val>= thresh(th));
                   [tp(th),fp(th),T,GT_not_detected, r, r500]=get_tpfp(GT_mass,indexorg,index(idx_over,:),mass,q,GTALL,GT_ALLmass);
                   GT_not_det{th}=unique(GT_not_detected','rows');
                   fn(th)=length(indexorg(:,1))-tp(th);
                   R(th) = r;
                   R500(th) = r500;

                   if tp(th)+fp(th)~=0
                        precision(th)=tp(th)/(tp(th)+fp(th));
                   else
                       precision(th)=0;
                   end
                   if (tp(th)+fn(th))~=0
                        recal(th)=tp(th)/(tp(th)+fn(th));
                   else
                       recal(th)=0;
                   end
                   thresh(th+1)=thresh(th)+mm;   

        end


end


