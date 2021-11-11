%%%%%%%% Computation of noise correlation with On/Off rate sampled from recording data

%%%% 1. input of  simulation result: population On/Off phases: <S>
% <S>= t_att_1 (attention condition)
% <S>= t_out_1 (control condition)


att_1=rate_att(:,:,:,end-num_bin:end);
out_1=rate_out(:,:,:,end-num_bin:end);

normal_att = heaviside(att_1);
normal_out = heaviside(out_1);


t_att_1 = squeeze( mean(normal_att,4) );
t_out_1 = squeeze( mean(normal_out,4) );


%%%% 2. input of On/Off rate distriubtion from recording data:
%  roff_cont: r_{off} in control condition 
%  roff_att:  r_{off} in attention condition 
%  ron_cont   r_{on}  in control condition 
%  ron_att    r_{on}  in attention condition 



% compute noise correlation with On/Off rate sampled from recording data, with fixed On/Off population phase <S>
% To aviod insufficent sampling of On/Off rate, we reapat 10 times (num_resample=10)



%%%%%%%%%%%%%%%%% noise correlation with resampled On/Off rate  across cortical column  %%%%%%%%%%%%

% iteration number of resampling
num_resample=10;

% noise correlation with resampled On/Off rate:
% attention conidtion 
NC_att_resample=zeros(61,num_resample);
% control condition 
NC_out_resample=zeros(61,num_resample);


% iteration of resampling
for i_resample=1:10



%%%%%% generate random indices, and sample On/Off rate based on these indices 
r_g1_off_cont=rand(11);
r_g1_off_att=rand(11);
r_g1_on_cont=rand(11);
r_g1_on_att=rand(11);

r_g2_off_cont=ceil(r_g1_off_cont.*(31*8*16));
r_g2_off_att=ceil(r_g1_off_att.*(31*8*16));
r_g2_on_cont=ceil(r_g1_on_cont.*(31*8*16));
r_g2_on_att=ceil(r_g1_on_att.*(31*8*16));


r_g3_off_cont=zeros(11,11);
r_g3_off_att=zeros(11,11);
r_g3_on_cont=zeros(11,11);
r_g3_on_att=zeros(11,11);


for rgi=1:11
    for rgi2=1:11
    
        r_g3_off_cont(rgi,rgi2) = roff_cont ( r_g2_off_cont(rgi,rgi2));
        r_g3_off_att(rgi,rgi2)  = roff_att  ( r_g2_off_att(rgi,rgi2));
        r_g3_on_cont(rgi,rgi2)  = ron_cont  ( r_g2_on_cont(rgi,rgi2));
        r_g3_on_att(rgi,rgi2)   = ron_att   ( r_g2_on_att(rgi,rgi2));
        
        
        
    end
    
    
    
end


rtrial_delta_r_att=zeros(num_trial,11,11);
rtrial_delta_r_out=zeros(num_trial,11,11);
rtrial_r_off_att=zeros(num_trial,11,11);
rtrial_r_off_out=zeros(num_trial,11,11);


for nti=1:num_trial
    
rtrial_delta_r_att(nti,:,:) = r_g3_on_att -  r_g3_off_att;
rtrial_delta_r_out(nti,:,:) = r_g3_on_cont -  r_g3_off_cont;
rtrial_r_off_att(nti,:,:)   = r_g3_off_att  ;
rtrial_r_off_out(nti,:,:)   = r_g3_off_cont ;
    
    
end




% Compuate Poisson rate based on rampled On/Off rates
t_att_2 = t_att_1.*rtrial_delta_r_att + rtrial_r_off_att;
t_out_2 = t_out_1.*rtrial_delta_r_out + rtrial_r_off_out;


t_att = t_att_2*deltaT;
t_out = t_out_2*deltaT;



% Compuate Poisson spikes based on rampled On/Off rates
p_att = poissrnd(t_att);
p_out = poissrnd(t_out);




% Compuate noise correlation across cortical column 
CorrDist_att=zeros(15000,2);
CorrDist_out=zeros(15000,2);
Corrind_att=1;
Corrind_out=1;

for ox=1:11
    for oy=1:11
        for tx=1:11
            for ty=1:11
                if  (((ox-tx)^2+(oy-ty)^2)>0)
       CorrDist_att(Corrind_att,1)=sqrt((ox-tx)^2+(oy-ty)^2);
       CorrDist_att(Corrind_att,2)=corr(p_att(:,ox,oy),p_att(:,tx,ty));
       Corrind_att=Corrind_att+1;
                end
            end
        end
    end
end
               
SortCorrDist_att=sortrows(CorrDist_att(1:(Corrind_att-1),:));



for ox=1:11
    for oy=1:11
        for tx=1:11
            for ty=1:11
                if  (((ox-tx)^2+(oy-ty)^2)>0)
       CorrDist_out(Corrind_out,1)=sqrt((ox-tx)^2+(oy-ty)^2);
       CorrDist_out(Corrind_out,2)=corr(p_out(:,ox,oy),p_out(:,tx,ty));
       Corrind_out=Corrind_out+1;
                end
            end
        end
    end
end
               
SortCorrDist_out=sortrows(CorrDist_out(1:(Corrind_out-1),:));




MeanCorrDist_att=zeros(length(unique(SortCorrDist_att(:,1))),2);
MeanCorrDist_out=zeros(length(unique(SortCorrDist_out(:,1))),2);


cumattn=SortCorrDist_att(1,2);
cumout=SortCorrDist_out(1,2);
subo=1;
suboo=2;
subt=1;
subtt=2;

for ko=1:length(SortCorrDist_att(:,1))
    if (ko==length(SortCorrDist_att(:,1)))
        MeanCorrDist_att(suboo,1)=SortCorrDist_att(ko,1);
        MeanCorrDist_att(suboo,2)=cumattn/subo; 
    elseif (SortCorrDist_att(ko,1)==SortCorrDist_att(ko+1,1))
        cumattn = cumattn + SortCorrDist_att(ko+1,2);
        subo=subo+1;
    else
        MeanCorrDist_att(suboo,1)=SortCorrDist_att(ko,1);
        MeanCorrDist_att(suboo,2)=cumattn/subo;
        suboo=suboo+1;
        subo=1;
        cumattn=SortCorrDist_att(ko+1,2);
    end
end


for kt=1:length(SortCorrDist_out(:,1))
    if (kt==length(SortCorrDist_out(:,1)))
        MeanCorrDist_out(subtt,1)=SortCorrDist_out(kt,1);
        MeanCorrDist_out(subtt,2)=cumout/subt; 
    elseif (SortCorrDist_out(kt,1)==SortCorrDist_out(kt+1,1))
        cumout = cumout + SortCorrDist_out(kt+1,2);
        subt=subt+1;
    else
        MeanCorrDist_out(subtt,1)=SortCorrDist_out(kt,1);
        MeanCorrDist_out(subtt,2)=cumout/subt;
        subtt=subtt+1;
        subt=1;
        cumout=SortCorrDist_out(kt+1,2);
    end
end



NC_att_resample(:,i_resample)=MeanCorrDist_att(:,2);
NC_out_resample(:,i_resample)=MeanCorrDist_out(:,2);



end






%%%%%%%%%%%%%%%%% noise correlation with resampled On/Off rate  with a single cortical column  %%%%%%%%%%%%


%%%%%%%%% On/Off population phases: 

single_att_1=rate_att(:,:,:,end-num_bin:end);
single_out_1=rate_out(:,:,:,end-num_bin:end);



normal_single_att = heaviside(single_att_1);
normal_single_out = heaviside(single_out_1);



single_att = squeeze( mean(normal_single_att,4) );
single_out = squeeze( mean(normal_single_out,4) );






% iteration of resampling
for i_resample=1:num_resample



num_single=10;

psingle_att=zeros(num_trial,n_sample,n_sample,num_single);
psingle_out=zeros(num_trial,n_sample,n_sample,num_single);

SinDist=1;
SinDist_att=zeros(1,(num_single*num_single-num_single)/2*n_sample*n_sample);
SinDist_out=zeros(1,(num_single*num_single-num_single)/2*n_sample*n_sample);



    F_single_att=zeros(num_trial,n_sample,n_sample);
    F_single_out=zeros(num_trial,n_sample,n_sample);
    
    

for sk=1:num_single
    
 
 
 
 %%%%%% generate random indices, and sample On/Off rate based on these indices    

   sample  = rand(n_sample);
    
    r1_t_on  = ron_cont(ceil(sample *31*8*16));
    r1_t_off  = roff_cont(ceil(sample *31*8*16));
    
    
    
    r2_t_on  = ron_att(ceil(sample *31*8*16));
    r2_t_off  = roff_att(ceil(sample *31*8*16));
    
    

    
    
    dr1_t = (r2_t_on - r2_t_off);
    
   
    
   
   for i=1:n_sample
     
      for j=1:n_sample
         
          
          
        for t=1:num_trial    
          
          
         F_single_out(t,i,j)=(single_out(t,i,j)*dr1_t(i,j)+r1_t_off(i,j))*deltaT;
         F_single_att(t,i,j)=(single_att(t,i,j)*dr1_t(i,j)+r1_t_off(i,j))*deltaT; 
         
         
         
         
         
        end
      
      end
  end 
  
  
    
    
% generate Poisson spikes 
 
psingle_att(:,:,:,sk) = poissrnd(F_single_att);
psingle_out(:,:,:,sk) = poissrnd(F_single_out);


end


for sx=1:11
    for sy=1:11
        


for si=1:num_single-1
    for sj=si+1:num_single
    

SinDist_att(SinDist) = corr(psingle_att(:,sx,sy,si),psingle_att(:,sx,sy,sj));
SinDist_out(SinDist) = corr(psingle_out(:,sx,sy,si),psingle_out(:,sx,sy,sj));
SinDist=SinDist+1;

    end
end


    end
end



NC_att_resample(1,i_resample)= mean(SinDist_att);
NC_out_resample(1,i_resample)= mean(SinDist_out);



end


%%%%%%%%%%%% average noise correlations obtained from each iteration of resampling %%%%%%%%%%%%%%

MeanCorrDist_att(:,2)=mean(NC_att_resample,2);
MeanCorrDist_out(:,2)=mean(NC_out_resample,2);






%%%%%%%%%%%%% plot noise correlation as a function of distance %%%%%%%%%%%%%%%
figure(2)


plot(MeanCorrDist_out(1:30,1),MeanCorrDist_out(1:30,2),MeanCorrDist_out(1:30,1),MeanCorrDist_att(1:30,2))
legend('control','attention')




