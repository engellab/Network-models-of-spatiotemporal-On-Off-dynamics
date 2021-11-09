
n=256;

%stimulus and attention condition

stim_alpha_1 = zeros(n,n);
stim_alpha_2 = zeros(n,n);
stim_beta_1 = zeros(n,n);
stim_beta_2 = zeros(n,n);

stim_alpha_1(50:100,50:100) = 1; 
stim_alpha_2(50:100,50:100) = -2;
stim_beta_1(50:100,50:100) = -1;
stim_beta_2(50:100,50:100) = -1;


stim_alpha_1(n-100:n-50,n-100:n-50) = 1; 
stim_alpha_2(n-100:n-50,n-100:n-50) = -2;
stim_beta_1(n-100:n-50,n-100:n-50) = -1;
stim_beta_2(n-100:n-50,n-100:n-50) = -1;



attn_alpha_1 = zeros(n,n);
attn_alpha_2 = zeros(n,n);
attn_beta_1 = zeros(n,n);
attn_beta_2 = zeros(n,n);

%attn_alpha_1(50:100,50:100) = 1; 
%attn_alpha_2(50:100,50:100) = -2;
attn_beta_1(50:100,50:100) = -2;
attn_beta_2(50:100,50:100) = -2;


attn_alpha_1(50:100,50:100) = 0; 
attn_alpha_2(50:100,50:100) = 0;

%parameters

alpha_1=2*(5*ones(n)+stim_alpha_1+attn_alpha_1);
alpha_2=2*(5*ones(n)+stim_alpha_2+attn_alpha_2);
beta_1=2*(5*ones(n)+stim_beta_1+attn_beta_1);
beta_2=2*(5*ones(n)+stim_beta_2+attn_beta_2);
gamma=0.2*ones(n);
%k=2;
k=0;


%load inital spin configuration

initial = load('SPIN_init.mat');
spins = initial.spin;


nearby_spin=zeros(n,n);






ind1 = 69:1:79;
ind2 = 175:1:185;
n_sample = length(ind1);


num_trial = 100;
num_iter = 1000;
num_time=1000;


%Tracing the time sequence of 11*11 units in attention/out conditions

rate_att = zeros(num_trial,n_sample,n_sample,num_time);
rate_out = zeros(num_trial,n_sample,n_sample,num_time);




x = 2:n-1;
y = 2:n-1;



dt=0.001;




for iTrial = 1:num_trial
    
    for i=1:num_iter
        
        
        %Summation of nearby units
        
         nearby_spin(x,y) = (spins(x,y-1) - spins(x,y)) + (spins(x,y+1) - spins(x,y)) + ...
            (spins(x-1, y) - spins(x, y)) + (spins(x+1,y) - spins(x,y));

        nearby_spin(1,y) = (spins(1,y-1) - spins(1,y)) + (spins(1,y+1) - spins(1,y)) + ...
            (spins(n, y) - spins(1, y)) + (spins(2,y) - spins(1,y));

        nearby_spin(n,y) = (spins(n,y-1) - spins(n,y)) + (spins(n,y+1) - spins(n,y)) + ...
            (spins(n-1, y) - spins(n, y)) + (spins(1,y) - spins(n,y));

        nearby_spin(x,1) = (spins(x,n) - spins(x,1)) + (spins(x,2) - spins(x,1)) + ...
            (spins(x-1, 1) - spins(x, 1)) + (spins(x+1,1) - spins(x,1));

        nearby_spin(x,n) = (spins(x,n-1) - spins(x,n)) + (spins(x,1) - spins(x,n)) + ...
            (spins(x-1, n) - spins(x, n)) + (spins(x+1,n) - spins(x,n));

    
        %transition-rate matrix (n*n) for each spin
        %Update transition rate
        R_spins=alpha_1+beta_1.*(nearby_spin-k*gamma)+spins.*((alpha_2-beta_2.*(nearby_spin-k*gamma))-(alpha_1+beta_1.*(nearby_spin-k*gamma)));
        
        

        %Traceing the sequence before the a new transition
    rate_att(iTrial,:,:,i) = spins(ind1, ind1);% stimulus, attened in
    rate_out(iTrial,:,:,i) = spins(ind2, ind2);% stimulus, attened out
       
        
        
        %Decide if unit on each site will make a transiton
        for aind = 1:n
            for bind= 1:n
                probind=rand;
                  if probind >= exp(-R_spins(aind,bind)*dt)
                      spins(aind,bind)=1-spins(aind,bind);
                  end
            end
        end
        
    end
    
end
    
        

%Spike count generation

delta_r_att = unifrnd(5,200,[11,11]);
delta_r_out = unifrnd(5,200,[11,11]);

r_off_att = unifrnd(5,50,[11,11]);
r_off_out = unifrnd(5,50,[11,11]);

num_bin = 200;
%num_bin=50;

deltaT =num_bin*dt;
%deltaT =0.2;

t_att = squeeze( mean(rate_att(:,:,:,end-num_bin:end),4) );
t_out = squeeze( mean(rate_out(:,:,:,end-num_bin:end),4) );

rtrial_delta_r_att=zeros(num_trial,11,11);
rtrial_delta_r_out=zeros(num_trial,11,11);
rtrial_r_off_att=zeros(num_trial,11,11);
rtrial_r_off_out=zeros(num_trial,11,11);

for rtrial_delta_r_att_num=1:num_trial
    rtrial_delta_r_att(rtrial_delta_r_att_num,:,:)=delta_r_att;
end

for rtrial_delta_r_out_num=1:num_trial
    rtrial_delta_r_out(rtrial_delta_r_out_num,:,:)=delta_r_out;
end

for rtrial_r_off_att_num=1:num_trial
    rtrial_r_off_att(rtrial_r_off_att_num,:,:)=r_off_att;
end

for rtrial_r_off_out_num=1:num_trial
    rtrial_r_off_out(rtrial_r_off_out_num,:,:)=r_off_out;
end



t_att=t_att*100+25;
t_out=t_out*100+25;

%t_att = t_att.*rtrial_delta_r_att + rtrial_r_off_att;
%t_out = t_out.*rtrial_delta_r_out + rtrial_r_off_out;


t_att = t_att*deltaT;
t_out = t_out*deltaT;


p_att = poissrnd(t_att);
p_out = poissrnd(t_out);

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

%histogram(SortCorrDist_att(:,2));
%mean(SortCorrDist_att(:,2));

%histogram(SortCorrDist_out(:,2));
%mean(SortCorrDist_out(:,2));


%Sort spike-count correlation as function of distance

MeanCorrDist_att=zeros(1+length(unique(SortCorrDist_att(:,1))),2);
MeanCorrDist_out=zeros(1+length(unique(SortCorrDist_out(:,1))),2);


cumattn=SortCorrDist_att(1,2);
cumout=SortCorrDist_out(1,2);
subo=1;
suboo=2;
subt=1;
subtt=2;


%Average all of piarwise correlations with the same distance 


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correlation within single cortical column




single_att = squeeze( mean(rate_att(:,:,:,end-num_bin:end),4) );
single_out = squeeze( mean(rate_out(:,:,:,end-num_bin:end),4) );



single_att=single_att*100+25;
single_out=single_out*100+25;

single_att = single_att*deltaT;
single_out = single_out*deltaT;


num_single=10;

psingle_att=zeros(num_trial,n_sample,n_sample,num_single);
psingle_out=zeros(num_trial,n_sample,n_sample,num_single);

SinDist=1;
SinDist_att=zeros(1,(num_single*num_single-num_single)/2*n_sample*n_sample);
SinDist_out=zeros(1,(num_single*num_single-num_single)/2*n_sample*n_sample);


for sk=1:num_single
    
psingle_att(:,:,:,sk) = poissrnd(single_att);
psingle_out(:,:,:,sk) = poissrnd(single_out);

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

MeanCorrDist_att(1,2)= mean(SinDist_att);
MeanCorrDist_out(1,2)= mean(SinDist_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(MeanCorrDist_att(:,1),MeanCorrDist_att(:,2))
%plot(MeanCorrDist_out(:,1),MeanCorrDist_out(:,2))
%plot(MeanCorrDist_out(1:50,1),MeanCorrDist_out(1:50,2)-MeanCorrDist_att(1:50,2))

%f = fit(MeanCorrDist_out(1:50,1),MeanCorrDist_out(1:50,2)-MeanCorrDist_att(1:50,2),'exp1');
%plot(f,MeanCorrDist_out(1:50,1),MeanCorrDist_out(1:50,2)-MeanCorrDist_att(1:50,2));

plot(MeanCorrDist_att(1:30,1),MeanCorrDist_att(1:30,2),MeanCorrDist_out(1:30,1),MeanCorrDist_out(1:30,2))
title('Distance dependence of Spike-count correlation (Markov-Simulation)')
xlabel('distance=|i-j|')
ylabel('Spike-count correlation')
legend('Within attention','Out of attention')
        
        






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sequence
%ax1 = subplot(2,1,1); % top subplot
%plot(ax1,1:1000,squeeze(rate_att(1,2,3,1:1000)),'color','red','Linewidth',1.2)
%title(ax1,'On-off transitin sequence (Simulation) (Attention)(alpha1=12, alpha2=6, beta1=4)')
%ylabel(ax1,'Binary States')
%xlabel(ax1,'Time (ms)')
%ax2 = subplot(2,1,2); % bottom subplot
%plot(ax2,1:1000,squeeze(rate_out(1,2,3,1:1000)),'color','blue','Linewidth',1.2)
%ylabel(ax2,'Binary States')
%xlabel(ax2,'Time (ms)')
%title(ax2,'On-off transitin sequence (Simulation) (Out)(alpha1=12, alpha2=6, beta1=8)')
%axis(ax1,[0 1000 -0.5 1.5])
%axis(ax2,[0 1000 -0.5 1.5])
%legend(ax1,'Attention, tau-off ~1/alpha1=83ms, tau-on ~1/alpha2=167ms ')
%legend(ax2,'Out, tau-off ~1/alpha1=83ms, tau-on ~1/alpha2=167ms ')
