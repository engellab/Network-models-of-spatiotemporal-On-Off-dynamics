
%This is simulation code of reduced binary-units network model

% The size of network model is given by a n*n 2-d square lattice. Each node represents a binary unit.
n=256;


%%%%%%%%%% Part I: Model Simulation %%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%


% Transition rate of each unit:       r_{off to on}= alpha_1 + beta_1 ( sum of nearby units )
%                                     r_{on to off}= alpha_2 - beta_2 ( sum of nearby units ) 
%
%
%                                     r_{off to on} (control condition)= stim_alpha_1 + stim_beta_1 ( sum of nearby units )
%                                     r_{on to off} (control condition)= stim_alpha_2 - stim_beta_2 ( sum of nearby units )
%                                     r_{off to on} (attention condition)= attn_alpha_1 + attn_beta_1 ( sum of nearby units )
%                                     r_{on to off} (attention condition)= attn_alpha_2 - atten_beta_2 ( sum of nearby units )
%
% units under attention condition: (50:100,50:100)
% units under control   condition: (n-100:n-50,n-100:n-50)

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



k=0;






%load inital binary-unit configuration for simulation. initial.spin= n*n  matrix, where each element represent the initial binary-state(0 or 1) of one unit.
% spins = n*n matrix of binary unit variable in the model. The element has discrete value:0,1.

initial = load('SPIN_init.mat');
spins = initial.spin;


% define summation of states of nearby units
nearby_spin=zeros(n,n);






% Select samples in a 11*11 square for attention condition 
ind1 = 69:1:79;
% Select samples in a 11*11 square for control condition 
ind2 = 175:1:185;
n_sample = length(ind1);

% number of trials in each simulation
num_trial = 100;
% number of iteration in each trial
num_iter = 1000;
num_time=1000;


%Tracing the time sequence of 11*11 units in attention/out conditions

rate_att = zeros(num_trial,n_sample,n_sample,num_time);
rate_out = zeros(num_trial,n_sample,n_sample,num_time);




x = 2:n-1;
y = 2:n-1;


% size of time-step 
dt=0.001;




% Update activity of network in each trial:

for iTrial = 1:num_trial

% In a given trial= iTrial, update activity of network in each iteration i:       
    for i=1:num_iter
        
        
        %Summation of states of nearby units
        
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

    
        %transition-rate matrix (n*n) for each unit
        %Update transition rate
        R_spins=alpha_1+beta_1.*(nearby_spin-k*gamma)+spins.*((alpha_2-beta_2.*(nearby_spin-k*gamma))-(alpha_1+beta_1.*(nearby_spin-k*gamma)));
        
        

        %Traceing the sequence before the a new transition
    rate_att(iTrial,:,:,i) = spins(ind1, ind1);% stimulus, attened in
    rate_out(iTrial,:,:,i) = spins(ind2, ind2);% stimulus, attened out
       
        
        
        %Decide if unit on each site will make a transiton based on transition matrix 
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
    
        

%%%%%%%%%%%%%%%%%%Spike count generation from given On/Off population phases <S> %%%%%%%%%%%%%%%%%%



% number of bins for a given Time-window
num_bin = 200;


deltaT =num_bin*dt;
%deltaT =0.2; Time-window for computing noise correlation

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





% constant On/Off firing rates: r_{on}=125 Hz, r_{off}=25 Hz 

t_att=t_att*100+25;
t_out=t_out*100+25;



% Poisson rate = (r_{on} <S> + r_{off})*deltaT
t_att = t_att*deltaT;
t_out = t_out*deltaT;


% spike counts from Poisson distribution 
p_att = poissrnd(t_att);
p_out = poissrnd(t_out);




%%%%%%%%%% Part II: Computation of noise correlation  %%%%%%%%%%%%%


%%%%%%%%%%  Computation of noise correlation across cortical column %%%%%%%%%%%%%


%computate all pairs of correlation among 11*11 units in attention condition and control condition 



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



%compuate average noise for each distance in attention and control condition 
% Define variables: MeanCorrDist_att(:,1)=spatial distance, MeanCorrDist_att(:,2)=average noise correlation in attention condition at specific distance 
%                   MeanCorrDist_out(:,1)=spatial distance, MeanCorrDist_out(:,2)=average noise correlation in control condition at specific distance 


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

%%%%%%%%%% Computation of noise correlation within a single column %%%%%%%%%%%%%

%compuate average noise for each distance in attention and control condition 
% Define variables: SinDist_att(:,1)= all pairs of noise correlation within a single cortical column in attention condition 
%                   SinDist_out(:,1)= all pairs of noise correlation within a single cortical column in control condition 





single_att = squeeze( mean(rate_att(:,:,:,end-num_bin:end),4) );
single_out = squeeze( mean(rate_out(:,:,:,end-num_bin:end),4) );



single_att=single_att*100+25;
single_out=single_out*100+25;

single_att = single_att*deltaT;
single_out = single_out*deltaT;


%%%%%%% number of units within single cortical column
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



% mean of mean(SinDist_att) is avarage noise correlation at distance=0 in atttention condition
% mean of mean(SinDist_out) is avarage noise correlation at distance=0 in control condition

MeanCorrDist_att(1,2)= mean(SinDist_att);
MeanCorrDist_out(1,2)= mean(SinDist_out);




%%%%%%%%%%% plot noise correlation as a function of distance %%%%%%%%%%%%%%%



plot(MeanCorrDist_att(1:30,1),MeanCorrDist_att(1:30,2),MeanCorrDist_out(1:30,1),MeanCorrDist_out(1:30,2))
title('Distance dependence of Spike-count correlation (Markov-Simulation)')
xlabel('distance=|i-j|')
ylabel('Spike-count correlation')
legend('Within attention','Out of attention')
        
        







