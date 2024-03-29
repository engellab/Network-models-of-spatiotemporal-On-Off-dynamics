%
%This is simulation code of dynamical system model

% The size of network model is given by a n*n 2-d square lattice. Each node represents a unit, including activity variable and adaptation variable.
n=256;


%%%%%%%%%% Part I: Model Simulation %%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%
% time-scale of adaptation variable  
epsilon = 0.05;
% Strength of interaction between nearby units 
 W = 5*0.015;
% Strength of gaussian noise 
D =8*0.0025; 
% constant term: a in the dynamical equation 
a = 0.15*ones(n);
% constant term: gamma in the dynamical equation 
gamm = 0.55;

% load initial condition of activity configuration of simulation. init.cells =  n*n  matrix, where each element represent the initial activity of one unit. init.inh =
% n*n  matrix, where each element represent the initial adaptation variable activity of one unit.
init = load('FHN_init.mat');
cells = init.cells;
inh = init.inh;


% define variables used for each iteration
sum = zeros(n);
sum2=zeros(n);
act_new=zeros(n,n);


x = 2:n-1;
y = 2:n-1;



% Define size of time-step in simulation
dt = 0.01*0.5;


 
% Select samples in a 11*11 square for attention condition 
ind1 = 69:1:79;
% Select samples in a 11*11 square for control condition 
ind2 = 175:1:185;

% size of sample
n_sample = length(ind1);



% Strength of stimulus
stim=0.15*ones(n,n);


% Strength of stimulus attentional input
att = zeros(n,n);
att(50:100,50:100) = 0.04;
s = stim + att;

% number of trials in one-simulation
num_trial = 50;

% number of iterations in each trial
num_iter = 3000;



% Define time-series variable for sample unit in attetion condition 
rate_att = zeros(num_trial,n_sample,n_sample,num_iter);
% Define time-series variable for sample unit in control condition 
rate_out = zeros(num_trial,n_sample,n_sample,num_iter);


trace_att = zeros(n_sample,n_sample,num_iter);
trace_out = zeros(n_sample,n_sample,num_iter);



% Update activity of network in each trial:

for iTrial = 1:num_trial


% In a given trial= iTrial, update activity of network in each iteration:    
    for i=1:num_iter
        % interaction terms
        sum(x,y) = (cells(x,y-1) - cells(x,y)) + (cells(x,y+1) - cells(x,y)) + ...
            (cells(x-1, y) - cells(x, y)) + (cells(x+1,y) - cells(x,y));

        sum(1,y) = (cells(1,y-1) - cells(1,y)) + (cells(1,y+1) - cells(1,y)) + ...
            (cells(n, y) - cells(1, y)) + (cells(2,y) - cells(1,y));

        sum(n,y) = (cells(n,y-1) - cells(n,y)) + (cells(n,y+1) - cells(n,y)) + ...
            (cells(n-1, y) - cells(n, y)) + (cells(1,y) - cells(n,y));

        sum(x,1) = (cells(x,n) - cells(x,1)) + (cells(x,2) - cells(x,1)) + ...
            (cells(x-1, 1) - cells(x, 1)) + (cells(x+1,1) - cells(x,1));

        sum(x,n) = (cells(x,n-1) - cells(x,n)) + (cells(x,1) - cells(x,n)) + ...
            (cells(x-1, n) - cells(x, n)) + (cells(x+1,n) - cells(x,n));

        sum2 = W*sum;

     % update activity variable 
     %  act_new [cells(t+1)] =  cells  + [ F( cells ) - inh + W*( interaction ) +  s ]*(dt/epsilon);
     %  F(.) piecewise linear function   
    z1=find(cells>0.5);
    act_new(z1)=cells(z1)+(1-cells(z1)-inh(z1)+sum2(z1)+s(z1))*(dt/epsilon);
    z2=find(cells>-0.5 & cells<=0.5);
    act_new(z2)=cells(z2)+(cells(z2)-inh(z2)+sum2(z2)+s(z2))*(dt/epsilon);
    z3=find(cells<=-0.5);
    act_new(z3)=cells(z3)+(-1-cells(z3)-inh(z3)+sum2(z3)+s(z3))*(dt/epsilon);
  
  
  
        % update adaptation variable 
        inh2 =  inh + (gamm*cells - inh + a )*dt + sqrt(2*D*dt)*randn(n);
        cells = act_new;
        inh = inh2;
        
        % store activity time-series 
        trace_att(:,:,i) = cells(ind1, ind1); % stimulus, attened in
        trace_out(:,:,i) = cells(ind2, ind2); % stimulus, attened out
        
    end
    
    rate_att(iTrial,:,:,:) = trace_att(:,:,:);
    rate_out(iTrial,:,:,:) = trace_out(:,:,:);
    
    
end


%plot(1:size(trace_out,3),squeeze(trace_out(2,2,:)))




%%%%%%%% number of time-bins associated with time-window=200 ms, scaled by 100ms On/Off durations.
num_bin = 1300;

%%%%%%%% time-window
deltaT =0.2;


att_1=rate_att(:,:,:,end-num_bin:end);
out_1=rate_out(:,:,:,end-num_bin:end);

normal_att = heaviside(att_1);
normal_out = heaviside(out_1);



%%%%%%%%%% mean On-Off population phases in given time-window  <S>
t_att = squeeze( mean(normal_att,4) );
t_out = squeeze( mean(normal_out,4) );




%%%%%%%%%% constant On/Off firing rates: r_{on}=130 Hz, r_{off}=50 Hz 

rtrial_delta_r_att=50*ones(num_trial,11,11);
rtrial_delta_r_out=50*ones(num_trial,11,11);
rtrial_r_off_att=80*ones(num_trial,11,11);
rtrial_r_off_out=80*ones(num_trial,11,11);




%%%%%%%%%% Poisson firing rate  (r_{on} <S> + r_{off})*deltaT

t_att = t_att.*rtrial_delta_r_att + rtrial_r_off_att;
t_out = t_out.*rtrial_delta_r_out + rtrial_r_off_out;


t_att = t_att*deltaT;
t_out = t_out*deltaT;


%%%%%%%%%% Generate Poisson spikes
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





%%%%%%%%%% Computation of noise correlation within a single column %%%%%%%%%%%%%

%compuate average noise for each distance in attention and control condition 
% Define variables: SinDist_att(:,1)= all pairs of noise correlation within a single cortical column in attention condition 
%                   SinDist_out(:,1)= all pairs of noise correlation within a single cortical column in control condition 


single_att_1=rate_att(:,:,:,end-num_bin:end);
single_out_1=rate_out(:,:,:,end-num_bin:end);



normal_single_att = heaviside(single_att_1);
normal_single_out = heaviside(single_out_1);



single_att = squeeze( mean(normal_single_att,4) );
single_out = squeeze( mean(normal_single_out,4) );




% Poisson rate
single_att=single_att*50+80;
single_out=single_out*50+80;

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
figure(1)

plot(MeanCorrDist_out(1:30,1),MeanCorrDist_out(1:30,2),MeanCorrDist_out(1:30,1),MeanCorrDist_att(1:30,2))
legend('control','attention')





