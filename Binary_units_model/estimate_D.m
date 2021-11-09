fun = @(x,D) erfc(x).*exp(x.^2);


q = integral(@(x) fun(x,5),-2,2);




D=0.0025/1.55;

x=0.0001:0.0001:0.001;


q1=(-0.5+x)/sqrt(2*D);
q2=( 0.5+x)/sqrt(2*D);

y=zeros(1,length(x));

for i=1:length(x)
     
    y(i) = integral(@(x) fun(x,5),q1(i),q2(i));
   
    
end



Eqn='a*exp(-b*x)+c';

startPoint=[100,1000,1000];

f1 = fit(x',y','exp1') 