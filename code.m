%%%% Computational Assignment by Amit Kumar 200109 %%%%

p=563.00; % Isobaric !!! p -> constant pressure (in mmhg)
T=zeros(1,101); % T -> tempreture data
y1 =zeros(1,101); % y1 -> mole fraction of water in gase phase
x1 =zeros(1,101); % x1 -> mole fraction of water in liquid phase
j=1; % j -> an index

for i=0:0.01:1 % i -> index for value of x1
    T1 = 1730.630/(8.07131 - log10(p)) - 233.426; % by antoine's equation
    T2 = 1447.700/(7.16030 - log10(p)) - 210.000; % by antoine's equation
    Ti = T1*i + T2*(1-i);   % avarage temperature for initial guess
    g1 = exp((0.1776-0.7932*i)*(1-i)^2); % g1 -> gamma(1)
    g2 = exp((-0.219+0.7932*(1-i))*i^2); % g2 -> gamma(2)
    
    e=100; % e -> error
    
    while e>0.05
        p1 = 10^(8.07131 - 1730.630/(Ti+233.426)); % p1 -> saturated pressure for water
        p2 = 10^(7.16030 - 1447.700/(Ti+210.000)); % p2 -> saturated pressure for morpholine
        a  = p1/p2; % alfa
        p1 = p/(i*g1 + (1-i)*g2/a);
        Temp  = 1730.630/(8.07131 - log10(p1)) - 233.426;
        e  = abs(Temp-Ti)/Ti;
        Ti = Temp;
    end
    x1(j) = i;
    y1(j) = i*g1*p1/p;
    T(j) = Ti;
    j=j+1;
end

 % Tabulation of theoretical results
 disp(x1)
 disp(y1);
 disp(T);
 
 %equation for straigt line for better understandig in graph
 t=-26.4608*x1+118.2922;
 
 figure(1) % Isobaric T vs x, y diagram
 hold on
 plot(x1,T,y1,T);
 title('Isobaric T vs x, y diagram');
 xlabel('x1 and y1 --->');
 ylabel('T(degree celcius) --->');
 plot(x1,t,'--');
 hold off
 
 figure(2) % y vs x diagram
 plot(x1,y1);
 title('y vs x diagram');
 xlabel('x --->');
 ylabel('y --->');
 
% Experimental results
T_exp = [113.67, 109.10, 107.35, 105.23, 102.78, 101.67, 100.19, 97.41, 96.67, 95.27, 93.86, 93.75, 92.92, 92.69, 92.42];
x1_exp = [0.0811 0.1834 0.2242 0.2883 0.3736 0.4287 0.5135 0.6899 0.7173 0.8084 0.8792 0.8983 0.9385 0.9521 0.9816 ];
y1_exp = [0.1758 0.3427 0.4458 0.5115 0.5928 0.6370 0.6949 0.8211 0.8398 0.9042 0.9457 0.9554 0.9747 0.9803 0.9923 ];
 
 figure(3) % Comparision with experimental results
 hold on
 plot(x1,T,y1,T);
 plot(x1_exp,T_exp,y1_exp,T_exp);
 title('Comparision with experimental results');
 xlabel('x1 and y1 --->');
 ylabel('T(degree celcius) --->');
 plot(x1,t,'--');
 hold off
 
