function f = findbetaSEIR(x,theta)

% original one 
R0_target = theta(6) ;

% thetas
% theta(1) = N  
% theta(2) = alpha
% x = beta
% theta(4) = mu
% theta(5) = gamma
% based on R0 = (alpha*beta)/((alpha + mu)*(gamma + mu))

% conditions of the equation - frequency dependent transmission bSI/N and
% population size NOT considered 

%new one with new b
R0_i = (theta(2)*x)/((theta(2) + theta(4))*(theta(5) + theta(4)));
 
f = abs(R0_target-R0_i);

% theta()

