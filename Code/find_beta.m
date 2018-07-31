function [beta_city,beta_cityt,R0_city] = find_beta(n,popsize,alpha,mu,gamma,R0_max,season_temp)
beta_city = ones(n,1);
beta_cityt = ones(n,365);
R0_city = ones(n,1);
theta = zeros(6,1);
for i=1:n
    
    R0_city(i) = R0_max*(mean(season_temp(i,:))./max(max(season_temp))) ; % how does the city scale with the global value?
    theta(1) = 1 ;%popsize(i);
    theta(2) = alpha(i);
    %theta(3) = beta(i);
    theta(4) = mu(i);
    theta(5) = gamma(i);
    theta(6) = R0_city(i);
    fun = @(x)findbetaSEIR(x,theta);
    x0 = [0.3];
    beta_city(i) = fminsearch(fun,x0);
    tmp = beta_city(i)./mean(season_temp(i,:));
    if(tmp>1000)
        beta_cityt(i,:) = zeros(1,365);
    else
        beta_cityt(i,:) = tmp*season_temp(i,:);
    end
end
end
