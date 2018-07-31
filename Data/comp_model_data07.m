function[Y_out,II,AI,Ydat,Ydat_incid,Ydat_peak_ever,Ymod_peak_ever,Ydat_weight] = comp_model_data07(n,t_out,MaxTime,pop_out,country_paho,country_who,intro_admn1,intro_admn1ALT,season_country,season_admn1,Xdat_epid,Xdat_incid,Xpop,index_Y,index_I)
% [Y_out,II,AI,Ydat,Ydat_incid,Xdat_incid,Ydat_epid] = comp_model_data(n,t,MaxTime,pop,intro,intro_mth,intro_admn1,season_country,season_admn1,Xdat_epid);
    
Y_out = zeros(length(t_out),n);                       % *** infected
II = zeros(length(t_out),n);                          % *** culmulative incidence
tt = ceil((1:MaxTime)/365);
AI = zeros(n,max(tt));                            % *** total cases for each year the model is run

%Ymod_peak = zeros(n,max(tt));                     % timing of the peak for each year (including 0 for no epidemic)
Ymod_peak_ever = zeros(n,1);
%AI2 = zeros(n,max(tt));                           % *** total incidence for each year the model is run
    %      index_Y = diag(reshape(2*(n*n)+[1:(n*n)],n,n)); diagnoals are
    %      infecteds in the city rather than infecteds elsewhere
Y_out=1000*pop_out(:,index_Y);  % REAL POPULATION SCALE
II=1000*pop_out(:,index_I);     % incidence
for i=1:n
    for j = 1:max(tt)   % loop through model for each year
        oo = find(tt==j);
        if(j==1)
            AI(i,1) = max(II(oo,i));              % cases in first year
        else
            AI(i,j) = max(II(oo,i)) - sum(AI(i,1:j-1));  % cases in subsequent years
        end
    end
    % just overall model peak
    tmp1 = max(Y_out(:,i));
    aa = find(Y_out(:,i)==tmp1);
    if length(aa)==1 
        Ymod_peak_ever(i,1) = t_out(aa);
    end
    if length(aa)>1
        Ymod_peak_ever(i,1) = mean(t_out(aa));
    end
end

% incidence for these cities now needs to be converted into incidence at a
% Country / State level to enable comparison
% data are Xdat
% model output is Ydat
Ydat = zeros(length(country_who),max(tt));   % need to match the country and admin1, and be in the same order!!  *** total cases ***
Ydat_peak_ever = zeros(length(country_who),1);
Ydat_incid = zeros(length(country_who),3);
Ydat_weight = zeros(length(country_who),1);   % number of cities per GU
Ycheck = zeros(length(country_who),1);
simcheck = zeros(n,1);
for i = 1:length(country_who)
    % LA countries
    k = strmatch(country_who(i),season_country,'exact');
    if(length(k)>0)  % we have some obs that fit the country
        k2 = strmatch(country_who(i),'BRAZIL','exact');  % states of brazil
        if(length(k2)>0) % brazil States
            m = strmatch(intro_admn1(i),season_admn1);
            if(length(m)>0)
                Ydat(i,:) = sum(AI(m,:),1);
                %Ydat_peaktmp(i,:) = mean(Ymod_peak(m,:),1);  % mean value??
                Ydat_peak_ever(i,1) = mean(Ymod_peak_ever(m,1));
                Ydat_incid(i,:) = sum(AI(m,1:3)); % should sum on rows?
                Ydat_weight(i,1) = length(m);
            end
        else % not brazil
            k3 = strmatch(country_who(i),'MEXICO','exact');
            if(length(k3)>0) % brazil States
                m = strmatch(intro_admn1(i),season_admn1);
                if(length(m)>0)
                    Ydat(i,:) = sum(AI(m,:),1);
                    %Ydat_peaktmp(i,:) = mean(Ymod_peak(m,:),1);  % mean value??
                    Ydat_peak_ever(i,1) = mean(Ymod_peak_ever(m,1));
                    Ydat_incid(i,:) = sum(AI(m,1:3)); % should sum on rows?
                    Ydat_weight(i,1) = length(m);
                end
            else % other countries
                 Ydat(i,:) = sum(AI(k,:),1);
                 %Ydat_peaktmp(i,:) = mean(Ymod_peak(k,:),1); % mean value??
                 Ydat_peak_ever(i,1) = mean(Ymod_peak_ever(k,1));   % mean of the different cities
                 Ydat_incid(i,:) = sum(AI(k,1:3)); % should sum on rows - yes
                 Ydat_weight(i,1) = length(k);
            end
        end
    end
end
% oo = find(char(intro{:,25})=='1');
% Xdat_epid2 = Xdat_epid(oo,:);
% Xdat_incid2 = Xdat_incid(oo,:);
% Ydat_incid = ((Ydat(oo,:)*0.1)./Xpop(oo))*1e+06;  % incidence per million population (of aggregated cities), accounting for asmptomatic rate
% Ydat_epid = +Ydat_incid>0.5;
% Ydat_peak = mod(Ydat_peaktmp,365)./365; 

end
