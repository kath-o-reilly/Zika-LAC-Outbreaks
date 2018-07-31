function [Xdat_cases,Xdat_incid,Xdat_peak_ever,Xpop] = reshape_data07(country_who,country_paho,incid,season_country,season_admn1,popsize,intro_admn1)
% outputs: Xdat_cases,Xdat_incid,Xdat_epid,Xdat_onset,Xdat_peak,Xdat_peak_ever,Xpop
% requirements: match incidence to country/state level  "intro_country"
Xdat_cases = zeros(length(country_paho),1)  ;
Xdat_incid = zeros(length(country_paho),3)  ;
Xdat_peak_ever = zeros(length(country_paho),1)  ;
Xpop = zeros(length(country_paho),1)  ;
%Xdat_epid = zeros(length(country_paho),3)  ;
%Xdat_onset = zeros(length(country_paho),3)  ;  % onset
%Xdat_peak = zeros(length(country_paho),3)  ;  % onset

for i = 1:length(country_paho)
    k = strmatch(country_paho(i),upper(incid{:,1}),'exact');  % use upper
    if(length(k)==3) % 1 country entry
        % peak_ever
        Xdat_peak_ever(i,1) = mean(incid{k,7});
        % cases
        Xdat_cases(i,1) = sum(incid{k,4});
        % incid - for each year
        oo = find(incid{k,2}==2015);
        Xdat_incid(i,1) = incid{k(oo),5}; % per 100,000
        oo = find(incid{k,2}==2016);
        Xdat_incid(i,2) = incid{k(oo),5}; % per 100,000
        oo = find(incid{k,2}==2017);
        Xdat_incid(i,3) = incid{k(oo),5}; % per 100,000
        % pop size
        Xpop(i,1) = max(incid{k,3});
    end
    m = strmatch(intro_admn1(i),upper(incid{:,8}),'exact'); % brazil & mexico States
    if(length(m)>0)
        % peak_ever
        Xdat_peak_ever(i,1) = mean(incid{m,7});
        % cases
        Xdat_cases(i,1) = sum(incid{m,4});
        % incid - for each year
        oo = find(incid{m,2}==2015);
        if(length(oo)==1)
            Xdat_incid(i,1) = incid{m(oo),5}; % per 100,000
        end
        oo = find(incid{m,2}==2016);
        if(length(oo)==1)
            Xdat_incid(i,2) = incid{m(oo),5}; % per 100,000
        end
        oo = find(incid{m,2}==2017);
        if(length(oo)==1)
            Xdat_incid(i,3) = incid{m(oo),5}; % per 100,000
        end
        % pop size
        Xpop(i,1) = max(incid{m,3});
    end
end
    

end % end function
