function [Out_peak_everw,Out_incidw]=fit_funct07(Xdat_incid,Xdat_peak_ever,Ydat_incid,Ydat_peak_ever,Ydat_weight)
    % 2 outcomes 
    % - timing of the peak                     "Out_peak_everw"
    % - whether incidence was >1 per 100,000   "Out_incidw"
    
    oo = find(Ydat_weight>0 & Xdat_peak_ever>0);      % select relevant observations (there are Y and relevant data)
    Wtmp = sqrt(Ydat_weight(oo));
    Xpeak = Xdat_peak_ever(oo)-2015;        % years post 2015
    Ypeak = Ydat_peak_ever(oo)/365;    % change to years post 2015 - and actually this then makes the scale of peak and incidence similar.
    % peak ever
    Out_peak_everw = sum(((Xpeak-Ypeak).^2).*Wtmp,'omitnan') ; 
    % incidence
    if(size(Ydat_incid,2)>3)
        Xtmp = Xdat_incid(oo,:)>1;   % change to logical
        Ytmp = Ydat_incid(oo,1:3)>1;
    else
        Xtmp = Xdat_incid(oo,:)>1;   % change to logical
        Ytmp = Ydat_incid(oo,:)>1;
    end
    
    
    Out_incidw = sum(sum(((Xtmp-Ytmp).^2)').*Wtmp','omitnan') ;
end
