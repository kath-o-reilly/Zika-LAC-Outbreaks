
% this is an example fiel for the parameter estimation
% only the first section, selecting suitable parameter sets for the ABC,
% are shown
% parameter fitting is done by selecting different movement models
% the parallel toolbox is required for "parfor"

% feature('DefaultCharacterSet','UTF8')
% SECTION 1 - load in data
 
season = readtable('/Users/kathleen/Documents/GitHub/Zika-LAC-Outbreaks/Data/dat_cities_checked_R0_Nov17.csv','FileEncoding','UTF-8');
gmat = readtable('/Users/kathleen/Documents/GitHub/Zika-LAC-Outbreaks/Data/distmat_Nov17.csv','FileEncoding','UTF-8');
rmat  = readtable('/Users/kathleen/Documents/GitHub/Zika-LAC-Outbreaks/Data/rad_matrix_24Nov17.csv','FileEncoding','UTF-8');
% <data not included> fmatavg = readtable('/Users/kathleen/Dropbox (VERG)/Dropbox/Zika/ZikaData/Outputs/flights_avg_24Nov17.csv','FileEncoding','UTF-8');
fmatavg = rmat;
intro = readtable('/Users/kathleen/Documents/GitHub/Zika-LAC-Outbreaks/Data/country_onset_2Oct2017.csv','FileEncoding','UTF-8');
incid = readtable('/Users/kathleen/Documents/GitHub/Zika-LAC-Outbreaks/Data/all_countries_incidence_18Feb18.csv','FileEncoding','UTF-8');  % I don't think this needed updating for GLP, etc
 
ref = 1:90;   % 70 = the full lot
n=length(ref);%height(scale);
popsize = round(season{ref,9},0); % population size - was called "model_pop" before
intro_yr = str2double(intro{:,14});
intro_mth = str2double(intro{:,13});
country_who = upper(intro{:,3});
country_paho = upper(intro{:,27});
intro_admn1 = upper(intro{:,5});
intro_admn1ALT = intro{:,7};  % inc. latin letters
season_country = upper(season{ref,5});
season_admn1 = upper(season{ref,11});
season_city = upper(season{ref,1});  
season_temp = season{ref,13:377};
%              
% % data - order the same as "intro_country"
% % requires re-formatting of incid files
[Xdat_cases,Xdat_incid,Xdat_peak_ever,Xpop] = reshape_data07(country_who,country_paho,incid,season_country,season_admn1,popsize,intro_admn1);
Xdat_epid = Xdat_incid>0.1;  % logical
% 
% % SECTION 2 - define parameters (independent of varying parameters)
% 
% % NH parameters - these don't change
%  
N0=popsize;
X0=popsize;             % susceptible
W0=0.0*ones(n,1);             % exposed
Y0=0.0*ones(n,1); 
Ystart=0.0*ones(n,1); 
if n==90
    Ystart(61)=1/100; %    % infected for proper sims in Sao Paulo 
    Ystart(45)=1/100; %    % infected for proper sims in Matto Grosso
    Ystart(50)=1/100; %    % infected for proper sims in Maranao
    Ystart(52)=1/100; %    % infected for proper sims in Pariaba
    Ystart(85)=1/100; %    % infected for proper sims in Sergipe
end
if n~=90
    Ystart(2)=10; %X0(2)=X0(2)-Y0(2);    % infected
end
X=zeros(n,n); W=zeros(n,n); Y=zeros(n,n); NN=zeros(n,n);II=zeros(n,n);
X=diag(X0); W=diag(W0); Y=diag(Y0); NN=diag(N0); All=reshape([X W Y NN II],5*n*n,1);  % added in II so we have incidence too
  
years = 8;
MaxTime=365*years;
options = odeset('RelTol', 1e-2);
% create dummy variables
c_ones = ones(1,n); % a column of ones
r_ones = ones(n,1); % a row of ones
% beta included to pre-allocate vector
rr = 0.001;          % return rate per day
TIntro = zeros(1,n) ; %floor(365+(365*3/12)); % introduce 1 year into run and some point within that year
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ABC section - round 1
% % random profile - extract covariance structure
% % as we are examining all model options, we can have model as a loop
% % use acceptance rate as sampler for particles in the next round
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% 
Iter = 5; % number of parameter sets for each model
% 
% % uniform : a + (b-a).*rand(Iter,1)
p_R0 = 2 + (20-5).*rand(Iter,5);            %  5-11  
p_lgs = (-10 + (-1--10).*rand(Iter,5));  % -10 to -3  "gravity-simple"
p_lrs = (-10 + (-1--10).*rand(Iter,1));  % -10 to -3  "radiation-simple"
p_lfs = (-10 + (-3--10).*rand(Iter,1));  % -10 to -3  "flight-simple"
p_ga = 2 + (15-2).*rand(Iter,1);            % 2-15
p_in = 2 + (15-2).*rand(Iter,1);            % 2-15
p_mo = 70 + (80-70).*rand(Iter,1) ;         % 70-80
p_e1 = 1 + (5-1).*rand(Iter,1);             % 1-5
p_e2 = 0.01 + (1-0.01).*rand(Iter,1)  ;     % 0.01-1
p_e3 = 0.1 + (1-0.1).*rand(Iter,1);         % 0.1 - 1
p_pr = rand(Iter,1);         % 0.1 - 1
% % [2.2 0.1 0.9];
% 
% % *** now, to avoid something going wrong, we need param_block
% % *** update param_block as the iterations loop thorugh
param_block = zeros(Iter,11+5,3);
for i = 1:3
    param_block(:,1,i) = p_R0(:,i);
    param_block(:,2,i) = p_lgs(:,i);
    param_block(:,3,i) = p_lrs;
    param_block(:,4,i) = p_lfs;
    param_block(:,5,i) = p_ga;
    param_block(:,6,i) = p_in;
    param_block(:,7,i) = p_mo;
    param_block(:,8,i) = p_e1;
    param_block(:,9,i) = p_e2;
    param_block(:,10,i) = p_e3;
    param_block(:,11,i) = p_pr;
end
% 
% 
mod_opts = ["grav-simple","grav-exp","rad-simple"]; %<flight analysis excluded>,"flight-simple","flight-rad"]  % STRING (note double quotation marks)
% % at each iteration store the relvant param_block
% %Store1 =  [p_mod p_R0 p_lgs p_lrs p_lfs p_ga p_in p_mo p_e1 p_e2 p_e3 p_pr zeros(Iter,5) ] ;
%  
% % store the output below in Out_paramsI and Out_paramsP...
% 
Out_paramsI = zeros(Iter,3);
Out_paramsP = zeros(Iter,3);

% run through the loops - note we're using parfor to make use of the cores
% parallel computing toolbox required!!!

parfor mm = 1:3
%for mm = 1
    Out_incidwL = zeros(Iter,1);
    Out_peak_everwL = zeros(Iter,1);
    % make all sims grav-exp
    method = mod_opts(mm);  % STRING
    for ii = 1:Iter

        R0_max  = param_block(ii,1,mm) ;% Store(ii,1);
        if(strcmp(method,"grav-simple"))
            leave   = 10.^p_lgs(ii,mm)  %Store(ii,2); NOTE THIS IS MM!!!
            params_exp = [ 2.2 0.1 0.9] ; 
            pr = 1 ;
        end
        if(strcmp(method,"grav-exp"))
            leave   = 10.^p_lgs(ii,1)  %Store(ii,2);
            pr = 1 ;
            params_exp = [ p_e1(ii) p_e2(ii) p_e3(ii)];
        end
        if(strcmp(method,"rad-simple"))
            leave   = 10.^p_lrs(ii,1)  %Store(ii,2);
            params_exp = [ 2.2 0.1 0.9] ; 
            pr = 1 ;
        end
%         if(strcmp(method,"flight-simple"))
%             leave   = 10.^p_lfs(ii,1)  %Store(ii,2);
%             params_exp = [ 2.2 0.1 0.9] ; 
%             pr = 1 ;
%         end
%         if(strcmp(method,"flight-rad"))
%             leave   = 10.^p_lfs(ii,1)  %Store(ii,2);
%             params_exp = [ 2.2 0.1 0.9] ; 
%             pr = p_pr(ii,1) ;
%         end
        %gamma   =(1/p_ga(ii,1))*ones(n,1);
        %alpha   =(1/p_in(ii,1))*ones(n,1);
        %mu      =(1/(p_mo(ii,1)*365))*ones(n,1);
        
        gamma=(1/20)*ones(n,1);
        alpha=(1/5)*ones(n,1);
        mu=(1/(70*365))*ones(n,1); % 0*ones(n,1); %
        
        % SECTION 4 - find parameters dependent on varying parameters
        %R0 = (R0_max-1)*season{ref,2}+1;
        %                                     n,popsize,alpha,mu,gamma,R0_max,season
        [beta_city,beta_cityt,R0_city] = find_beta(n,popsize,alpha,mu,gamma,R0_max,season_temp);
        [l2,r2] = est_mig04(n,popsize,gmat,rmat,fmatavg,ref,method,leave,rr,params_exp,pr,mod_opts);
        
        dX=zeros(n,n);dW=zeros(n,n);dY=zeros(n,n);dNN=zeros(n,n);%dI=zeros(n,n);
        start = 30; % when do we introduce infection - 30 days into 2015? might be better to run for a few years so we have stable populations?
        
        % SECTION 5 - run the model
        % so we need to specifically stop the model at initiation and then
        % re-start, (having a if statement isn't good enough)
        [t, pop]=ode45(@Diff_7_1_SEIRtincid_introv7,[0 : 1 : start],All,options,n, beta_cityt, gamma, alpha, mu, l2, r2, c_ones, r_ones);
        % store model output
        t_out = t(1:end-1);
        pop_out = pop(1:end-1,:);
        % initiate infection
        % we need the correct index here!
        index_Y = diag(reshape(2*(n*n)+[1:(n*n)],n,n));
        index_I = diag(reshape(4*(n*n)+[1:(n*n)],n,n));
        pop_out(end,index_Y) = Ystart;
        [t, pop]=ode45(@Diff_7_1_SEIRtincid_introv7,[0 : 1 : MaxTime-start],pop_out(end,:),options,n, beta_cityt, gamma, alpha, mu, l2, r2, c_ones, r_ones);
        t_out = [t_out ; t+start ];
        pop_out = [ pop_out ; pop];
        
        % SECTION 6 - compare model to data
        % Y_out - daily infections (prevalence)
        % Ydat  - infections per year for GU
        % Ydat2  - incidence per year for "relevant" GU
        
        [Y_out,II,AI,Ydat,Ydat_incid,Ydat_peak_ever,Ymod_peak_ever,Ydat_weight] = comp_model_data07(n,t_out,MaxTime,pop_out,country_paho,country_who,intro_admn1,intro_admn1ALT,season_country,season_admn1,Xdat_epid,Xdat_incid,Xpop,index_Y,index_I);
        % define model fit - based on timing of peak and whether incidence
        % was above 0.1 per 100,000
        %                                
        [Out_peak_everw,Out_incidw] = fit_funct07(Xdat_incid,Xdat_peak_ever,Ydat_incid,Ydat_peak_ever,Ydat_weight);   % by just assessing how simiar 0's and 1's are...
        Out_peak_everwL(ii) =  Out_peak_everw;
        Out_incidwL(ii) =  Out_incidw;
        
        disp(ii)
    end
 
    Out_paramsI(:,mm) = Out_incidwL;
    Out_paramsP(:,mm) = Out_peak_everwL;
end

Out_total = Out_paramsI + Out_paramsP;

% % how well did the sims do?
% % from these we need to decide which particles to hold onto....
% 
% Distance = [115 ];  %  first 100 sims suggest that values <115 look like they are worth following 
% 
% p_mod = zeros(1,5);
% Store1keep = zeros(Iter*4,11+3);
% row = 1;
% for mm = 1 : 4
%     tmp = find(Out_total(:,mm)<Distance(1));
%     p_mod(mm) = length(tmp)./Iter;
%     Store1keep(row:(row+length(tmp)-1),1:11) = param_block(tmp,1:11,mm);
%     % for the fuck-up
%     Store1keep(row:(row+length(tmp)-1),2) = param_block(tmp,2,1);
%     Store1keep(row:(row+length(tmp)-1),12) = Out_total(tmp,mm);  % fit 
%     Store1keep(row:(row+length(tmp)-1),13) = mm;                  % model
%     row = row+length(tmp);
% end
% p_mod2 = p_mod./sum(p_mod); 
% toc
% 
% plot(Store1keep(:,1),Store1keep(:,2),'kx')

% save('est_params_ABCv7ALT_30Jul18mcmc.mat')
% 
% M = array2table(Store1keep);%
% filename = 'params_ABCv7ALT_30Jul2018.csv' ;  % to add varibale to filename use %s  
% writetable(M,filename,'Delimiter',',')

%load('est_params_ABCv7INIT_18Jul18mcmc.mat')

% ------------- SECTION 2 --------------------------
% take these parameter sets and run into 2018 - how many cases are we
% expecting?
% ----------------------------------------------------

% years = 5;
% MaxTime=365*years;
% 
% 
% simm = 70;
% method = mod_opts(1);  % STRING
% 
% % AI are the number of infections each year. Store 2018 values for each sim
% AI2016 = zeros(90,simm);
% AI2017 = zeros(90,simm);
% AI2018 = zeros(90,simm);
% 
% for ii = 1:simm
% 
%     R0_max  = Store1keep(ii,1);
%     leave   = 10.^Store1keep(ii,2);
%     params_exp = [ 2.2 0.1 0.9] ;
%     pr = 1 ;
% 
%     gamma=(1/20)*ones(n,1);
%     alpha=(1/5)*ones(n,1);
%     mu=(1/(70*365))*ones(n,1); % 0*ones(n,1); %
% 
%     % SECTION 4 - find parameters dependent on varying parameters
%     %R0 = (R0_max-1)*season{ref,2}+1;
%     %                                     n,popsize,alpha,mu,gamma,R0_max,season
%     [beta_city,beta_cityt,R0_city] = find_beta(n,popsize,alpha,mu,gamma,R0_max,season_temp);
%     [l2,r2] = est_mig04(n,popsize,gmat,rmat,fmatavg,ref,method,leave,rr,params_exp,pr,mod_opts);
% 
%     dX=zeros(n,n);dW=zeros(n,n);dY=zeros(n,n);dNN=zeros(n,n);%dI=zeros(n,n);
%     start = 30; % when do we introduce infection - 30 days into 2015? might be better to run for a few years so we have stable populations?
% 
%     % SECTION 5 - run the model
%     % so we need to specifically stop the model at initiation and then
%     % re-start, (having a if statement isn't good enough)
%     [t, pop]=ode45(@Diff_7_1_SEIRtincid_introv7,[0 : 1 : start],All,options,n, beta_cityt, gamma, alpha, mu, l2, r2, c_ones, r_ones);
%     % store model output
%     t_out = t(1:end-1);
%     pop_out = pop(1:end-1,:);
%     % initiate infection
%     % we need the correct index here!
%     index_Y = diag(reshape(2*(n*n)+[1:(n*n)],n,n));
%     index_I = diag(reshape(4*(n*n)+[1:(n*n)],n,n));
%     pop_out(end,index_Y) = Ystart;
%     [t, pop]=ode45(@Diff_7_1_SEIRtincid_introv7,[0 : 1 : MaxTime-start],pop_out(end,:),options,n, beta_cityt, gamma, alpha, mu, l2, r2, c_ones, r_ones);
%     t_out = [t_out ; t+start ];
%     pop_out = [ pop_out ; pop];
% 
%     % SECTION 6 - compare model to data
%     % Y_out - daily infections (prevalence)
%     % Ydat  - infections per year for GU
%     % Ydat2  - incidence per year for "relevant" GU
% 
%     [Y_out,II,AI,Ydat,Ydat_incid,Ydat_peak_ever,Ymod_peak_ever,Ydat_weight] = comp_model_data07(n,t_out,MaxTime,pop_out,country_paho,country_who,intro_admn1,intro_admn1ALT,season_country,season_admn1,Xdat_epid,Xdat_incid,Xpop,index_Y,index_I);
%     % define model fit - based on timing of peak and whether incidence
%     % was above 0.1 per 100,000
%     %
%     AI2016(:,ii) = AI(:,2);
%     AI2017(:,ii) = AI(:,3);
%     AI2018(:,ii) = AI(:,4);
%     
%     %[Out_peak_everw,Out_incidw] = fit_funct07(Xdat_incid,Xdat_peak_ever,Ydat_incid,Ydat_peak_ever,Ydat_weight);   % by just assessing how simiar 0's and 1's are...
%     %Out_peak_everwL(ii) =  Out_peak_everw;
%     %Out_incidwL(ii) =  Out_incidw;
% 
%     disp(ii)
% end
% 
%  M = array2table(AI2018);%
%  filename = 'incid2018_ABCv7ALT_30Jul2018.csv' ;  % to add varibale to filename use %s  
%  writetable(M,filename,'Delimiter',',')
%  
% 
% 
% 
