function [l2,r2] = est_mig04(n,popsize,gmat,rmat,fmatavg,ref,method,leave,rr,params_exp,pr,mod_opts)
% inclusion of flight data
% params_exp [2.2 0.1 0.9];
function checkInputs(n,popsize,gmat,rmat,fmatavg,ref,method,leave,rr,params_exp,pr,mod_opts,varargin)
    minArgs=12;
    maxArgs=12;
    narginchk(minArgs,maxArgs)
    fprintf('est_mig04 checked ok\n', length(varargin));
end
checkInputs(n,popsize,gmat,rmat,fmatavg,ref,method,leave,rr,params_exp,pr,mod_opts)

if(sum(strcmp(method,mod_opts))<1)
    error('Error. Check method for migration - not matching')
end

l1 = zeros(n,n);     % leave matrix
r1 = l1;
l2 = zeros(n,n);     % leave matrix
r2 = l2;             % return matrix
if(isequal(method,"grav-simple"))
    disp(method)
    for i=1:n  % rows
        for j=1:n
            if i~=j
                l1(i,j) = ((popsize(i)*popsize(j))./(gmat{i+(ref(1)-1),j+ref(1)-1}+0.01)) ;
                r1(i,j) = l1(i,j) ;%((popsize(i)*popsize(j))./(gmat{i+(ref(1)-1),j+ref(1)-1}+0.01)) ;
            end
        end
    end
    % then normalise and scale
    tmpsum = sum(sum(l1,2)); % sum the columns to return 
    l2 = leave.*(l1./tmpsum);     % sum the columns, normalise each column
    r2 = rr.*(r1./tmpsum);
end
if(isequal(method,"grav-exp"))
    disp(method)
    for i=1:n  % rows
        for j=1:n
            if i~=j
                l1(i,j) = ((popsize(i)^params_exp(1))*(popsize(j)^params_exp(2)))./((gmat{i+(ref(1)-1),j+ref(1)-1}+0.01)^params_exp(3)) ;
                r1(i,j) = l1(i,j) ;%(((popsize(i)^params_exp(1))*popsize(j)))./((gmat{i+ref(1)-1,j+ref(1)-1}+0.01)^params_exp(1)) ;
            end
        end
    end
    % then normalise and scale
    tmpsum = sum(sum(l1,1));
    l2 = leave.*(l1./tmpsum) ;     % sum the columns, normalise each column
    r2 = rr.*(r1./tmpsum) ;
end
if(isequal(method,"flight-rad"))
    disp(method)
    % composite of flightand radiation model
    % then normalise and scale
    tmp = pr.*table2array(fmatavg) + (1-pr).*table2array(rmat);
    tmpsum = sum(sum(tmp,1));
    frad = tmp./tmpsum;
    l2 = leave.*(frad);     % sum the columns, normalise each column
    r2 = rr.*(frad);
end
if(isequal(method,"flight-simple"))
    disp(method)
    % then normalise and scale
    tmp = table2array(fmatavg);
    tmpsum = sum(sum(tmp,1));
    fmatavgp = tmp./tmpsum;     % this sums to 1.00
    l2 = leave.*(fmatavgp);     % sum the columns, normalise each column
    r2 = rr.*(fmatavgp);
end
if(isequal(method,"rad-simple"))
    % then normalise and scale
    % this one is fine because it sums to 1.00
    l2 = leave.*(table2array(rmat));     % sum the columns, normalise each column
    r2 = rr.*(table2array(rmat));
end


end
