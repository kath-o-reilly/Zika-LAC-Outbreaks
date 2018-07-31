function dPop=Diff_7_1_SEIRtincid_introv7(t,pop, n, beta_cityt, gamma, alpha, mu, l, r, c_ones, r_ones)

if sum(pop<0)>0
    pop(pop<0) = 0;
    %stop
end
% bring in what was estimated from the previous time-point
X=reshape(pop(1:(n*n)),n,n); W=reshape(pop((n*n)+[1:(n*n)]),n,n); Y=reshape(pop(2*(n*n)+[1:(n*n)]),n,n); 
NN=reshape(pop(3*(n*n)+[1:(n*n)]),n,n); II=reshape(pop(4*(n*n)+[1:(n*n)]),n,n);

% set up for estimation of the current one
dX=zeros(n,n); dW=zeros(n,n); dY=zeros(n,n); dNN=zeros(n,n); dII=zeros(n,n);

% Note the different use of .* and *
% if round(mod(t,1),4)==0
%     disp(t)
% end
% from beta_cityt make a beta for each city for this timepoint
tmp = floor(mod(t,365));
beta = beta_cityt(:,tmp+1);  % so between 0-0.9999 use 1, etc...

sumY=sum(Y')'; sumNN=sum(NN')';  % for calculations, to improve efficiency....
diaY=diag(Y); diaX=diag(X);

% First do the off diagonals.
%dX = -X.*((beta.*sumY./sumNN)*c_ones) + l.*(r_ones*diaX')   - r.*X - X.*(mu*c_ones) + NN.*(mu*c_ones)  ;
dX = -X.*((beta.*sumY./sumNN)*c_ones) - X.*(mu*c_ones) + NN.*(mu*c_ones)  ;
dW = X.*((beta.*sumY./sumNN)*c_ones) - W.*(alpha*c_ones) + l.*(r_ones*diag(W)') - r.*W - W.*(mu*c_ones);
%dY = W.*(alpha*c_ones) - Y.*(gamma*c_ones) + l.*(r_ones*diaY') - r.*Y - Y.*(mu*c_ones);
dY = W.*(alpha*c_ones) - Y.*(gamma*c_ones) - Y.*(mu*c_ones);
dNN = l.*(r_ones*diag(W)') - r.*W  ;        %l.*(r_ones*diag(NN)') - r.*NN;
dII = X.*((beta.*sumY./sumNN)*c_ones);

% Now do diagonals
%DX = -diaX.*(beta.*sumY./sumNN) - sum(l)'.*diaX + sum(r.*X)' - diaX.*mu + diag(NN).*mu ;
DX = -diaX.*(beta.*sumY./sumNN) - diaX.*mu + diag(NN).*mu ;
DW = diaX.*(beta.*sumY./sumNN) - diag(W).*alpha - sum(l)'.*diag(W) + sum(r.*W)' - diag(W).*mu;
%DY = diag(W).*alpha - diaY.*gamma - sum(l)'.*diaY + sum(r.*Y)' - diaY.*mu;
DY = diag(W).*alpha - diaY.*gamma - diaY.*mu;
DNN = - sum(l)'.*diag(W) + sum(r.*W)' ;            %- sum(l)'.*diag(NN) + sum(r.*NN)';
DII = diaX.*(beta.*sumY./sumNN);

dX = dX - diag(diag(dX)) + diag(DX);
dW = dW - diag(diag(dW)) + diag(DW);
dY = dY - diag(diag(dY)) + diag(DY);
dNN = dNN - diag(diag(dNN)) + diag(DNN);
dII = dII - diag(diag(dII)) + diag(DII); 

dPop=reshape([dX dW dY dNN dII],5*n*n,1);

