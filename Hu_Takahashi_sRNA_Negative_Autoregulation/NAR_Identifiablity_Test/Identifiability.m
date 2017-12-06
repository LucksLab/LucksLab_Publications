function [pset] = Identifiability (SSM)
%this function uses SSM to to determine identifiability of parameters
[r,c]=size(SSM); % r=time steps, c=15
X=[];
pset=[];


% McAuley procedure doi:10.1081/PRE-120024426
R=SSM; % the first time, use SSM to do column sumComputeIdentifiability

for j=1:c
     for i=1:c
         M(i)=R(:,i)'*R(:,i); % the square sum of each column
     end
     
     [a,pos]=max(M);          % finds the indices of the maximum values of M, and returns them in output vector pos. 
     
     if a>1.0e-08                 % this is the tolerance
         X=[X SSM(:,pos)];      % the colomn that has the largest SS magnitude
         pset=[pset; pos];      % give the index of the parameter
         Shat=X*inv(X'*X)*X'*SSM; % Find the prediction SSM
         R=SSM-Shat;            % residual mtrx, now the residual matrix is the new mtrx that we find the next identifiable parameter, return to the top of the J loop  
     end
 end

pset=unique(pset);                % identifiable parameters


end

