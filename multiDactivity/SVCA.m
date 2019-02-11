% computes a cross-validated form of PCA: Shared Variance Component
% Analysis. Components are extracted from the covariance between neuron
% sets "ntrain" and "ntest" on training timepts "itrain". The variance of
% these covariance components is computed on testing timepts "itest". 
% This variance is the amount of reliable variance of that component 
% (because it's consistent across timepts). 
% Why compute it in the covariance space?
%   we assume the covariance btw/ neurons is consistent across time 
%   - we have no other metric to track because we are not presenting stimuli
% INPUTS: 
%     Ff (neurons x timepts)
%     npc (number of PCs to compute)
%     ntrain (one half of neurons)
%     ntest (other half of neurons)
%     itrain (one half of timepts)
%     itest (other half of timepts)
% OUTPUTS:
%     sneur (shared variance of each covariance component)
%     vneur (total variance of each covariance component)
%     u (left eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)
%     v (right eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)
function [sneur, varneur, u, v, noise] = SVCA(Ff, npc, ntrain, ntest, itrain, itest)

cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
[u,~,v] = svdecon(cov);
u = u(:,1:npc);
v = v(:,1:npc);
s1 = u' * Ff(ntrain,itest);
s2 = v' * Ff(ntest,itest);
sneur = sum(s1 .* s2, 2);
varneur = sum(s1.^2 + s2.^2,2)/2;




