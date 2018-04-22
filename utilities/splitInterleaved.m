function [itrain, itest] = splitInterleaved(NT, Lblock, fractrain, rseed)

indx = ceil((1:NT)/Lblock);
Nblocks = max(indx(:));
rng(rseed);
irand = randperm(Nblocks);
Ntrain = ceil(fractrain * Nblocks);
Ntest  = Nblocks - Ntrain;

itrain = ismember(indx, irand(1:Ntrain));
itest  = ismember(indx, irand((1+Ntrain):Nblocks));
    