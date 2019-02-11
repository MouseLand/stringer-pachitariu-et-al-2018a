% decode responses to 32 stimuli (natimgs or orientations)
% SVM decoder
% itrain and itest are stim IDs
% rtrain and rtest are neural responses
function [accuracy, Model] = decoder(itrain, itest, rtrain, rtest, ori)

if nargin < 5
	ori = 0;
end

Model = fitcecoc(rtrain, itrain, 'Learners', 'linear');
y = predict(Model, rtest);

nstims = length(unique(itrain));
if ori>0
	ydiff = abs(y-itest);
	ydiff = ydiff - (ydiff>nstims/2) .* (ydiff-nstims/2);
	accuracy = mean(ydiff);
else
	accuracy = 1 - mean(abs(y-itest)>0);
end