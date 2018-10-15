function [varexp, gain] = fitAffine(R, istims, estimate_additive)

%% fit affine model to stimulus dims
%R = R ./ sqrt(sum(R.^2,2));
[ntrials,nstim] = size(R);

offset = 1 * ones(ntrials, 1);
gain   = ones(ntrials,1);

% initialize stimuli with mean response
sm = [];
for i = 1:nstim
	sm(i,:) = mean(R(istims==i,:),1);
end

% R = gain * sm + offset * soff
soff = ones(1, nstim);
sm = [sm; soff];

gain0 = gain;
offset0 = offset;

Rfit = sm(istims, :);
cost = mean((R - Rfit).^2,1);
%disp(mean(cost./var(R,1,1)));

Rrez = R;
for t = 1:10
	for i = 1:nstim
		if estimate_additive
			goff = ((sm([i nstim+1],:)*sm([i nstim+1],:)')\(sm([i nstim+1],:) * R(istims==i,:)'))';
			offset(istims==i) = goff(:,2);
		else
			goff = ((sm(i,:)*sm(i,:)') \ (sm(i,:) * Rrez(istims==i,:)'))';		
		end
		gain(istims==i) = goff(:,1);		
	end
	
	gdesign = zeros(ntrials,nstim+1);
	gdesign(:,nstim+1) = offset;
	for n = 1:nstim		
		for i = 1:nstim
			gdesign(istims==i,i) = gain(istims==i);
		end		
		xtx = gdesign' * gdesign/ntrials + 1e-4 * eye(33);
		xty = gdesign' * R(:,n)/ntrials;
		sm(:,n) = xtx \ xty;		 
		Rfit(:,n) = gdesign * sm(:,n);		
	end
	if ~estimate_additive
		Rrez = R - sm(end,:);
	end
	
	cost = mean((R - Rfit).^2,1);
	
	sm = normr(sm);
end

varexp = 1 - mean(cost./var(R,1,1));
%disp(varexp)

Rtrain = [];
Rtest = [];
RtrainFit = [];
RtestFit = [];
for isti = 1:32
    isa = find(istims==isti);
	iss = randperm(numel(isa));
	iss = isa(iss);
	ni = numel(iss);
	RtrainFit = cat(1,RtrainFit,Rfit(iss(1:floor(ni/2)),:));
	RtestFit = cat(1,RtestFit,Rfit(iss(floor(ni/2)+[1:floor(ni/2)]),:));
	Rtrain = cat(1,Rtrain,R(iss(1:floor(ni/2)),:));
	Rtest = cat(1,Rtest,R(iss(floor(ni/2)+[1:floor(ni/2)]),:));
end
vsignal = mean(mean((Rtrain-RtrainFit).*(Rtest-RtestFit)));
vsignal = mean(mean((RtrainFit).*(RtestFit)));
% %%
% ypred = gain - mean(gain);
% % ypred = offset - mean(offset);
% 
% a = (x'*x/ntrials+eye(size(x,2))*1e-2)\(x'/ntrials*ypred);
% 
% corr(x*a,ypred)
% 
% clf
% plot(ypred)
% hold all
% plot(x*a)
% hold off
% title(sprintf('correlation = %2.2f',corr(x*a,ypred)));
% legend('gain','predicted gain from face-only dimensions');
% 
% %%
% 
% x = projstim{4}-mean(projstim{4});
% 
% a = (x'*x/ntrials+eye(size(x,2))*1e-2)\(x'/ntrials*ypred);
% 
% corr(x,ypred)
% 
% clf
% plot(ypred)
% hold all
% plot(x*a)
% hold off
% title(sprintf('correlation = %2.2f',corr(x*a,ypred)));
% legend('gain','predicted gain from shared dimension');
% 
% %%
% clf;
% for k =2:32
% 	isti = [1 k];
% 	subplot(6,5,k-1);
% 	for j = 1:2
% 		plot(R(istims==isti(j),isti(1)),R(istims==isti(j),isti(2)),'.');
% 		hold all;
% 		plot(Rfit(istims==isti(j),isti(1)),Rfit(istims==isti(j),isti(2)),'k')
% 	end
% 	if k==2
% 		xlabel('stim-only 1');
% 	end
% 	ylabel(sprintf('stim-only %d',k));
% 	set(gca,'xtick',[],'ytick',[]);
% 	box off;
% 	axis square;
% end