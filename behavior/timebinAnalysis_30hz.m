function results = timebinAnalysis_30hz(dataroot, matroot, useGPU)

load(fullfile(dataroot, 'dbspont_30hz.mat'));
dall.db = dbspont_30hz;

ndims0 = [1 2 4 8 16 32 64];
lams = [1e-4 1e-3 5e-3 0.01 0.05 0.1 0.15 0.3 0.5 .8 1.0];
tlag0 = [-8:-3 -2 -1.5 -1:.1:1 1.5 2 3:8];
tbins = [1 2 4 16 36 128];
npc = 128;

clear results

%%
results.cov_neur = NaN*zeros(512,length(tbins),length(dall.db));
results.var_neur = NaN*zeros(512,length(tbins),length(dall.db));
results.cov_res_beh = NaN*zeros(128,length(ndims0),length(tbins),length(lams),length(dall.db));
results.cov_res_beh_t = NaN*zeros(128,length(ndims0),length(tlag0),length(lams),length(dall.db));
rng('default');
%%
for d = [1:length(dall.db)]
	%%
	dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
	
	%%
	
	if isfield(dat.stat,'redcell')
		Ff = dat.Fsp(~logical([dat.stat(:).redcell]), :);
		stat = dat.stat(~logical([dat.stat(:).redcell]));
	else
		Ff = dat.Fsp;
		stat = dat.stat;
	end
	med = reshape([stat.med], 2, length(stat))';
	
	Ff = Ff(sum(Ff,2)>0, :);
	
	
	[NN, NT] = size(Ff);
	
	%% divide X and Y into checkerboard and use every other square
	y = round(med(:,1));
	ymax=max(med);
	ymax = ymax(1);
	nby = floor(ymax / 16);
	ytrain = ([1:2:16]-1) * nby + [1:nby-10]';
	ytrain = ytrain(:)';
	ytrain = repmat(ytrain,NN,1);
	nt= size(ytrain,2);
	ntrain = find(sum(repmat(y,1,nt) == ytrain, 2)>0);
	ytest = ([2:2:16]-1) * nby + [1:nby-10]';
	ytest = ytest(:)';
	ytest = repmat(ytest,NN,1);
	nt = size(ytest,2);
	ntest = find(sum(repmat(y,1,nt) == ytest, 2)>0);
	
	%% bin spikes in 1.2 second bins
	
	for td = 1:length(tbins)
		tbin=tbins(td);
		[NN, NT] = size(Ff);
		clear y;
		y    = squeeze(mean(reshape(Ff(:,1:floor(NT/tbin)*tbin),...
			NN, tbin, []),2));
		y = (y - mean(y,2));
		
		% apply time delay
		[NN NT] = size(y);
		fprintf('\nrecording %d\n',d);
		disp([NN NT]);
		
		% divide time in half
		Lblock = max(1,20*30/tbin);
		fractrain = 0.5;
		[itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
		tic;
		
		% SVCA
		[sneur, varneur, u, v] = SVCA(y, min(length(ntrain),length(ntest)), ntrain, ntest, itrain, itest);
		u = u(:,1:npc);
		v = v(:,1:npc);
		
		y0 = y;
		itrain0 = itrain;
		itest0 = itest;
		
		clf;
		subplot(1,2,1);
		semilogx(sneur./varneur);
		results.cov_neur(1:length(sneur),td,d) = gather_try(sneur);
		results.var_neur(1:length(sneur),td,d) = gather_try(varneur);
		hold all;
		drawnow;
		
		x = dat.beh.face.motionSVD;
		x = x - mean(x,1);
		x = x / std(x(:,1));
		x    = x * 10;
			
		x0 = x;
		
		%%
		if td == 1
			tlag = tlag0;
		else
			tlag = 0;
		end
		
		%%
		for j = 1:numel(tlag)
			
			k = round(tlag(j) / .03);
			if k > 0
				y = y0(:, k+1:end);
				itrain = itrain0(k+1:end);
				itest  = itest0(k+1:end);
				x = x0(1:end-k,:);
			else
				itrain = itrain0(1:end-abs(k));
				itest  = itest0(1:end-abs(k));
				y = y0(:, 1:end-abs(k));
				x = x0(abs(k)+1:end,:);
			end
			if useGPU
				y = gpuArray(single(y));
			end
		
			x    = bin2d(x, tbin, 1);
			x    = x';
			ndims1 = ndims0(ndims0<=size(x,1));
			
			
			%%
			for l = 1:length(lams)
				[atrain, btrain] = CanonCor2(y(ntrain,itrain)'*u, x(:,itrain)', lams(l));
				[atest, btest] = CanonCor2(y(ntest,itrain)'*v, x(:,itrain)', lams(l));
				k=0;
				
				xtrain = btrain' * x(:,itest);
				xtest = btest' * x(:,itest);
				ytrain =  atrain' * u' * y(ntrain,itest);
				ytest =   atest' * v' * y(ntest,itest);
				
				%cc = 0.5 * (ccf(ytrain(1:16,:)', xtrain(1:16,:)', 1000) + ...
				%	ccf(ytest(1:16,:)', xtest(1:16,:)', 1000));
				
				%results.ccA(:,:,td,l,d) = cc;
				
				for n = ndims1
					k=k+1;
					vp_train     = atrain(:,1:n) * btrain(:,1:n)' * x(:,itest);
					vp_test      = atest(:,1:n) * btest(:,1:n)' * x(:,itest);
					
					
					s1 = u' * y(ntrain,itest) - vp_train;
					s2 = v' * y(ntest,itest) - vp_test;
					sout = sum(s1 .* s2, 2);
					%sout = max(0, sout);
					vars = sum((u' * y(ntrain,itest)).^2 + (v' * y(ntest,itest)).^2,2)/2;
					if k==length(ndims1)
						subplot(1,2,1),
						semilogx((results.cov_neur(1:npc,td,d)-sout)./vars)
						ypred = [u*vp_train; v*vp_test];
						vpred = vp_test;
						hold all;
						drawnow;
					end
					
					if td==1
						results.cov_res_beh_t(:,k,j,l,d) = gather_try(sout);
						if j==ceil(length(tlag)/2)
							results.cov_res_beh(:,k,td,l,d) = gather_try(sout);
						end
					else
						results.cov_res_beh(:,k,td,l,d) = gather_try(sout);
					end
				end
			end
		end
			subplot(1,2,2),
			bx=cumsum(results.cov_neur(1:npc,td,d) - results.cov_res_beh(:,:,td,l,d));
			semilogx(ndims0,bx(npc,:)/sum(results.var_neur(1:npc,td,d)))
			disp(max(bx(npc,:)/sum(results.var_neur(1:npc,td,d))))
			ylim([0 max(bx(npc,:)/sum(results.var_neur(1:npc,td,d)))]);
			drawnow;
		
	end
end

%%
results.ndims0 = ndims0;
results.lams = lams;
results.tbins = tbins;
results.npc = npc;

%%
save(fullfile(matroot,'timebins_30Hz_spont.mat'),'-struct', 'results');


%%
results = load(fullfile(matroot,'timebins_30Hz_spont.mat'));

%%
%%
clf;
mm=[];
for d=1:3
for td = 1
	for j = 1:length(tlag0)
		cov_neur0 = squeeze(results.cov_neur(1:128,td,d));
		var_neur0 = squeeze(results.var_neur(1:128,td,d));
		% put the lambda as the last index
		cov_res_beh0 = squeeze(results.cov_res_beh_t(:,end-1,j,:,d));
		[m(j),ilam] = max(nanmean(nansum((cov_neur0 - cov_res_beh0)./var_neur0,1),2));
		mm(j,d) = sum(cov_neur0 - cov_res_beh0(:,ilam)) / sum(var_neur0);
	end
end
end
plot(tlag0,mean(mm,2));
xlabel('time from behavior (s)');
ylabel('variance explained');
title('SVC 1-128 (30 ms bins)');
drawnow;
%xlim([-2 2]);
%%


%%
clf;
cc=[];
for d = 1:3
	for td = 1:length(tbins)
		cov_neur = squeeze(results.cov_neur(1:128,td,d));
		var_neur = squeeze(results.var_neur(1:128,td,d));
		% put the lambda as the last index
		cov_res_beh = squeeze(results.cov_res_beh(:,end,td,:,d));
		[~,ilam] = max(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2));
		cv;
	end
end

cm=colormap('jet');
cm=cm(1:5:end,:);
for n = 1:6
	my_subplot(2,3,n);
	for td = 1:13
		semilogx([1:1000]*tbins(td)*0.03, mean(cc(:,n,td,:),4),'color',cm(td,:))
		hold all
		if n==1
			text(0.05,0.7-.05*td,sprintf('%1.2fs',tbins(td)*0.03),'color',cm(td,:),'Units','normalized');
		end
	end
	xlim([0 300]);
	set(gca,'xtick',10.^[-1:2]);
	xlabel('time (s)');
	ylabel('cross-correlation');
	ylim([0 .5]);
	grid on;
	grid minor;
	grid minor;
	title(sprintf('RR comp %d',n));
end


%%

clf
td = 1;
for n=2:10
	semilogx([1:1000]*tbins(td)*0.03, mean(cc(:,n,td,:), 4),'color',cm(n,:))
	hold all;
end
xlim([0 1000]);

%%
clf;
plot([1:1000]*tbins(td)*0.03, mean(mean(cc(:,:,td,:), 4),2),'color',cm(n,:))
xlim([0 1])


