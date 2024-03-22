1% runs all the analyses for the EPHYS recordings!
function masterAnalysis(matroot)

mstr = {'Krebs','Waksman','Robbins'};

tdelay = [-8:-3 -2 -1.5 -1:.1:1 1.5 2 3:8];
tsmooth = 1;
tbin=1.2;
Lam = [50 150 200];
useGPU=1;

%%
tscale  = 6; % number of bins at 30Hz. 3 = 100ms binning

behAnalysis = 1;
behAnalysisShort = 1;
neurAnalysis = 1;
pcAnalysis = 1;
svcAnalysis = 1;

%iupc = [];
allLoc=[];
for d = [1:3]
	%%
    mouse_name = mstr{d};

    dat = load(fullfile(matroot, sprintf('%swithFaces_KS2.mat',mouse_name)));
    
    %
    y = bin2d(single(dat.stall),tbin*dat.srate,2);
	allLoc=[allLoc;dat.brainLoc];
    
	%%
   if pcAnalysis
       [pc1,upc,mot1,corr_pc1] = pc1Analysis(y,dat.motSVD,dat.tspont,dat.tVid);
       results.pc1{d} = pc1;
       results.upc{d} = upc;
       results.mot1{d} = mot1;
	   results.corr_pc1{d} = corr_pc1;
   end
    
   % SVC + behavior analysis - lump the whole recording together instead of
   % predicting single neurons (at 30 ms time bins)
   if svcAnalysis
	   tlag = [-8:-2 -1:2/30:1 2:8];
	  [cov_neur_t, var_neur, cov_res_beh_t] = ...
		  timeDelayAnalysis(dat,tlag);
	  results.cov_neur_t(:,:,d) = cov_neur_t;
	  results.var_neur(:,d) = var_neur;
	  results.cov_res_beh_t(:,:,:,:,d) = cov_res_beh_t;
	  results.tlag = tlag;
   end
   
    % predict from faces
    if behAnalysis
        [expv,resid,totvar,ndims0,isort,ytest,ypred] = faceAnalysis(dat, tbin*dat.srate, tdelay, Lam(d));
        results.isort{d} = isort;
		results.ytest{d} = ytest;
		results.ypred{d} = ypred;
        results.expv_behavior{d} = expv;
        results.resid_behavior{d} = resid;
        results.totvar_behavior{d} = totvar;
        results.ndims0 = ndims0;
        results.tdelay = tdelay;
    end
    % high-pass filter faces
    if behAnalysisShort     
        [expv,resid,totvar,ndims0] = faceAnalysis(dat, tscale, tdelay, Lam(d));
        
        results.expv_bshort{d} = expv;
        results.resid_bshort{d} = resid;
        results.totvar_bshort{d} = totvar;
        results.tscale(d) = tscale;
        results.Lam(d) = Lam(d);
    end
    %
    if neurAnalysis
        dmin = 5;
        [expv, resid, totvar, nPCs] = peerAnalysis(y, dat.iprobe, dat.Wh, dmin);
        results.expv_all{d} = expv;
        results.nPCs = nPCs;
        results.resid{d} = resid;
        results.totvar{d} = totvar;
	end
    

	results.spks{d} = gather_try(y);
	results.Wh{d} = dat.Wh;

	results.iprobe{d} = dat.iprobe;
	results.brainLoc{d} = dat.brainLoc;

end
%

for j = 1:14
    results.allLoc(j)=sum(allLoc==j);
end

%%
save(fullfile(matroot,'ephysResults_KS2.mat'),'results');