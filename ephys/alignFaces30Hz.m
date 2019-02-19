function alignFaces30Hz(ephysroot, matroot)

load(fullfile(ephysroot, 'probeLocations.mat'));

areaLabels = {'FrCtx','FrMoCtx','SomMoCtx','SSCtx','V1','V2','RSP',...
    'CP','LS','LH','HPF','TH','SC','MB'};

% map from KS2 result order to Nick's probe location order
premap = [1 2 7 8 3 4 5 6];

mstr = {'Krebs','Waksman','Robbins'};
tstart = [3811 3633 3323]; % start of spontaneous activity

%%
% which mouse
for imouse = [1:3]
    probeLoc = probeLocations(imouse).probe;
    probeLoc = probeLoc(premap); % remap probes in same order as spks
    mouse_name = mstr{imouse};
        
    
    
    % load files
    beh = load(fullfile(ephysroot, sprintf('%s_face_proc.mat',mouse_name)));
	motSVD = beh.motionSVD;
	tVid = beh.times; % times of movie frames in spike reference frame
	
    load(fullfile(ephysroot, sprintf('spks%s_Feb18.mat',mouse_name)));
    
    %% extract spikes
    stall = zeros(5e3,5500,'uint8');
    ij = 0;
    maxt=0;
    Wh = [];
    iprobe=[];
    brainLoc=[];
    srate = 30; % sampling rate in Hz
    for k = 1:numel(spks)
        %%
        clu = spks(k).clu;
        st  = spks(k).st;
        st = round(st*srate);
        
        S = sparse(st, clu, ones(1, numel(clu)));
        S = uint8(full(S))';
        S = S(sort(unique(clu)),:);
        goodcells = mean(S,2)*srate > 0;
        S = S(goodcells,:);
        stall(ij+[1:size(S,1)],1:size(S,2)) = S;
        ij = ij + size(S,1);
        maxt = max(maxt, size(S,2));
        
        % height of spikes
        whp = spks(k).Wheights(sort(unique(clu)));
        whp = whp(goodcells);
        
        Wh = [Wh; whp];
        
        % where is the probe
        loc = zeros(numel(whp),1);
        lowerBorder = probeLoc(k).borders.lowerBorder;
        upperBorder = probeLoc(k).borders.upperBorder;
        acronym     = probeLoc(k).borders.acronym;
        %
        for j = 1:numel(acronym)
            whichArea = find(strcmp(areaLabels, acronym{j}));
            loc(whp >= lowerBorder(j) & whp < upperBorder(j)) = whichArea;
        end
        brainLoc = [brainLoc; loc];
        iprobe=[iprobe; k * ones(size(S,1),1)];
    end
    %%
    stall = stall(1:ij, 1:maxt);
    %
    tspont = tstart(imouse)*srate : min(floor(tVid(end)*srate), size(stall,2)-4);
    stall = stall(:,tspont);
    
    tspont = tspont / srate;
    
    % exclude drifters
    drift_ratio  = -1;
    stbin = bin2d(single(stall),srate,2);
    Slow = my_conv2(stbin, 500, 2);
    rat = min(Slow, [], 2) ./max(Slow, [],2);
    %
    
    stall = stall(rat>drift_ratio, :);
    Wh = Wh(rat>drift_ratio);
    iprobe = iprobe(rat>drift_ratio);
    brainLoc = brainLoc(rat>drift_ratio);
    size(stall)
    %
    save(fullfile(matroot, sprintf('%swithFaces_KS2.mat',mouse_name)), 'stall','Wh','iprobe',...
        'motSVD','tspont','tVid','srate','brainLoc','areaLabels');
    
end