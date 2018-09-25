
dataroot = '/media/carsen/DATA2/grive/10krecordings/spont_paper';
load(fullfile(dataroot,'dbstimspont.mat'));

d = 1;
db = dbs(d);
dat = load(fullfile(dataroot,...
    sprintf('stimspont_%s_%s.mat',db.mouse_name,db.date)));
if isfield(dat.stat, 'redcell')
    redcell = logical([dat.stat.redcell]);
else
    redcell = false(numel(dat.stat), 1);
end
gcell = ~redcell(:);

%%
A = compute_means(dat.stim.istim, dat.stim.resp(:,gcell(:)),2);
sigvar = corr(A(1:end-1,:,1),A(1:end-1,:,2));
sigvar = diag(sigvar);
Fsp = dat.Fsp(gcell(:),:);
%%
[~,ix] = sort(sigvar,'descend');

Fex = single(Fsp(ix, :));

save(fullfile(matroot,'exResponses.mat'),'Fex');