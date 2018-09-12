
dataroot = 'D:/grive/10krecordings/spont_paper';
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
A = compute_means(dat.stim.istim, dat.stim.resp(:,gcell(:)),70);
%%
a=tsne(dat.stim.resp(:,gcell(:)),dat.stim.istim);

%%
clf;
hold all;
cm = colormap('jet');
cm = cm(1:2:end,:);
for j = 1:32
	plot3(a(dat.stim.istim==j,1), a(dat.stim.istim==j,2),j*ones(sum(dat.stim.istim==j),1),'.','color',cm(j,:));
end

%%

sigvar = corr(A(1:end-1,:,1),A(1:end-1,:,2));
sigvar = diag(sigvar);
Fsp = dat.Fsp(gcell(:),:);
%%
[~,ix] = sort(sigvar,'descend');

Fex = single(Fsp(ix(1:2000), :));

save('exResponses.mat','Fex');