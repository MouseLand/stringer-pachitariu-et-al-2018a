function fig_stimon(stimtpt,ip,pos)
if nargin==1
    ip=1;
end
if nargin<3
	pos=0.5;
end
    
plot(stimtpt,'k','linewidth',1);
text(pos,1,'visual stimuli','horizontalalignment','center','fontsize',6);
hold all;
if ip
plot(1 + [0 60*5*3],-.75 * [1 1],'k','linewidth',2);
text(0.05,-.1,'5 min','horizontalalignment','center','fontsize',6,'fontangle','normal');
end
axis tight;
axis off;
