function [xrange, yrange] = synth_plot(X, W, Z, What, What0, bhat, fname,mytitle)

figure('units','normalized','outerposition',[0 0 1 1]);
plot(X, W,'o');
hold on 
plot(X, What,'x');
plot(X, What0,'*');
subjIDs = get_subjIDs_cellstr(Z);
text(X, W, subjIDs,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(X, What, subjIDs,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(X, What0, subjIDs,'VerticalAlignment','bottom','HorizontalAlignment','right')
ylabel('PC1 of Y in SPD(3)')
xlabel('X in R')

xrange = [min(X), max(X)];
yrange = [min([min(What0), min(What), min(W)]),...
    max([max(What0), max(What), max(W)])];

y_offset = xrange*bhat;
if y_offset(:,2) < (min(What0) + max(What0))/2
    y_offset = y_offset + (min(What0) + max(What0))/2;
elseif (y_offset(:,1) > (min(What0) + max(What0))/2)
    y_offset = y_offset - (min(What0) + max(What0))/2;
end

plot(xrange, y_offset,'--b');
plot(xrange, [min(What0), max(What0)],'--r');

Y = What;
nsubjects = size(Z,2);

for i=1:nsubjects
    Xi = X(Z(:,i) ==1,1);
    Yi = Y(Z(:,i) ==1,1);
    plot(Xi, Yi, '--','Color', [0.1,0.1,0.1]);
end
%'MEM individual pattern' => MEM translated
legend('Data', 'MEM prediction','MMGLM prediction','MEM population pattern','MMGLM', 'MEM individual pattern');

title(mytitle)
hold off
saveas(gcf,fname)

% save calculated data in .mat file
% filename = 'plot_synth_data.mat';
% save(filename, 'X', 'W', 'What', 'What0', 'xrange', 'yrange', 'bhat', 'Y', 'Z');
