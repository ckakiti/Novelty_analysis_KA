# Notes to self, tips/tricks

## MATLAB
- to change legend for already existing figure (if you have no access to replot it from scratch)
```
f=get(gca,'Children');
legend([f(2),f(6)],'second graph','sixth graph')
```
- to copy contents from a subplot to its own figure
```
open('figurename.fig')

subplot(1,2,2); % subplot that you want to copy
ax1 = gca;
fig1 = get(ax1,'children');

figure(2);
s1 = subplot(1,1,1);
copyobj(fig1,s1)
```
