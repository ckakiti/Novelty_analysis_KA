# Notes to self, tips/tricks

## MATLAB
- to change legend for already existing figure (if you have no access to replot it from scratch)
```
f=get(gca,'Children');
legend([f(2),f(6)],'second graph','sixth graph')
```
- to copy contents from a subplot to its own figure (based on [this](https://www.mathworks.com/matlabcentral/answers/101806-how-can-i-insert-my-matlab-figure-fig-files-into-multiple-subplots) post)
```
open('figurename.fig')

subplot(1,2,2); % subplot that you want to copy
ax1 = gca;
fig1 = get(ax1,'children');

figure(2);
s1 = subplot(1,1,1);
copyobj(fig1,s1)
```
- to delete a subplot (based on [this](https://www.mathworks.com/matlabcentral/answers/213341-is-it-possible-to-delete-subplots) post)
```
curr_sub = subplot(1,2,1);
delete(curr_sub)
```
- to change font for all text of all subplots

`set(findall(gcf,'-property','FontSize'),'FontSize',20)`

- to change width of all lines in current plot (based on [this](https://www.mathworks.com/matlabcentral/answers/217993-how-can-i-change-linewidth-of-all-lines-in-a-printed-figure-from-simulink) post)

```
open('untitled.fig')
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
end
```

- to change color of existing plotted data

```
h=get(gca,'children');
set(h(1),'color','r')
```

- when SVGs don't save properly:

` set(gcf, 'Renderer', 'painters'); `

- when MATLAB crashes shortly after opening .fig file:

` set(0, 'DefaultFigureRenderer', 'painters'); `

## Terminal

- to quickly get path to all .dat files nested under folder, enter this into terminal:

`find -type d -printf '%d\t%P\n' | sort -r -nk1 | cut -f2-`

- to get all .mp4 files:

`find -type f -regex '.*.mp4' -printf '%d\t%P\n' | sort -r -nk1 | cut -f2-`

- to delete all .mp4 files under current folder:

`find -type f -name '*.mp4' -delete`

- to find all .txt files with specific path

`find -type f -path '*combine*' -name '*.txt'`

- keyboard shortcut to go to beginning of line: CTRL-A


## Word

- to lock a field: CMD + fn + F11
- to unlock a field: CMD + SHIFT + fn + F11

## Acrobat

- [to resize PDF](http://khkonsulting.com/2017/03/scaling-page-content-in-adobe-acrobat-pro-dc/)
