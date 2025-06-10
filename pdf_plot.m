function [a1 b1]=pdf_plot(Z,o)

f=max(abs(Z));
g=min(Z):o:max(Z);
h1=figure;
h=histogram(Z,g,'Normalization','pdf');
a1=h.Values;
b1=h.BinEdges;
b1=(b1(1:end-1)+b1(2:end))/2;
close(h1)

end