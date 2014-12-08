#require('tcltk')
require('tcsam2015')

res<-getModelResults();

plotParameters(res);
plotModelResults(res,out.pdf='tcsam2015.pdf');
