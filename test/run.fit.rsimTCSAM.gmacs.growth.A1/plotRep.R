#require('tcltk')
require('tcsam2015')

setwd('~/StockAssessments-Crab/AssessmentModelDevelopment/TCSAM2015/test/run.fit.rsimTCSAM.gmacs.growth.A1')
res.A1.0<-getModelResults();
res.A1.1<-getModelResults();

plotParameters(res);
plts.A1.0<-plotModelResults(res.A1.0,out.pdf='tcsam2015.A1.0.pdf');
plts.A1.1<-plotModelResults(res.A1.1,out.pdf='tcsam2015.A1.1.pdf');
