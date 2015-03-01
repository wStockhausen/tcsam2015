#require('tcltk')
require('tcsam2015')

setwd('~/StockAssessments-Crab/AssessmentModelDevelopment/TCSAM2015/test/run.fit.rsimTCSAM.gmacs.growth.A2')
res.A2.0<-getModelResults();

plotParameters(res);
plts.A2.0<-plotModelResults(res.A2.0,out.pdf='tcsam2015.A2.0.pdf');
