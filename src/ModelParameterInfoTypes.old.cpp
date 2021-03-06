/*------------------------------------------------------------------------------
 *  Includes:
 *      NumberInfo
 *      BoundedNumberInfo
 *      NumberVectorInfo
 *      BoundedNumberVectorInfo
 *----------------------------------------------------------------------------*/
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"
#include "ModelConfiguration.hpp"
#include "ModelParameterInfoTypes.hpp"

int NumberInfo::debug              = 0;
int BoundedNumberInfo::debug       = 0;
int VectorInfo::debug              = 1;
int BoundedVectorInfo::debug       = 1;
int NumberVectorInfo::debug        = 0;
int BoundedNumberVectorInfo::debug = 0;
int VectorVectorInfo::debug        = 1;
int BoundedVectorVectorInfo::debug = 1;
/*------------------------------------------------------------------------------
 *  NumberInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   Draw a random sample from the prior.                       * 
*   If phase<0, return initVal rather than resampling.         *
***************************************************************/
double NumberInfo::drawInitVal(random_number_generator& rng, double vif){
    if (debug) cout<<"starting NumberInfo::drawSample(random_number_generator& rng, vif)"<<this<<endl;
    double smpl  = initVal;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        smpl = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<tb<<initVal<<tb<<smpl<<endl;
        cout<<"finished NumberInfo::drawSample(random_number_generator& rng,vif)"<<this<<endl;
    }
    return smpl;
}

/******************************************************************
*   Calculate log prior probability.                              *
******************************************************************/
dvariable NumberInfo::calcLogPrior(prevariable& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting NumberInfo::calcLogPrior(prevariable& x)"<<this<<endl;
    dvariable val = pMPI->calcLogPDF(x,priorParams,priorConsts);
    if (debug) {
        if (pMPI->getNumParams()) cout<<"priorParams = "<<priorParams<<tb;
        if (pMPI->getNumConsts()) cout<<"priorConsts = "<<priorConsts<<tb;
        cout<<endl;
        cout<<"x = "<<x<<"; ln(pdf(x)) = "<<val<<endl;
        cout<<"finished NumberInfo::calcLogPrior(prevariable& x)"<<this<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    RETURN_ARRAYS_DECREMENT();
    return val;
}

/***************************************************************
*   Set prior type.                                            *
***************************************************************/
void NumberInfo::setPriorType(adstring & prior){
    if (debug) cout<<"starting NumberInfo::setPriorType(adstring prior)"<<this<<endl;
    if (pMPI) delete pMPI;
    pMPI = ModelPDFInfo::getInfo(prior);
    if (pMPI) {
        if (pMPI->getNumParams()>0) priorParams.allocate(1,pMPI->getNumParams());
        if (pMPI->getNumConsts()>0) priorConsts.allocate(1,pMPI->getNumConsts());
    }
    if (debug) cout<<"finished NumberInfo::setPriorType(adstring prior)"<<this<<endl;
}

/*********************************************\n
 * Read from the cifstream object.
 * Read order is:
 *  initVal, phase, resample, priorWgt, 
 *  priorType, priorParams, priorConsts
**********************************************/
void NumberInfo::read(cifstream & is){
    if (debug) cout<<"Starting NumberInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>initVal;
    is>>phase;
    is>>str; resample=wts::getOnOffType(str);
    is>>priorWgt;
    is>>priorType;//prior type
    if (debug){
        cout<<initVal<<tb<<"#initVal"<<endl;
        cout<<phase  <<tb<<"#phase"<<endl;
        cout<<wts::getOnOffType(resample)<<tb<<"#resample"<<endl;
        cout<<priorWgt<<tb<<"#priorWgt"<<endl;
        cout<<priorType<<tb<<"#prior type"<<endl;
    }
    setPriorType(priorType);
    if (pMPI->getNumParams()) {
        is>>priorParams;
        if (debug) cout<<priorParams<<tb<<"#prior params"<<endl;
    }
    if (pMPI->getNumConsts()) {
        is>>priorConsts;
        if (debug) cout<<priorConsts<<tb<<"#prior consts"<<endl;
    }
    if (debug) cout<<"Done NumberInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
 * Write to ofstream object.\n
 * Write order is:\n
 *   initVal, phase, resample, priorWgt, \n
 *   priorType, priorParams, priorConsts\n
 ***************************************************************/
void NumberInfo::write(ostream & os){
    os<<initVal<<tb;
    os<<phase<<tb;
    os<<wts::getOnOffType(resample)<<tb;
    os<<priorWgt<<tb;
    os<<priorType<<tb;
    if (pMPI->getNumParams()) os<<priorParams<<tb<<tb;
    if (pMPI->getNumConsts()) os<<priorParams<<tb;
    adstring_array pNames = pMPI->getNamesForParams();
    adstring_array cNames = pMPI->getNamesForConsts();
    os<<"#";
    for (int p=pNames.indexmax();p<=pNames.indexmax();p++) os<<pNames(p)<<tb;
    for (int c=cNames.indexmax();c<=cNames.indexmax();c++) os<<cNames(c)<<tb;
}

/***************************************************************
*   writeToR.
***************************************************************/
void NumberInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os);
    os<<")";
}

/****************************************************************\n
*   writeToR1.
***************************************************************/
void NumberInfo::writeToR1(ostream& os){
    os<<"initVal="<<initVal<<cc;
    os<<"phase="<<phase<<cc;
    os<<"resample="<<wts::getOnOffType(resample)<<cc;
    os<<"priorWgt="<<priorWgt<<cc;
    os<<"pdfType=list(type='"<<priorType<<"',";
    if (pMPI->getNumParams()) {
        int N=pMPI->getNumParams();
        adstring_array names=pMPI->getNamesForParams();
        os<<"params=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorParams(n)<<",";os<<names(N)<<"="<<priorParams(N)<<"),";
    } else {os<<"params=NULL,";}
    if (pMPI->getNumConsts()) {
        int N=pMPI->getNumConsts();
        adstring_array names=pMPI->getNamesForConsts();
        os<<"consts=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorConsts(n)<<",";os<<names(N)<<"="<<priorConsts(N)<<")";
    } else {os<<"consts=NULL)";}
}
////////////////////////////NumberInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedNumberInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   Draw a random sample from the prior.                       * 
*   If phase<0, return initVal rather than resampling.         *
***************************************************************/
void BoundedNumberInfo::setInitVal(double x){
    if (debug) cout<<"starting BoundedNumberInfo::setInitVal(double x)"<<this<<endl;
    if (x<lower) {initVal = lower+(upper-lower)/1000000.0;} else
    if (x>upper) {initVal = upper-(upper-lower)/1000000.0;} else
    {initVal=x;}
    if (debug) cout<<"finished BoundedNumberInfo::setInitVal(double x)"<<this<<endl;
}

/***************************************************************
*   Draw a random sample from the prior.                       * 
*   If phase<0, return initVal rather than resampling.         *
***************************************************************/
double BoundedNumberInfo::drawInitVal(random_number_generator& rng, double vif){
    if (debug) cout<<"starting BoundedNumberInfo::drawSample(random_number_generator& rng)"<<this<<endl;
    double smpl  = initVal;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        smpl = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<tb<<initVal<<tb<<smpl<<endl;
        cout<<"finished BoundedNumberInfo::drawSample(random_number_generator& rng)"<<this<<endl;
    }
    return smpl;
}

/***************************************************************
* Read from cifstream object.\n
 * Read order is:\n
 *  lower, upper, jitter + NumberInfo::read(is)
***************************************************************/
void BoundedNumberInfo::read(cifstream & is){
    if (debug) cout<<"Starting BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>lower;
    is>>upper;
    is>>str; jitter=wts::getOnOffType(str);
    if (debug){
        cout<<lower<<tb<<"#lower"<<endl;
        cout<<upper<<tb<<"#upper"<<endl;
        cout<<wts::getOnOffType(jitter)<<tb<<"#jitter"<<endl;
    }
    NumberInfo::read(is);//finish reading
    if (debug) cout<<"Done BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
 * Write to ostream object.\n
 * Write order is:\n
 *  lower, upper, jitter + NumberInfo::write(is)\n
***************************************************************/
void BoundedNumberInfo::write(ostream & os){
    os<<lower<<tb;
    os<<upper<<tb;
    os<<wts::getOnOffType(jitter)<<tb;
    NumberInfo::write(os);//finish writing
}

/***************************************************************
*   writeToR.
***************************************************************/
void BoundedNumberInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os); os<<cc;
    NumberInfo::writeToR1(os);
    os<<")";
}

/***************************************************************
*   writeToR1                                                    *
***************************************************************/
void BoundedNumberInfo::writeToR1(ostream& os){
    os<<"lower="<<lower<<cc;
    os<<"upper="<<upper<<cc;
    os<<"jitter="<<wts::getOnOffType(jitter);
}
////////////////////////////BoundedNumberInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  VectorInfo\n
 *----------------------------------------------------------------------------*/
/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector VectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting VectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector smpl  = initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=mni;i<=mxi;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<endl;
        cout<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        cout<<"finished VectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    }
    return smpl;
}

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
dvar_vector VectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting VectorInfo::calcLogPrior(pv)"<<this<<tb<<name<<endl;
    dvar_vector lps(pv.indexmin(),pv.indexmax());
    lps.initialize();
    dvariable x;
    for (int i=pv.indexmin();i<=pv.indexmax();i++) {
        x = pv(i);//have to copy pv(i) )
        lps(i) = NumberInfo::calcLogPrior(x);
    }
    if (debug) {
        cout<<"logPrior: "<<lps<<endl;
        cout<<"finished VectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + NumberInfo::read(is)
***************************************************************/
void VectorInfo::read(cifstream & is){
    if (debug) cout<<"Starting VectorInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>mni;
    is>>mxi;
    is>>str;
    readVals = wts::getBooleanType(str);
    NumberInfo::read(is);
    int imn = mni; if (imn<0) imn = ModelConfiguration::mnYr;
    int imx = mxi; if (imx<0) imx = ModelConfiguration::mxYr;
    initVals.allocate(imn,imx);
    initVals = initVal;
    if (debug) cout<<"Done VectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   write                                                      *
***************************************************************/
void VectorInfo::write(ostream & os){
    os<<mni<<tb;
    os<<mxi<<tb;
    os<<wts::getBooleanType(readVals);
    NumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void VectorInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os); os<<cc;
    NumberInfo::writeToR1(os); os<<cc;
    wts::writeToR(os,initVals);
    os<<")";
}

/***************************************************************
*   writeToR1                                                    *
***************************************************************/
void VectorInfo::writeToR1(ostream& os){
    os<<"min.idx="<<initVals.indexmin()<<cc;
    os<<"max.idx="<<initVals.indexmax();
}
////////////////////////////VectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************\n
*   Sets initial values.       \n
***************************************************************/
void BoundedVectorInfo::setInitVals(dvector& x){
    if (debug) cout<<"starting BoundedVectorInfo::setInitVals(dvector& x)"<<this<<endl;
    initVals = x;
    for (int i=mni;i<=mxi;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
    if (debug) {
        cout<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        cout<<"finished BoundedVectorInfo::setInitVals(dvector& x)"<<this<<endl;
    }
}

/***************************************************************
*   read initial values. \n
***************************************************************/
void BoundedVectorInfo::readInitVals(cifstream & is){
    is>>initVals;
    for (int i=mni;i<=mxi;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
}

/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector BoundedVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting BoundedVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector smpl  = initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=mni;i<=mxi;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<endl;
        cout<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        cout<<"finished BoundedVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    }
    return smpl;
}

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
dvar_vector BoundedVectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting BoundedVectorInfo::calcLogPrior(pv)"<<this<<tb<<name<<endl;
    dvar_vector lps(pv.indexmin(),pv.indexmax());
    lps.initialize();
    dvariable x;
    for (int i=pv.indexmin();i<=pv.indexmax();i++) {
        x = pv(i);//have to copy pv(i) )
        lps(i) = BoundedNumberInfo::calcLogPrior(x);
    }
    if (debug) {
        cout<<"logPrior: "<<lps<<endl;
        cout<<"finished BoundedVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + BoundedNumberInfo::read(is)
***************************************************************/
void BoundedVectorInfo::read(cifstream & is){
    if (debug) cout<<"Starting BoundedVectorInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>mni;
    is>>mxi;
    is>>str;
    readVals = wts::getBooleanType(str);
    BoundedNumberInfo::read(is);
    int imn = mni; if (imn<0) imn = ModelConfiguration::mnYr;
    int imx = mxi; if (imx<0) imx = ModelConfiguration::mxYr;
    initVals.allocate(imn,imx);
    initVals = initVal;
    if (debug) cout<<"Done BoundedVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   write                                                      *
***************************************************************/
void BoundedVectorInfo::write(ostream & os){
    os<<mni<<tb;
    os<<mxi<<tb;
    os<<wts::getBooleanType(readVals);
    BoundedNumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void BoundedVectorInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os); os<<cc;
    BoundedNumberInfo::writeToR1(os); os<<cc;
    wts::writeToR(os,initVals);
    os<<")";
}

/***************************************************************
*   writeToR1                                                    *
***************************************************************/
void BoundedVectorInfo::writeToR1(ostream& os){
    os<<"min.idx="<<initVals.indexmin()<<cc;
    os<<"max.idx="<<initVals.indexmax();
}
////////////////////////////BoundedVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  NumberVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   deallocation                                               *
***************************************************************/
void NumberVectorInfo::deallocate(void){
    if (debug) cout<<"starting NumberVectorInfo::deallocate(void) "<<this<<endl;
    if (ppNIs) {
        for (int p=0;p<nNIs;p++) if (ppNIs[p]!=0) delete ppNIs[p];
        delete[] ppNIs;
        ppNIs = 0;
    }
    if (debug) cout<<"finished NumberVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector NumberVectorInfo::getPhases(void){
    if (debug) cout<<"starting NumberVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nNIs);
    phases.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) phases(i) = ppNIs[i-1]->getPhase();
    if (debug) cout<<"finished NumberVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector NumberVectorInfo::getPriorWgts(){
    if (debug) cout<<"starting NumberVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nNIs);
    for (int i=1;i<=nNIs;i++) wts(i) = ppNIs[i-1]->getPriorWgt();
    if (debug) cout<<"finished NumberVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}

/***************************************************************
*   Get initial parameter values.                              *
***************************************************************/
dvector NumberVectorInfo::getInitVals(){
    if (debug) cout<<"starting NumberVectorInfo::getInitVals()"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->getInitVal();
    if (debug) cout<<"finished NumberVectorInfo::getInitVals()"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector NumberVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting NumberVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->drawInitVal(rng,vif);
    if (debug) cout<<"finished NumberVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Calculate log prior probability for all parameters.        *
***************************************************************/
dvar_vector NumberVectorInfo::calcLogPriors(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting NumberVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
    dvar_vector lps(pv.indexmin(),pv.indexmax());
    lps.initialize();
    if (ppNIs) {
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmin();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = ppNIs[i-1]->calcLogPrior(x);
        }
    }
    if (debug) {
        cout<<"logPriors: "<<lps<<endl;
        cout<<"finished NumberVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void NumberVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting NumberVectorInfo::read(cifstream & is) "<<this<<endl;
    is>>nNIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        ppNIs = new NumberInfo*[nNIs];
        int idx;
        for (int p=0;p<nNIs;p++) {
            is>>idx;
            ppNIs[idx] = new NumberInfo();
            is>>(*ppNIs[idx]);
        }
        if (debug) {
            for (int p=0;p<nNIs;p++) cout<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished NumberVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void NumberVectorInfo::write(ostream & os){
    os<<tb<<nNIs<<"  #number of parameters"<<endl;
    os<<"#id init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void NumberVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nNIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////NumberVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedNumberVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector BoundedNumberVectorInfo::getLowerBounds(void){
    if (debug) cout<<"starting BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nNIs);
    lbs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) lbs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getLowerBound();
    if (debug) cout<<"finished BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector BoundedNumberVectorInfo::getUpperBounds(void){
    if (debug) cout<<"starting BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nNIs);
    ubs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) ubs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getUpperBound();
    if (debug) cout<<"finished BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng,vif)"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->drawInitVal(rng,vif);
    if (debug) cout<<"finished BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng,vif)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void BoundedNumberVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting BoundedNumberVectorInfo::read(cifstream & is) "<<this<<endl;
    is>>nNIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        int idx;
        ppNIs = new NumberInfo*[nNIs];
        for (int p=0;p<nNIs;p++) {
            is>>idx;
            BoundedNumberInfo* pBNI = new BoundedNumberInfo();
            is>>(*pBNI);
            ppNIs[idx-1] = pBNI;
        }
        if (debug) {
            for (int p=0;p<nNIs;p++) cout<<p+1<<tb<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished BoundedNumberVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedNumberVectorInfo::write(ostream & os){
    os<<tb<<nNIs<<"  #number of bounded parameters"<<endl;
    os<<"#id lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void BoundedNumberVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nNIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////BoundedNumberVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  VectorVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   deallocation                                               *
***************************************************************/
void VectorVectorInfo::deallocate(void){
    if (debug) cout<<"starting VectorVectorInfo::deallocate(void) "<<this<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) cout<<"finished VectorVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMinIndices(void){
    if (debug) cout<<"starting VectorVectorInfo::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getMinIndex();
    if (debug) cout<<"finished VectorVectorInfo::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMaxIndices(void){
    if (debug) cout<<"starting VectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getMaxIndex();
    if (debug) cout<<"finished VectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector VectorVectorInfo::getPhases(void){
    if (debug) cout<<"starting VectorVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) cout<<"finished VectorVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector VectorVectorInfo::getPriorWgts(){
    if (debug) cout<<"starting VectorVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) cout<<"finished VectorVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}

/***************************************************************
*   Get initial parameter values.                              *
***************************************************************/
dvector VectorVectorInfo::getInitVals(){
    if (debug) cout<<"starting VectorVectorInfo::getInitVals()"<<this<<endl;
    dvector initVals(1,nVIs);
    for (int i=1;i<=nVIs;i++) initVals(i) = ppVIs[i-1]->getInitVal();
    if (debug) cout<<"finished VectorVectorInfo::getInitVals()"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector VectorVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting VectorVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector initVals(1,nVIs);
    for (int i=1;i<=nVIs;i++) initVals(i) = ppVIs[i-1]->drawInitVal(rng,vif);
    if (debug) cout<<"finished VectorVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Calculate log prior probability for all parameters.        *
***************************************************************/
dvar_matrix VectorVectorInfo::calcLogPriors(dvar_matrix & pm){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting VectorVectorInfo::calcLogPriors(pm)"<<this<<tb<<name<<endl;
    dvar_matrix lps(pm.indexmin(),pm.indexmax());
    if (ppVIs) {
        for (int i=pm.indexmin();i<=pm.indexmin();i++) {
            dvar_vector x = pm(i);//have to copy pv(i) )
            lps(i) = ppVIs[i-1]->calcLogPrior(x);
        }
    }
    if (debug) {
        cout<<"logPriors: "<<lps<<endl;
        cout<<"finished VectorVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void VectorVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting VectorVectorInfo::read(cifstream & is) "<<this<<endl;
    is>>nVIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        ppVIs = new VectorInfo*[nVIs];
        int idx;
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            ppVIs[idx] = new VectorInfo();
            is>>(*ppVIs[idx]);
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) cout<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished VectorVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void VectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of parameters"<<endl;
    os<<"#id imin imax read? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nVIs){
        for (int p=0;p<(nVIs-1);p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<nVIs<<tb<<(*ppVIs[nVIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void VectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nVIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////VectorVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedVectorVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector BoundedVectorVectorInfo::getLowerBounds(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nVIs);
    lbs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) lbs(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getLowerBound();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector BoundedVectorVectorInfo::getUpperBounds(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nVIs);
    ubs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) ubs(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getUpperBound();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void BoundedVectorVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting BoundedVectorVectorInfo::read(cifstream & is) "<<this<<endl;
    is>>nVIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        int idx;
        ppVIs = new BoundedVectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            BoundedVectorInfo* pBNI = new BoundedVectorInfo();
            is>>(*pBNI);
            ppVIs[idx-1] = pBNI;
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) cout<<p+1<<tb<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished BoundedVectorVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of bounded parameters"<<endl;
    os<<"#id imin imax read? lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nVIs){
        for (int p=0;p<(nVIs-1);p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<nVIs<<tb<<(*ppVIs[nVIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void BoundedVectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nVIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}
/***************************************************************
*   deallocation                                               *
***************************************************************/
void BoundedVectorVectorInfo::deallocate(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::deallocate(void) "<<this<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) cout<<"finished BoundedVectorVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector BoundedVectorVectorInfo::getMinIndices(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getMinIndex();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector BoundedVectorVectorInfo::getMaxIndices(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getMaxIndex();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector BoundedVectorVectorInfo::getPhases(void){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector BoundedVectorVectorInfo::getPriorWgts(){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}

/***************************************************************
*   Get initial parameter values.                              *
***************************************************************/
dvector BoundedVectorVectorInfo::getInitVals(){
    if (debug) cout<<"starting BoundedVectorVectorInfo::getInitVals()"<<this<<endl;
    dvector initVals(1,nVIs);
    for (int i=1;i<=nVIs;i++) initVals(i) = ppVIs[i-1]->getInitVal();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getInitVals()"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector BoundedVectorVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting BoundedVectorVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector initVals(1,nVIs);
    for (int i=1;i<=nVIs;i++) initVals(i) = ppVIs[i-1]->drawInitVal(rng,vif);
    if (debug) cout<<"finished BoundedVectorVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Calculate log prior probability for all parameters.        *
***************************************************************/
dvar_matrix BoundedVectorVectorInfo::calcLogPriors(dvar_matrix & pm){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting BoundedVectorVectorInfo::calcLogPriors(pm)"<<this<<tb<<name<<endl;
    dvar_matrix lps(pm.indexmin(),pm.indexmax());
    if (ppVIs) {
        for (int i=pm.indexmin();i<=pm.indexmin();i++) {
            dvar_vector x = pm(i);//have to copy pv(i) )
            lps(i) = ppVIs[i-1]->calcLogPrior(x);
        }
    }
    if (debug) {
        cout<<"logPriors: "<<lps<<endl;
        cout<<"finished BoundedVectorVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}
////////////////////////////BoundedVectorVectorInfo/////////////////////////////////
