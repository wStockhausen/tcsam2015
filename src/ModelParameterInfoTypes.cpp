/*------------------------------------------------------------------------------
 *  Includes:
 *      NumberInfo
 *      BoundedNumberInfo
 *      VectorInfo
 *      BoundedVectorInfo
 *      DevsVectorInfo
 *      NumberVectorInfo
 *      BoundedNumberVectorInfo
 *      VectorVectorInfo
 *      BoundedVectorVectorInfo
 *      DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"
#include "ModelConfiguration.hpp"
#include "ModelParameterInfoTypes.hpp"

int NumberInfo::debug              = 0;
int BoundedNumberInfo::debug       = 0;
int VectorInfo::debug              = 0;
int BoundedVectorInfo::debug       = 0;
int DevsVectorInfo::debug          = 0;
int NumberVectorInfo::debug        = 0;
int BoundedNumberVectorInfo::debug = 0;
int VectorVectorInfo::debug        = 0;
int BoundedVectorVectorInfo::debug = 0;
int DevsVectorVectorInfo::debug    = 0;
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
    os<<"finalVal="<<finlVal<<cc;
    os<<"phase="<<phase<<cc;
    os<<"resample="<<qt<<wts::getOnOffType(resample)<<qt<<cc;
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
    writeToR1(os);
    os<<")";
}

/***************************************************************
*   writeToR1                                                    *
***************************************************************/
void BoundedNumberInfo::writeToR1(ostream& os){
    os<<"lower="<<lower<<cc;
    os<<"upper="<<upper<<cc;
    os<<"jitter="<<qt<<wts::getOnOffType(jitter)<<qt<<cc;
    NumberInfo::writeToR1(os);
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
    if (debug) cout<<"starting VectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=1;i<=N;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<endl;
        cout<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        cout<<"finished VectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
dvar_vector VectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting VectorInfo::calcLogPrior(pv) for "<<name<<endl;
    dvar_vector lps;
    if (pMPI->canCalcLogPDF(pv)) {
        lps = pMPI->calcLogPDF(pv,priorParams,priorConsts);
    } else {
        lps.allocate(pv.indexmin(),pv.indexmax());
        lps.initialize();
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = NumberInfo::calcLogPrior(x);
        }
    }
    if (debug) {
        cout<<"logPrior: "<<lps<<endl;
        cout<<"finished VectorInfo::calcLogPriors(pv) for "<<name<<endl;
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
    if (debug) cout<<"Starting VectorInfo::read(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>idxType;
    int imn; int imx; 
    tcsam::getIndexLimits(idxType,imn,imx);
    ptrIB = new IndexBlock(imn,imx);
    is>>(*ptrIB);
    N = ptrIB->getSize();
    initVals.allocate(1,N);
    is>>str; readVals = wts::getBooleanType(str);
    NumberInfo::read(is);
    initVals = initVal;
    if (debug) {
        cout<<"idxType = "<<idxType<<endl;
        cout<<"IndexBlock = "<<(*ptrIB)<<endl;
        cout<<"N     = "<<N<<endl;
        cout<<"init values = "<<initVals<<endl;
        cout<<"Finished VectorInfo::read(cifstream & is) for "<<name<<endl;
        cout<<"Enter 1 to continue>> ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
}

/***************************************************************
*   write                                                      *
***************************************************************/
void VectorInfo::write(ostream & os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readVals);
    NumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void VectorInfo::writeToR(ostream& os){
    if (debug) cout<<"VectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    NumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector
 * @param os - the output stream.
 */
void VectorInfo::writeFinalValsToR(ostream& os){
    if (debug) cout<<"VectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}
////////////////////////////VectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************\n
*   Sets initial values.       \n
***************************************************************/
void BoundedVectorInfo::setInitVals(dvector& x){
    if (debug) {
        cout<<"starting BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
        cout<<"input  x index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        cout<<"initVals index limits: "<<initVals.indexmin()<<cc<<initVals.indexmax()<<endl;
    }
    initVals = x;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
    if (debug) {
        cout<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        cout<<"finished BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}

/***************************************************************
*   read initial values. \n
***************************************************************/
void BoundedVectorInfo::readInitVals(cifstream & is){
    is>>initVals;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
}

/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector BoundedVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting BoundedVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=1;i<=N;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<endl;
        cout<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        cout<<"finished BoundedVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
dvar_vector BoundedVectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting BoundedVectorInfo::calcLogPrior(pv) for "<<name<<endl;
    dvar_vector lps;
    if (pMPI->canCalcLogPDF(pv)) {
        lps = pMPI->calcLogPDF(pv,priorParams,priorConsts);
    } else {
        lps.allocate(pv.indexmin(),pv.indexmax());
        lps.initialize();
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = BoundedNumberInfo::calcLogPrior(x);
        }
    }
    if (debug) {
        cout<<"logPrior: "<<lps<<endl;
        cout<<"finished BoundedVectorInfo::calcLogPriors(pv) for "<<name<<endl;
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
    if (debug) cout<<"Starting BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>idxType;
    int imn; int imx; 
    tcsam::getIndexLimits(idxType,imn,imx);
    ptrIB = new IndexBlock(imn,imx);
    is>>(*ptrIB);
    N = ptrIB->getSize();
    initVals.allocate(1,N);
    is>>str; 
    readVals = wts::getBooleanType(str);
    BoundedNumberInfo::read(is);
    initVals = initVal;
    if (debug) {
        cout<<"idxType    = "<<idxType<<endl;
        cout<<"IndexBlock = "<<(*ptrIB)<<endl;
        cout<<"N          = "<<N<<endl;
        cout<<"Done BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
        cout<<"Enter 1 to continue>> ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
}

/***************************************************************
*   write                                                      *
***************************************************************/
void BoundedVectorInfo::write(ostream & os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readVals);
    BoundedNumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void BoundedVectorInfo::writeToR(ostream& os){
    if (debug) cout<<"BoundedVectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    BoundedNumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector.
 * @param os - the output stream.
 */
void BoundedVectorInfo::writeFinalValsToR(ostream& os){
    if (debug) cout<<"BoundedVectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}

////////////////////////////BoundedVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************\n
*   Calculates devs.       \n
***************************************************************/
void DevsVectorInfo::calcDevs(void){
    initVals(N) = -sum(initVals(1,(N-1)));
}
/***************************************************************\n
*   Sets initial values.       \n
***************************************************************/
void DevsVectorInfo::setInitVals(dvector& x){
    if (debug) cout<<"starting DevsVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    BoundedVectorInfo::setInitVals(x);//use parent class
    calcDevs();//ensure devs
    if (debug) {
        cout<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        cout<<"finished DevsVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}
/**
 * Sets initial values 1:(N-1) to those of the vector x, but
 * sets the value for element N to -sum(initVals(1,N-1)) so
 * the sum over all elements is 0. x may have size N-1.
 * 
 * @param x - param_init_bounded_vector of initial values
 */
void DevsVectorInfo::setInitVals(param_init_bounded_vector & x){
    if (debug) {
        cout<<"starting DevsVectorInfo::setInitVals(param_init_bounded_vector & x) for "<<name<<endl;
        cout<<"input x  index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        cout<<"initVals index limits: "<<1<<cc<<N<<endl;
    }
    initVals(1,N-1) = value(x);
    calcDevs();
    if (debug) cout<<"finished DevsVectorInfo::setInitVals(param_init_bounded_vector & x) for "<<name<<endl;
}     

/**
 * Sets final values 1:(N-1) to those of the vector x, but
 * sets the value for element N to -sum(initVals(1,N-1)) so
 * the sum over all elements is 0. x may have size N-1.
 * 
 * @param x - param_init_bounded_vector of final values
 */
void DevsVectorInfo::setFinalVals(param_init_bounded_vector & x){
    if (debug) {
        cout<<"starting DevsVectorInfo::setFinalVals(param_init_bounded_vector & x) for "<<name<<endl;
        cout<<"input x  index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        cout<<"initVals index limits: "<<1<<cc<<N<<endl;
    }
    if (!finlVals.allocated()) finlVals.allocate(initVals.indexmin(),initVals.indexmax());
    finlVals(1,N-1) = value(x);
    finlVals(N) = -sum(finlVals(1,N-1));
    if (debug) cout<<"finished DevsVectorInfo::setFinalVals(param_init_bounded_vector & x) for "<<name<<endl;
}     

/***************************************************************
*   read initial values. \n
***************************************************************/
void DevsVectorInfo::readInitVals(cifstream & is){
    BoundedVectorInfo::readInitVals(is);
    calcDevs();
}

/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector DevsVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) cout<<"starting DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl = BoundedVectorInfo::drawInitVals(rng,vif);
    smpl(N) = -sum(smpl(1,(N-1)));
    if (debug) {
        cout<<phase<<tb<<pMPI->canSample()<<endl;
        cout<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        cout<<"finished DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + BoundedNumberInfo::read(is)
***************************************************************/
void DevsVectorInfo::read(cifstream & is){
    if (debug) cout<<"Starting DevsVectorInfo::read(cifstream & is) for "<<name<<endl;
    BoundedVectorInfo::read(is);
    initVals = 0.0;
    if (debug) {
        cout<<"idxType    = "<<idxType<<endl;
        cout<<"IndexBlock = "<<(*ptrIB)<<endl;
        cout<<"N          = "<<N<<endl;
        cout<<"Done DevsVectorInfo::read(cifstream & is) for "<<name<<endl;
        cout<<"Enter 1 to continue>> ";
        cin>>debug;
        if (debug<0) exit(-1);
    }
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void DevsVectorInfo::writeToR(ostream& os){
    if (debug) cout<<"DevsVectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    BoundedNumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector.
 * @param os - the output stream.
 */
void DevsVectorInfo::writeFinalValsToR(ostream& os){
    if (debug) cout<<"DevsVectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}

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
*   Set initial parameter values.                              *
***************************************************************/
void NumberVectorInfo::setInitVals(param_init_number_vector& x){
    if (debug) cout<<"starting NumberVectorInfo::setInitVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setInitVal(x(i));
    if (debug) cout<<"finished NumberVectorInfo::setInitVals(x)"<<this<<endl;
}

/***************************************************************
*   Set final parameter values.                              *
***************************************************************/
void NumberVectorInfo::setFinalVals(param_init_number_vector& x){
    if (debug) cout<<"starting NumberVectorInfo::setFinalVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setFinalVal(x(i));
    if (debug) cout<<"finished NumberVectorInfo::setFinalVals(x)"<<this<<endl;
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
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
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

/**
 * Writes final values to stream as an R list.
 * @param os - the stream to write to
 */
void NumberVectorInfo::writeFinalValsToR(ostream& os){
    if (nNIs){
        os<<"list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->getFinalVal(); os<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->getFinalVal(); os<<")";
    } else {
        os<<"NULL";
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
*   Set initial parameter values.                              *
***************************************************************/
void BoundedNumberVectorInfo::setInitVals(param_init_bounded_number_vector& x){
    if (debug) cout<<"starting BoundedNumberVectorInfo::setInitVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->setInitVal(x(i));
    if (debug) cout<<"finished BoundedNumberVectorInfo::setInitVals(x)"<<this<<endl;
}

/***************************************************************
*   Set final parameter values.                              *
***************************************************************/
void BoundedNumberVectorInfo::setFinalVals(param_init_bounded_number_vector& x){
    if (debug) cout<<"starting BoundedNumberVectorInfo::setFinalVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->setFinalVal(x(i));
    if (debug) cout<<"finished BoundedNumberVectorInfo::setFinalVals(x)"<<this<<endl;
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
    if (debug) cout<<"starting BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nNIs;
    if (debug) cout<<"nNIs ="<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        int idx;
        ppNIs = new NumberInfo*[nNIs];
        for (int p=0;p<nNIs;p++) {
            is>>idx;
            if (idx<=nNIs){
                BoundedNumberInfo* pBNI = new BoundedNumberInfo();
                is>>(*pBNI);
                ppNIs[idx-1] = pBNI;
            } else {
                cout<<"Error reading file "<<is.get_file_name()<<endl;
                cout<<"for bounded parameter '"<<name<<"' defined for "<<nNIs<<" values."<<endl;
                cout<<"Tried to read "<<idx<<"th bounded number";
            }
        }
        if (debug) {
            for (int p=0;p<nNIs;p++) cout<<p+1<<tb<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedNumberVectorInfo::write(ostream & os){
    if (debug) cout<<"Starting BoundedNumberVectorInfo::write(ostream & os) for "<<name<<endl;
    os<<tb<<nNIs<<"  #number of bounded parameters"<<endl;
    os<<"#id lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
    if (debug) cout<<endl<<"Finished BoundedNumberVectorInfo::write(ostream & os)"<<endl;
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
    if (debug) cout<<"starting VectorVectorInfo::deallocate(void) for "<<name<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) cout<<"finished VectorVectorInfo::deallocate(void) for "<<name<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMinIndices(void){
    if (debug) cout<<"starting VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) cout<<"finished VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMaxIndices(void){
    if (debug) cout<<"starting VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) cout<<"finished VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector VectorVectorInfo::getPhases(void){
    if (debug) cout<<"starting VectorVectorInfo::getPhases(void) for "<<name<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) cout<<"finished VectorVectorInfo::getPhases(void) for "<<name<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector VectorVectorInfo::getPriorWgts(){
    if (debug) cout<<"starting VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) cout<<"finished VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    return wts;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void VectorVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
    is>>nVIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        ppVIs = new VectorInfo*[nVIs];
        int idx;
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                ppVIs[idx-1] = new VectorInfo();
                is>>(*ppVIs[idx-1]);
            } else {
                cout<<"Error in VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
                cout<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                cout<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                cout<<"Aborting..."<<endl;
                exit(-1);
            }
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) cout<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void VectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of parameters"<<endl;
    os<<"#id idx.type  idx.block  read? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
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

/**
 * Writes final values to an output stream as an R list of vector.
 * 
 * @param os - the output stream
 */
void VectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
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
    if (ppVIs) for (int i=1;i<=nVIs;i++) lbs(i) = (ppVIs[i-1])->getLowerBound();
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
    if (ppVIs) for (int i=1;i<=nVIs;i++) ubs(i) = (ppVIs[i-1])->getUpperBound();
    if (debug) cout<<"finished BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void BoundedVectorVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nVIs;
    if (debug) cout<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        int idx;
        ppVIs = new BoundedVectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                BoundedVectorInfo* pBNI = new BoundedVectorInfo();
                is>>(*pBNI);
            ppVIs[idx-1] = pBNI;
            } else {
                cout<<"Error in BoundedVectorVectorInfo::read(cifstream & is)"<<endl;
                cout<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                cout<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                cout<<"Aborting..."<<endl;
                exit(-1);
            }
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) cout<<p+1<<tb<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of bounded parameters"<<endl;
    os<<"#id   idx.block   read? lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
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

/**
 * Writes final values to an output stream as an R list of vectors.
 * 
 * @param os - the output stream
 */
void BoundedVectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
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
    if (ppVIs) idxs=1;
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
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
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
////////////////////////////BoundedVectorVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector DevsVectorVectorInfo::getLowerBounds(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nVIs);
    lbs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) lbs(i) = (ppVIs[i-1])->getLowerBound();
    if (debug) cout<<"finished DevsVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector DevsVectorVectorInfo::getUpperBounds(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nVIs);
    ubs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) ubs(i) = (ppVIs[i-1])->getUpperBound();
    if (debug) cout<<"finished DevsVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void DevsVectorVectorInfo::read(cifstream & is){
    if (debug) cout<<"starting DevsVectorVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nVIs;
    if (debug) cout<<"reading info for devs vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        int idx;
        ppVIs = new DevsVectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                DevsVectorInfo* pDVI = new DevsVectorInfo();
                is>>(*pDVI);
            ppVIs[idx-1] = pDVI;
            } else {
                cout<<"Error in DevsVectorVectorInfo::read(cifstream & is)"<<endl;
                cout<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                cout<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                cout<<"Aborting..."<<endl;
                exit(-1);
            }
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) cout<<p+1<<tb<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) cout<<"finished DevsVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void DevsVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of devs vectors"<<endl;
    os<<"#id   idx.block   read? lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts"<<endl;
    if (nVIs){
        for (int p=0;p<(nVIs-1);p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<nVIs<<tb<<(*ppVIs[nVIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void DevsVectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nVIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}

/**
 * Writes final values to an output stream as an R list of vectors.
 * 
 * @param os - the output stream
 */
void DevsVectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
    }
}
/***************************************************************
*   deallocation                                               *
***************************************************************/
void DevsVectorVectorInfo::deallocate(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::deallocate(void) "<<this<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) cout<<"finished DevsVectorVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector DevsVectorVectorInfo::getMinIndices(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) cout<<"finished DevsVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector DevsVectorVectorInfo::getMaxIndices(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) cout<<"finished DevsVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector DevsVectorVectorInfo::getPhases(void){
    if (debug) cout<<"starting DevsVectorVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) cout<<"finished DevsVectorVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector DevsVectorVectorInfo::getPriorWgts(){
    if (debug) cout<<"starting DevsVectorVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) cout<<"finished DevsVectorVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}
////////////////////////////DevsVectorVectorInfo/////////////////////////////////
