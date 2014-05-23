//--------------------------------------------------------------------------------
//  Includes:
//      DevsVectorInfo
//      ModelDevsVector
//      ModelDevsVectorVector
//      ModelDevsVectorMatrix
//--------------------------------------------------------------------------------
#include <admodel.h>
#include "ModelConstants.hpp"
#include "admbFunctions.hpp"
#include "rFunctions.hpp"
#include "ModelFunctions.hpp"
#include "ModelDevs.hpp"

int DevsVectorInfo::debug        = 0;
int ModelDevsVector::debug       = 0;
int ModelDevsVectorVector::debug = 0;
int ModelDevsVectorMatrix::debug = 0;

//--------------------------------------------------------------------------------
//          DevsVectorInfo
//--------------------------------------------------------------------------------
/***************************************************************
*   destruction                                                *
***************************************************************/
DevsVectorInfo::~DevsVectorInfo(){delete pMPI; pMPI=0; delete pMTI; pMTI=0;}
/***************************************************************
*   Set prior type                                             *
***************************************************************/
void DevsVectorInfo::setPriorType(adstring prior){
    if (debug) cout<<"starting DevsVectorInfo::setPriorType(adstring name)"<<this<<endl;
    if (pMPI) delete pMPI;
    pMPI = ModelPDFInfo::getInfo(prior);
    if (pMPI) {
        if (pMPI->getNumParams()>0) priorParams.allocate(1,pMPI->getNumParams());
        if (pMPI->getNumConsts()>0) priorConsts.allocate(1,pMPI->getNumConsts());
    }
    if (debug) cout<<"finished DevsVectorInfo::setPriorType(adstring name)"<<this<<endl;
}

/***************************************************************
*   Set transform type.                                        *
***************************************************************/
void DevsVectorInfo::setTransformType(adstring type){
    if (debug) cout<<"starting DevsVectorInfo::setTransformType(adstring type)"<<this<<endl;
    if (pMTI) delete pMTI;
    pMTI = ModelTransformInfo::getInfo(type);
    if (debug) cout<<"finished DevsVectorInfo::setTransformType(adstring type)"<<this<<endl;
}

/***************************************************************
*   Get initial values for devs on parameter scale.            *
*   Resulting vector has indices (mnIndx,mxIndx).              *
***************************************************************/
dvector DevsVectorInfo::getInitParamVals(void){
    if (debug) cout<<"starting DevsVectorInfo::getInitParamVals(void)"<<this<<tb<<name<<endl;
    dvector vals = initVals;
    for (int i=mnIndx;i<=mxIndx;i++) vals(i) = pMTI->calcTransform(initVals(i));
    vals -= mean(vals);          //make sure these are devs with mean 0
    if (debug) {
        cout<<"initVals (natural scale) = "<<initVals<<endl;
        cout<<"initVals (param scale)   = "<<vals<<endl;
        cout<<"finished DevsVectorInfo::getInitParamVals(void)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return vals;
}

/***************************************************************
*   Draw random sample for the devs vector by sampling from    *
*   the prior, then transforming the values to param space.    *
*   If phase<0, just return params-space initVals rather than  *
*   resampling.                                                *
***************************************************************/
dvector DevsVectorInfo::drawInitParamVals(random_number_generator& rng){
    if (debug) cout<<"starting DevsVectorInfo::drawInitParamVals(random_number_generator& rng)"<<this<<tb<<name<<endl;
    dvector smpls(mnIndx,mxIndx);
    if (flgResample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=mnIndx;i<=mxIndx;i++) {
            smpls(i) = pMPI->drawSample(rng,priorParams,priorConsts);
            smpls(i) = pMTI->calcTransform(smpls(i));
        }
        smpls -= mean(smpls);          //make sure samples are devs with mean 0
        if (debug) cout<<"sampled devs = "<<smpls<<endl;
    } else {
        if (debug) cout<<"cannot sample"<<endl;
        return getInitParamVals();
    }
    if (debug) cout<<"finished DevsVectorInfo::drawInitParamVals(random_number_generator& rng)"<<this<<tb<<name<<endl;
    return smpls;
}

/***************************************************************
*   Calculate log prior probability.                           *
*   Here, xv is a vector on the natural scale.                 *
***************************************************************/
dvar_vector DevsVectorInfo::calcLogPriors(dvar_vector& xv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting DevsVectorInfo::calcLogPrior(dvar_vector& xv)"<<this<<tb<<name<<endl;
    int mn = xv.indexmin(); int mx = xv.indexmax();
    dvar_vector lps(mn,mx);
    if (debug) cout<<"lps indices = "<<lps.indexmin()<<tb<<lps.indexmax()<<endl;
    for (int i=mn;i<=mx;i++) {
	    dvariable xval = xv(i);
        lps(i) = pMPI->calcLogPDF(xval,priorParams,priorConsts);
    }
    if (debug) {
        if (pMPI->getNumParams()) cout<<"priorParams = "<<priorParams<<tb;
        if (pMPI->getNumConsts()) cout<<"priorConsts = "<<priorConsts<<tb;
        cout<<endl;
        cout<<"xv          = "<<xv<<endl;
        cout<<"ln(pdf(xv)) = "<<lps<<endl;
        cout<<"finished DevsVectorInfo::calcLogPrior(dvar_vector& xv)"<<this<<tb<<name<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*   read from stream                                           *
***************************************************************/
void DevsVectorInfo::read(cifstream & is){
    adstring str;
    is>>phase;
    is>>priorWgt;
    is>>mnIndx;
    is>>mxIndx;
    initVals.allocate(mnIndx,mxIndx);
    is>>str;str.to_lower();//resample
    if (str=="true") flgResample = 1; else flgResample = 0;
    is>>str;str.to_lower();//has init vals
    if (str=="true") flgInitVals = 1; else flgInitVals = 0;
    is>>str;//transform type
    setTransformType(str);
    if (pMTI->getNumConsts()) is>>pMTI->consts;
    is>>str;//prior type
    setPriorType(str);
    if (pMPI->getNumParams()) is>>priorParams;
    if (pMPI->getNumConsts()) is>>priorConsts;
    if (flgInitVals) {
        is>>initVals;
        initVals -= mean(initVals);//ensure this is a devs vector
    } else {
        initVals = 0.0;
    }
}

/***************************************************************
*   write                                                      *
***************************************************************/
void DevsVectorInfo::write(ostream & os){
    os<<tb<<name<<tb;
    os<<phase<<tb;
    os<<priorWgt<<tb;
    os<<mnIndx<<tb;
    os<<mxIndx<<tb;
    if (flgResample) os<<"TRUE"<<tb; else os<<"FALSE"<<tb;
    if (flgInitVals) os<<"TRUE"<<tb; else os<<"FALSE"<<tb;
    os<<pMTI->getTransformType()<<tb;
    if (pMTI->getNumConsts()) os<<pMTI->consts<<tb<<tb;
    os<<pMPI->getPDFType()<<tb;
    if (pMPI->getNumParams()) os<<priorParams<<tb<<tb;
    if (pMPI->getNumConsts()) os<<priorParams<<tb;
    if (flgInitVals) os<<tb<<initVals;
    os<<"#"<<pMPI->getStringForParamsNames()<<tb<<pMPI->getStringForConstsNames();
    if (flgInitVals) os<<tb<<"init_vals"; else os<<tb<<"no_init_vals";
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void DevsVectorInfo::writeToR(ostream& os){
    os<<"list(";
    os<<"name="<<name<<cc;
    os<<"phase="<<phase<<cc;
    os<<"priorWgt="<<priorWgt<<cc;
    os<<"mnIndex="<<mnIndx<<cc;
    os<<"mxIndex="<<mxIndx<<cc;
    os<<"initVals=c("; for (int i=initVals.indexmin();i<initVals.indexmax();i++) os<<initVals(i)<<cc; os<<initVals(initVals.indexmax())<<"),";
    os<<"transformType=list(type='"<<pMTI->getTransformType()<<"',";
    if (pMTI->getNumConsts()) {
        int N=pMTI->getNumConsts();
        adstring_array names=pMTI->getNamesForConsts();
        os<<"consts=list("; for (int n=1;n<N;n++) os<<names(n)<<"="<<pMTI->consts(n)<<",";os<<names(N)<<"="<<pMTI->consts(N)<<"),";
    } else {os<<"consts=NULL),";}
    os<<"pdfType=list(type='"<<pMPI->getPDFType()<<"',";
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
    os<<")";
}
////////////////////////////DevsVectorInfo/////////////////////////////////

//--------------------------------------------------------------------------------
//          ModelDevsVector
//--------------------------------------------------------------------------------
/******************************************************************
*   Set the devs parameter vector and creates the natural-scale   *
*   devs.                                                         *
*   The input param_init_vector should have                       *
*       min index = this->getMinIndex()                           *
*       max index = this->getMaxIndex())-1                        *
******************************************************************/
void ModelDevsVector::setDevParamsVector(param_init_vector& pars){
    int mn  = getMinIndex();
    int mx  = getMaxIndex();
    int mn1 = pars.indexmin();
    int mx1 = pars.indexmax();
    if ((mn!=mn1)|| ((mx-1)!=mx1)){
        cout<<"Error: mismatch between indices in ModelDevsVector::setDevParamsVector(pars)."<<endl;
        cout<<"Expected: "<<mn<<tb<<mx-1<<endl;
        cout<<"Got     : "<<mn1<<tb<<mx1<<endl;
        cout<<"Terminating program..."<<endl;
        exit(1);
    }
    pDvs = &pars;//take address
    createNaturalScaleDevs();
}

/******************************************************************
*   Create the natural-scale devs vector.                         *
*   Once calculated (use calcNaturalScaleDevs()), this vector     *
*   will be mean zero.                                            *
******************************************************************/
void ModelDevsVector::createNaturalScaleDevs(void){
    if (pASDvs) delete pASDvs;
    pASDvs = new dvar_vector(getMinIndex(),getMaxIndex());
    pASDvs->initialize();
}

/******************************************************************
*   Calcs the natural-scale devs vector corresponding to pDvs.    *
*   The vector runs                                               *
*       pDvs->indexmin():pDvs->indexmax()+1     or                *
*       getMinIndex() to getMaxIndex()).                          *
*   The vector is zero-mean, but variance DOES NOT = 1.           *
*   createNaturalScaleDevsVector() must be called prior to        *
*   calling this function.                                        *
******************************************************************/
void ModelDevsVector::calcNaturalScaleDevs(){
    if (debug) cout<<"starting ModelDevsVector::calcNaturalScaleDevs()"<<this<<endl;
    int mn = (*pDvs).indexmin();
    int mx = (*pDvs).indexmax();
    if (debug) cout<<"pDvs = "<<(*pDvs)<<endl;
    pASDvs->initialize();
    for (int i=mn;i<=mx;i++) {
        (*pASDvs)(i) = pDVI->pMTI->calcInvTransform((*pDvs)(i));
    }
    (*pASDvs)(mx+1) = -sum((*pASDvs)(mn,mx));
    dvariable var = variance(*pASDvs);
    if (debug) {
        cout<<"var = "<<var<<".  devs = "<<(*pASDvs)<<endl;
        cout<<"finished ModelDevsVector::calcNaturalScaleDevs()"<<this<<endl;
    }
}

/***************************************************************
*   Check that the given index correspond to a dev value.      *
***************************************************************/
int ModelDevsVector::checkIndex(int i){
    if (debug) {
        cout<<"starting ModelDevsVector::checkIndex(int i) "<<this<<tb<<i<<endl;
        cout<<"min, max indices = "<<pDVI->getMinIndex()<<tb<<pDVI->getMaxIndex()<<endl;
    }
    int res = (pDVI->getMinIndex()<=i)&&(i<=pDVI->getMaxIndex());
    if (debug) cout<<"results = "<<res<<endl;
    return res;
}

/***************************************************************
*   Returns the ith devs variable corresponding to pDvs.       *
*   Valid i run from                                           *
*       pDvs->indexmin():pDvs->indexmax()+1     or             *
*       pASDvs->indexmin():pASDvs->indexmax()   or             *
*       getMinIndex() to getMaxIndex()).                       *
*   Returns 0 if index is not valid.                           *
*   calcNaturalScaleDevs() should be called prior              *
*   to calling this function.                                  *
***************************************************************/
dvariable ModelDevsVector::getNaturalScaleDev(int i){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting ModelDevsVector::getNaturalScaleDev(int i)"<<this<<endl;
    dvariable val = 0;
    if ((pASDvs->indexmin()<=i)&&(i<=pASDvs->indexmax())) {
        val = (*pASDvs)(i);
    }
    if (debug) cout<<"finished ModelDevsVector::getNaturalScaleDev(int i)"<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return val;
}

/********************************************************************
*   Returns the natural-scale devs vector corresponding to pDvs.    *
*   The vector runs                                                 *
*       pDvs->indexmin():pDvs->indexmax()+1     or                  *
*       pASDvs->indexmin():pASDvs->indexmax()   or                  *
*       getMinIndex() to getMaxIndex()).                            *
*   The returned vector is zero-mean.                               *
*   calcNaturalScaleDevs() should be called prior                   *
*   to calling this function.                                       *
********************************************************************/
dvar_vector ModelDevsVector::getNaturalScaleDevsVector(){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting ModelDevsVector::getNaturalScaleDevsVector()"<<this<<endl;
    int mnIndx = pDvs->indexmin();
    int mxIndx = pDvs->indexmax();
    dvar_vector vec(mnIndx,mxIndx+1);
    vec = (*pASDvs);
    if (debug) cout<<"finished ModelDevsVector::getNaturalScaleDevsVector()"<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return vec;
}

/***************************************************************
*   Calculate log prior probability vector over all devs.      *
*   The resulting vector runs                                  *
*       pDvs->indexmin():pDvs->indexmax()+1     or             *
*       pASDvs->indexmin():pASDvs->indexmax()   or             *
*       getMinIndex() to getMaxIndex()).                       *
*   Prior is specified relative to natural-scale devs.         *
*   calcNaturalScaleDevs() should be called prior              *
*   to calling this function.                                  *
***************************************************************/
dvar_vector ModelDevsVector::calcLogPriors(){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"starting ModelDevsVector::calcLogPrior()"<<this<<endl;
    dvar_vector lps = pDVI->calcLogPriors((*pASDvs));
    if (debug) {
        cout<<"lps(indexmin,indexmax) = "<<lps.indexmin()<<" "<<lps.indexmax()<<endl;
        cout<<"finished ModelDevsVector::calcLogPrior()"<<this<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVector::writeInfoToR(ostream& os, adstring nm, int indent){
    os<<nm<<"="; pDVI->writeToR(os);
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVector::writeValuesToR(ostream& os, adstring nm, int indent){
    if (pASDvs){
        dvector vals = value((*pASDvs));
        os<<nm<<"=structure(c("; for (int i=vals.indexmin();i<vals.indexmax();i++) os<<vals(i)<<cc; os<<vals(vals.indexmax())<<"),";
        os<<"names="<<vals.indexmin()<<":"<<vals.indexmax()<<",dim=c("<<vals.indexmax()-vals.indexmin()+1<<"))";
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////ModelDevsVector/////////////////////////////////

//--------------------------------------------------------------------------------
//          ModelDevsVectorVector
//--------------------------------------------------------------------------------
/***************************************************************
*   deallocation                                               *
***************************************************************/
void ModelDevsVectorVector::deallocate(void){
    if (debug) cout<<"starting ModelDevsVectorVector::deallocate(void) "<<this<<endl;
    if (ppDvs!=0) {
        for (int i=0;i<nDvs;i++) if (ppDvs[i]!=0) delete ppDvs[i];
        delete[] ppDvs;
    }
    nDvs=0;
    if (debug) cout<<"finished ModelDevsVectorVector::deallocate(void) "<<this<<endl;
}


/***************************************************************
*   Gets the integer index corresponding to the input name.    *
*   Returns 0 if no match.                                     *
***************************************************************/
int ModelDevsVectorVector::getIndexForName(adstring & name){
    if (debug) cout<<"starting ModelDevsVectorVector::getIndexForName(adstring & name) "<<this<<endl;
    if (!ppDvs) {
        if (debug) cout<<"finished ModelDevsVectorVector::getIndexForName(adstring & name) "<<this<<tb<<0<<endl;
        return 0;
    }
    for (int i=1;i<=nDvs;i++) {
        if (name==ppDvs[i-1]->getName()) {
            if (debug) cout<<"finished ModelDevsVectorVector::getIndexForName(adstring & name) "<<this<<tb<<i<<endl;
            return i;
        }
    }
    if (debug) cout<<"finished ModelDevsVectorVector::getIndexForName(adstring & name) "<<this<<tb<<0<<endl;
    return 0;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector ModelDevsVectorVector::getPhases(void){
    if (debug) cout<<"starting ModelDevsVectorVector::getPhases(void) "<<this<<endl;
    ivector phases(1,nDvs);
    phases.initialize();
    if (ppDvs!=0) for (int i=1;i<=nDvs;i++) phases(i) = ppDvs[i-1]->getPhase();
    if (debug) cout<<"finished ModelDevsVectorVector::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector ModelDevsVectorVector::getPriorWgts(void){
    if (debug) cout<<"starting ModelDevsVectorVector::getPriorWgts()"<<this<<endl;
    dvector wts(1,nDvs);
    for (int i=1;i<=nDvs;i++) wts(i) = ppDvs[i-1]->getPriorWgt();
    if (debug) cout<<"finished ModelDevsVectorVector::getPriorWgts()"<<this<<endl;
    return wts;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector ModelDevsVectorVector::getMinIndices(void){
    if (debug) cout<<"starting ModelDevsVectorVector::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nDvs);
    idxs.initialize();
    if (ppDvs!=0) for (int i=1;i<=nDvs;i++) idxs(i) = ppDvs[i-1]->getMinIndex();
    if (debug) cout<<"finished ModelDevsVectorVector::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector ModelDevsVectorVector::getMaxIndices(void){
    if (debug) cout<<"starting ModelDevsVectorVector::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nDvs);
    idxs.initialize();
    if (ppDvs!=0) for (int i=1;i<=nDvs;i++) idxs(i) = ppDvs[i-1]->getMaxIndex();
    if (debug) cout<<"finished ModelDevsVectorVector::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

// /***************************************************************
// *   get initial parameter values (on parameter scale).         *
// ***************************************************************/
// dvector ModelDevsVectorVector::getInitParamVals(int i){
//     if (debug) cout<<"starting ModelDevsVectorVector::getInitParamVals(i)"<<this<<endl;
//     if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]->getInitParamVals();
//     dvector initVals;
//     if (debug) cout<<"finished ModelDevsVectorVector::getInitParamVals(i)"<<this<<endl;
//     return initVals;
// }

/***************************************************************
*   get initial parameter values (on parameter scale).         *
***************************************************************/
dmatrix ModelDevsVectorVector::getInitParamVals(){
    if (debug) cout<<"starting ModelDevsVectorVector::getInitParamVals()"<<this<<endl;
    dmatrix initVals(1,nDvs,getMinIndices(),getMaxIndices());
    if (ppDvs) for (int i=1;i<=nDvs;i++) initVals(i) = ppDvs[i-1]->getInitParamVals();
    if (debug) cout<<"finished ModelDevsVectorVector::getInitParamVals()"<<this<<endl;
    return initVals;
}

/***************************************************************
*   get initial parameter values (on parameter scale).         *
***************************************************************/
dvector ModelDevsVectorVector::drawInitParamVals(int i,random_number_generator& rng){
    if (debug) cout<<"starting ModelDevsVectorVector::drawInitParamVals(rng,i)"<<this<<endl;
    if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]->drawInitParamVals(rng);
    if (debug) cout<<"finished ModelDevsVectorVector::drawInitParamVals(rng,i)"<<this<<endl;
    dvector initVals;
    return initVals;
}

/***************************************************************
*   get initial parameter values (on parameter scale).         *
***************************************************************/
dmatrix ModelDevsVectorVector::drawInitParamVals(random_number_generator& rng){
    if (debug) cout<<"starting ModelDevsVectorVector::drawInitParamVals(random_number_generator& rng)"<<this<<endl;
    dmatrix initVals(1,nDvs,getMinIndices(),getMaxIndices());
    for (int i=1;i<=nDvs;i++) initVals(i) = ppDvs[i-1]->drawInitParamVals(rng);
    if (debug) cout<<"finished ModelDevsVectorVector::drawInitParamVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   set param_init_vector_vector                               *
***************************************************************/
void ModelDevsVectorVector::setDevParamsVectorVector(param_init_vector_vector& pv){
    if (debug) cout<<"starting ModelDevsVectorVector::setDevParamsVectorVector(param_init_vector_vector& pv)"<<this<<endl;
    if ((pv.indexmin()==1)&&(pv.indexmax()==nDvs)) {
        for (int i=1;i<=nDvs;i++) ppDvs[i-1]->setDevParamsVector(pv(i));//set ptr to devs vector
    } else {
        cout<<"Incompatible sizes in ModelDevsVectorVector::setDevParamsVectorVector(pv)"<<endl;
        cout<<"Valid row indices are "<<1<<" to "<<nDvs<<endl;
        cout<<"pv has row indices "<<pv.indexmin()<<" to "<<pv.indexmax()<<endl;
        cout<<"Terminating..."<<endl;
        exit(-1);
    }
    if (debug) cout<<"finished ModelDevsVectorVector::setDevParamsVectorVector(param_init_vector_vector& pv)"<<this<<endl;
}

/***************************************************************
*   Check that the given indices correspond to a dev value.    *
*   The index 'i' refers to the ith devs vector in the         *
*   underlying param_init_vector_vector.                       *
***************************************************************/
int ModelDevsVectorVector::checkIndices(int i,int j){
    if (debug) cout<<"starting ModelDevsVectorVector::checkIndices(int i, int j) "<<this<<tb<<i<<tb<<j<<endl;
    if ((ppDvs)&&(i<=nDvs)) {
        return ppDvs[i-1]->checkIndex(j);
    }
    if (debug) cout<<"done ModelDevsVectorVector::checkIndices(int i, int j) "<<this<<endl;
    return 0;
}

/***************************************************************
*   Calculate log prior probability for ith devs vector.       *
*   Note that priors are specified and calculated on the       *
*   natural scale.                                             *
*   Must call calcNaturalScaleDevs() first.                    *
***************************************************************/
dvar_vector ModelDevsVectorVector::calcLogPriors(int i){
    if (debug) cout<<"starting ModelDevsVectorVector::calcLogPriors()"<<this<<endl;
    if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]->calcLogPriors();
    if (debug) cout<<"finished ModelDevsVectorVector::calcLogPriors()"<<this<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector lps;
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

// /***************************************************************
// *   Calculate log prior probability for all devs vectors.      *
// *   Note that priors are specified and calculated on the       *
// *   natural scale.                                             *
// *   Must call calcNaturalScaleDevs() first.                    *
// ***************************************************************/
// dvar_matrix ModelDevsVectorVector::calcLogPriors(){
//     RETURN_ARRAYS_INCREMENT();
//     if (debug) cout<<"starting ModelDevsVectorVector::calcLogPriors()"<<this<<endl;
//     ivector mnIndxs(1,nDvs);
//     ivector mxIndxs(1,nDvs);
//     for (int i=1;i<=nDvs;i++) {
//         mnIndxs(i) = ppDvs[i-1]->pDVI->getMinIndex();
//         mxIndxs(i) = ppDvs[i-1]->pDVI->getMaxIndex();
//     }
//     dvar_matrix lps(1,nDvs,mnIndxs,mxIndxs);
//     for (int i=1;i<=nDvs;i++) lps(i) = ppDvs[i-1]->calcLogPriors();
//     if (debug) cout<<"finished ModelDevsVectorVector::calcLogPriors()"<<this<<endl;
//     RETURN_ARRAYS_DECREMENT();
//     return lps;
// }

/***************************************************************
*   input operator                                            *
***************************************************************/
void ModelDevsVectorVector::read(cifstream & is){
    if (debug) cout<<"starting ModelDevsVectorVector::read(cifstream & is) "<<this<<endl;
    is>>nDvs;
    if (ppDvs) deallocate();
    if (nDvs>0) {
        ppDvs = new ModelDevsVector*[nDvs];
        idxNames.allocate(1,nDvs);
        for (int i=0;i<nDvs;i++) {
            is>>idxNames(i+1);
            ppDvs[i] = new ModelDevsVector(idxNames(i+1));
            is>>(*ppDvs[i]);
            if (debug)  cout<<(*ppDvs[i]);
        }
    }
    if (ModelDevsVectorVector::debug) cout<<"finished ModelDevsVectorVector::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   output operator                                            *
***************************************************************/
void ModelDevsVectorVector::write(ostream & os){
    os<<name<<tb<<"#model devs parameter vector vector name"<<endl;
    os<<tb<<nDvs<<tb<<"#number of dev vectors"<<endl;
    os<<"#index  phs   priorWgt   minYr   maxYr   initVals?  transform   priorType   priorParams"<<endl;
    for (int i=0;i<nDvs;i++) os<<(*ppDvs[i])<<endl;
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVectorVector::writeInfoToR(ostream& os, adstring nm, int indent){
    if (ppDvs) {
        os<<nm<<"=list("<<endl;
        for (int i=1;i<nDvs;i++) {ppDvs[i-1]->writeInfoToR(os,"'"+idxNames(i)+"'",indent+1); os<<","<<endl;}
        int i=nDvs;               ppDvs[i-1]->writeInfoToR(os,"'"+idxNames(i)+"'",indent+1); os<<endl;
        os<<")";
    } else {
        os<<nm<<"=NULL";
    }
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVectorVector::writeValuesToR(ostream& os, adstring nm, int indent){
    if (ppDvs) {
        os<<nm<<"=list("<<endl;
        for (int i=1;i<nDvs;i++) {ppDvs[i-1]->writeValuesToR(os,"'"+idxNames(i)+"'",indent+1); os<<","<<endl;}
        int i=nDvs;               ppDvs[i-1]->writeValuesToR(os,"'"+idxNames(i)+"'",indent+1); os<<endl;
        os<<")";
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////ModelDevsVectorVector/////////////////////////////////

//--------------------------------------------------------------------------------
//          ModelDevsVectorMatrix
//--------------------------------------------------------------------------------
/***************************************************************
*   deallocation                                               *
***************************************************************/
void ModelDevsVectorMatrix::deallocate(void){
    if (debug) cout<<"starting ModelDevsVectorMatrix::deallocate(void) "<<this<<endl;
    if (ppDVVs!=0) {
        for (int i=0;i<nDVVs;i++) if (ppDVVs[i]!=0) delete ppDVVs[i];
        delete[] ppDVVs;
    }
    nDVVs = 0;
    nDvs  = 0;
    if (debug) cout<<"finished ModelDevsVectorMatrix::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   Gets the integer index corresponding to the input row name.*
*   Returns 0 if no match.                                     *
***************************************************************/
int ModelDevsVectorMatrix::getRowIndexForName(adstring & rowName){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getRowIndexForName(adstring & rowName) "<<this<<tb<<rowName<<endl;
    if (!ppDVVs) {
        if (debug) cout<<"finished ModelDevsVectorMatrix::getRowIndexForName(adstring & rowName) "<<this<<tb<<0<<endl;
        return 0;
    }
    for (int r=1;r<=nDVVs;r++) {
        if (rowName==rowNames(r)) {
            if (debug) cout<<"finished ModelDevsVectorMatrix::getRowIndexForName(adstring & rowName) "<<this<<tb<<r<<endl;
            return r;
        }
    }
    if (debug) {
        if (debug) cout<<"finished ModelDevsVectorMatrix::getRowIndexForName(adstring & rowName) "<<this<<tb<<0<<endl;
    }
    return 0;
}

/***************************************************************
*   Gets the integer index corresponding to the input row and  *
*   column names.  Returns 0 if no match.                      *
***************************************************************/
int ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName) "<<this<<tb<<rowName<<tb<<colName<<endl;
    if (!ppDVVs) {
        if (debug) cout<<"finished ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName) "<<this<<tb<<0<<endl;
        return 0;
    }
    int idx = 0;
    for (int r=1;r<=nDVVs;r++){
        if (!(rowName==rowNames(r))){
            idx += ppDVVs[r-1]->getNP();
        } else {
            int c = ppDVVs[r-1]->getIndexForName(colName);
            if (c>0) {
                if (debug) cout<<"finished ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName) "<<this<<tb<<idx+c<<endl;
                return idx+c;
            } else {
                if (debug) {
                    cout<<"finished ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName) "<<this<<tb<<0<<endl;
                }
                return 0;
            }
        }
    }
    if (debug) {
        cout<<"finished ModelDevsVectorMatrix::getIndexForNames(adstring & rowName, adstring & colName) "<<this<<tb<<0<<endl;
    }
    return 0;//shouldn't ever get here
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector ModelDevsVectorMatrix::getPhases(void){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getPhases(void) "<<this<<endl;
    ivector phases;
    if (ppDVVs) {
        phases.allocate(1,nDvs);
        int i = 1;
        for (int j=1;j<=nDVVs;j++) {
            ModelDevsVectorVector* pDVV = ppDVVs[j-1];
            ivector phases1(pDVV->getPhases());
            for (int k=1;k<=pDVV->getNP();k++) phases(i++) = phases1(k);
        }
    }
    if (debug) {
        if (phases.allocated()) {cout<<phases<<endl;} else {cout<<"unallocated"<<endl;}
        cout<<"finished ModelParameterMatrix::getPhases(void) "<<this<<endl;
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector ModelDevsVectorMatrix::getPriorWgts(void){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getPriorWgts()"<<this<<endl;
    dvector wts;
    if (ppDVVs) {
        wts.allocate(1,nDvs);
        int i = 1;
        for (int j=1;j<=nDVVs;j++) {
            ModelDevsVectorVector* pDVV = ppDVVs[j-1];
            dvector wts1(pDVV->getPriorWgts());
            for (int k=1;k<=pDVV->getNP();k++) wts(i++) = wts1(k);
        }
    }
    if (debug) {
        if (wts.allocated()) {cout<<wts<<endl;} else {cout<<"unallocated"<<endl;}
        cout<<"finished ModelDevsVectorMatrix::getPriorWgts()"<<this<<endl;
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getPriorWgts()"<<this<<endl;
    return wts;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector ModelDevsVectorMatrix::getMinIndices(void){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getMinIndices() "<<this<<endl;
    dvector idxs;
    if (ppDVVs) {
        idxs.allocate(1,nDvs);
        int i = 1;
        for (int j=1;j<=nDVVs;j++) {
            ModelDevsVectorVector* pDVV = ppDVVs[j-1];
            dvector idxs1(pDVV->getMinIndices());
            for (int k=1;k<=pDVV->getNP();k++) idxs(i++) = idxs1(k);
        }
    }
    if (debug) {
        if (idxs.allocated()) {cout<<idxs<<endl;} else {cout<<"unallocated"<<endl;}
        cout<<"finished ModelDevsVectorMatrix::getMinIndices()"<<this<<endl;
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getMinIndices() "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector ModelDevsVectorMatrix::getMaxIndices(void){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getMaxIndices() "<<this<<endl;
    dvector idxs;
    if (ppDVVs) {
        idxs.allocate(1,nDvs);
        int i = 1;
        for (int j=1;j<=nDVVs;j++) {
            ModelDevsVectorVector* pDVV = ppDVVs[j-1];
            dvector idxs1(pDVV->getMaxIndices());
            for (int k=1;k<=pDVV->getNP();k++) idxs(i++) = idxs1(k);
        }
    }
    if (debug) {
        if (idxs.allocated()) {cout<<idxs<<endl;} else {cout<<"unallocated"<<endl;}
        cout<<"finished ModelDevsVectorMatrix::getMaxIndices()"<<this<<endl;
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getMaxIndices() "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get initial parameter values (on parameter scale).         *
***************************************************************/
dmatrix ModelDevsVectorMatrix::getInitParamVals(){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getInitParamVals()"<<this<<endl;
    dmatrix initVals;
    if (ppDVVs) {
        initVals.allocate(1,nDvs,getMinIndices(),getMaxIndices());
        int j = 1;
        for (int i=1;i<=nDVVs;i++) {
            ModelDevsVectorVector* pDVV = ppDVVs[i-1];
            dmatrix initVals1(pDVV->getInitParamVals());
            for (int k=1;k<=pDVV->getNP();k++) initVals(j++) = initVals1(k);
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getInitParamVals()"<<this<<endl;
    return initVals;
}

// /***************************************************************
// *   get initial parameter values (on parameter scale).         *
// ***************************************************************/
// dvector ModelDevsVectorMatrix::drawInitParamVals(int i,random_number_generator& rng){
//     if (debug) cout<<"starting ModelDevsVectorMatrix::drawInitParamVals(rng,i)"<<this<<endl;
//     if ((ppDVVs)&&(i<=nDVVs)) return ppDVVs[i-1]->drawInitParamVals(rng);
//     if (debug) cout<<"finished ModelDevsVectorMatrix::drawInitParamVals(rng,i)"<<this<<endl;
//     dvector initVals;
//     return initVals;
// }

/***************************************************************
*   get initial parameter values (on parameter scale).         *
***************************************************************/
dmatrix ModelDevsVectorMatrix::drawInitParamVals(random_number_generator& rng){
    if (debug) cout<<"starting ModelDevsVectorMatrix::drawInitParamVals(random_number_generator& rng)"<<this<<endl;
    dmatrix initVals;
    if (ppDVVs) {
        initVals.allocate(1,nDvs,getMinIndices(),getMaxIndices());
        int j = 1;
        for (int i=1;i<=nDVVs;i++) {
            ModelDevsVectorVector* pDVV = ppDVVs[i-1];
            dmatrix initVals1(pDVV->drawInitParamVals(rng));
            for (int k=1;k<=pDVV->getNP();k++) initVals(j++) = initVals1(k);
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::drawInitParamVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   set param_init_vector_vector                               *
***************************************************************/
void ModelDevsVectorMatrix::setDevParamsVectorVector(param_init_vector_vector& pv){
    if (debug) cout<<"starting ModelDevsVectorMatrix::setDevParamsVectorVector(param_init_vector_vector& pv)"<<this<<endl;
    if ((pv.indexmin()==1)&&(pv.indexmax()==nDvs)) {
        int j=1;
        for (int i=1;i<=nDVVs;i++) {
            ModelDevsVectorVector* pDVV = ppDVVs[i-1];
            for (int k=1;k<=pDVV->getNP();k++) pDVV->ppDvs[k-1]->setDevParamsVector(pv(j++));//set ptr to devs vector
        }
    } else {
        cout<<"Incompatible sizes in ModelDevsVectorMatrix::setDevParamsVectorVector(pv)"<<endl;
        cout<<"Valid row indices are "<<1<<" to "<<nDvs<<endl;
        cout<<"pv has row indices "<<pv.indexmin()<<" to "<<pv.indexmax()<<endl;
        cout<<"Terminating..."<<endl;
        exit(-1);
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::setDevParamsVectorVector(param_init_vector_vector& pv)"<<this<<endl;
}

/***************************************************************
*   Check that the given indices correspond to a dev value.    *
*   The index 'i' refers to the ith devs vector in the         *
*   underlying param_init_vector_vector, NOT into the          *
*   ppDVVS array.                                              *
***************************************************************/
int ModelDevsVectorMatrix::checkIndices(int i,int j){
    if (debug) cout<<"starting ModelDevsVectorMatrix::checkIndices(int i, int j) "<<this<<tb<<i<<tb<<j<<endl;
    if ((ppDVVs)&&(0<i)&&(i<=nDvs)) {
        int m=0;
        for (int k=0;k<nDVVs;k++) {
            int dm = ppDVVs[k]->getNP();
            if (i<=(m+dm)) {
                return ppDVVs[k]->checkIndices(i-m,j);
            }
            m += dm;
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::checkIndices(int i, int j) "<<this<<endl;
    return 0; //false
}

/***************************************************************
*   Get the jth dev from the ith devs vector.                  *
*   The index 'i' refers to the ith devs vector in the         *
*   underlying param_init_vector_vector, NOT into the          *
*   ppDVVS array.                                              *
*   Note that priors are specified and calculated on the       *
*   natural scale.                                             *
*   Must call calcNaturalScaleDevs() first.                    *
***************************************************************/
dvariable ModelDevsVectorMatrix::getNaturalScaleDev(int i,int j){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getNaturalScaleDevsVector(int i, int j) "<<this<<endl;
    if ((ppDVVs)&&(0<i)&&(i<=nDvs)) {
        int m=0;
        for (int k=0;k<nDVVs;k++) {
            int dm = ppDVVs[k]->getNP();
            if (i<=(m+dm)) {
                return ppDVVs[k]->getNaturalScaleDev(i-m,j);
            }
            m += dm;
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getNaturalScaleDevsVector(int i, int j) "<<this<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable val;
    RETURN_ARRAYS_DECREMENT();
    return val;
}

/***************************************************************
*   Get the jth dev from the ith devs vector.                  *
*   The index 'i' refers to the ith devs vector in the         *
*   underlying param_init_vector_vector, NOT into the          *
*   ppDVVS array.                                              *
*   Note that priors are specified and calculated on the       *
*   natural scale.                                             *
*   Must call calcNaturalScaleDevs() first.                    *
***************************************************************/
dvar_vector ModelDevsVectorMatrix::getNaturalScaleDevsVector(int i){
    if (debug) cout<<"starting ModelDevsVectorMatrix::getNaturalScaleDevsVector(int i) "<<this<<endl;
    if ((ppDVVs)&&(0<i)&&(i<=nDvs)) {
        int m=0;
        for (int k=0;k<nDVVs;k++) {
            int dm = ppDVVs[k]->getNP();
            if (i<=(m+dm)) {
                return ppDVVs[k]->getNaturalScaleDevsVector(i-m);
            }
            m += dm;
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::getNaturalScaleDevsVector(int i) "<<this<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector vec;
    RETURN_ARRAYS_DECREMENT();
    return vec;
}

/***************************************************************
*   Calculate log prior probability for ith devs vector.       *
*   The index 'i' refers to the ith devs vector in the         *
*   underlying param_init_vector_vector, NOT into the          *
*   ppDVVS array.                                              *
*   Note that priors are specified and calculated on the       *
*   natural scale.                                             *
*   Must call calcNaturalScaleDevs() first.                    *
***************************************************************/
dvar_vector ModelDevsVectorMatrix::calcLogPriors(int i){
    if (debug) cout<<"starting ModelDevsVectorMatrix::calcLogPriors(int i)"<<this<<endl;
    if ((ppDVVs)&&(0<i)&&(i<=nDvs)) {
        int j=0;
        for (int k=0;k<nDVVs;k++) {
            int dj = ppDVVs[k]->getNP();
            if (i<=(j+dj)) {
//                if (debug) cout<<"i,j,k = "<<i<<tb<<j<<tb<<k<<endl;
                return ppDVVs[k]->calcLogPriors(i-j);
            }
            j += dj;
        }
    }
    if (debug) cout<<"finished ModelDevsVectorMatrix::calcLogPriors(int i)"<<this<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector lps;
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

// /***************************************************************
// *   Calculate log prior probability for all devs vectors.      *
// *   Note that priors are specified and calculated on the       *
// *   natural scale.                                             *
// *   Must call calcNaturalScaleDevs() first.                    *
// ***************************************************************/
// dvar_matrix ModelDevsVectorMatrix::calcLogPriors(){
//     RETURN_ARRAYS_INCREMENT();
//     if (debug) cout<<"starting ModelDevsVectorMatrix::calcLogPriors()"<<this<<endl;
//     ivector mnIndxs(1,nDVVs);
//     ivector mxIndxs(1,nDVVs);
//     for (int i=1;i<=nDVVs;i++) {
//         mnIndxs(i) = ppDVVs[i-1]->pDVI->getMinIndex();
//         mxIndxs(i) = ppDVVs[i-1]->pDVI->getMaxIndex();
//     }
//     dvar_matrix lps(1,nDVVs,mnIndxs,mxIndxs);
//     for (int i=1;i<=nDVVs;i++) lps(i) = ppDVVs[i-1]->calcLogPriors();
//     if (debug) cout<<"finished ModelDevsVectorMatrix::calcLogPriors()"<<this<<endl;
//     RETURN_ARRAYS_DECREMENT();
//     return lps;
// }

/***************************************************************
*   input operator                                            *
***************************************************************/
void ModelDevsVectorMatrix::read(cifstream & is){
    if (debug) cout<<"starting ModelDevsVectorMatrix::read(cifstream & is) "<<this<<endl;
    is>>nDVVs;
    if (ppDVVs) deallocate();
    if (nDVVs>0) {
        ppDVVs = new ModelDevsVectorVector*[nDVVs];
        rowNames.allocate(1,nDVVs);
        for (int i=0;i<nDVVs;i++) {
            is>>rowNames(i+1);
            ppDVVs[i] = new ModelDevsVectorVector(rowNames(i+1));
            is>>(*ppDVVs[i]);
            if (debug)  cout<<(*ppDVVs[i]);
            nDvs += ppDVVs[i]->getNP();
        }
    }
    if (ModelDevsVectorMatrix::debug) cout<<"finished ModelDevsVectorMatrix::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   output operator                                            *
***************************************************************/
void ModelDevsVectorMatrix::write(ostream & os){
    os<<name<<tb<<"#model devs vector matrix name"<<endl;
    os<<tb<<nDVVs<<tb<<"#number of dev vector vectors"<<endl;
    for (int i=0;i<nDVVs;i++) os<<(*ppDVVs[i])<<endl;
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVectorMatrix::writeInfoToR(ostream& os, adstring nm, int indent){
    if (ppDVVs){
        os<<nm<<"=list("<<endl;
        for (int i=1;i<nDVVs;i++) {ppDVVs[i-1]->writeInfoToR(os,"'"+rowNames(i)+"'",indent+1); os<<","<<endl;}
        int i=nDVVs;               ppDVVs[i-1]->writeInfoToR(os,"'"+rowNames(i)+"'",indent+1); os<<endl;
        os<<")";
    } else {
        os<<nm<<"=NULL";
    }
}

/***************************************************************
*   Write to stream in R format.                               *
***************************************************************/
void ModelDevsVectorMatrix::writeValuesToR(ostream& os, adstring nm, int indent){
    if (ppDVVs){
        os<<nm<<"=list("<<endl;
        for (int i=1;i<nDVVs;i++) {ppDVVs[i-1]->writeValuesToR(os,"'"+rowNames(i)+"'",indent+1); os<<","<<endl;}
        int i=nDVVs;               ppDVVs[i-1]->writeValuesToR(os,"'"+rowNames(i)+"'",indent+1); os<<endl;
        os<<")";
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////ModelDevsVectorMatrix/////////////////////////////////
