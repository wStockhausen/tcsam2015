//--------------------------------------------------------------------------------
//  Includes:
//      DevsVectorInfo
//      ModelDevsVector
//      ModelDevsVectorVector
//      ModelDevsVectorMatrix
//--------------------------------------------------------------------------------
#pragma once
#ifndef MODELDEVS_HPP
    #define MODELDEVS_HPP

#include <admodel.h>

    //forward declarations
    class ModelDevsVectorMatrix;
    class ModelPDFInfo;
    class ModelTransformInfo;

//--------------------------------------------------------------------------------
//          DevsVectorInfo
//  Defines characteristics for a single devs vector.
//  Vector indices run from mnIndex to mxIndex.
//  Values are such that val(mxIndx) = -sum(val(mnIndx:mxIndx-1) so that
//  the total vector has a mean of 0.  As a consequence, the final value
//  in initVals is recalculated to ensure this after reading from input.
//--------------------------------------------------------------------------------
    class DevsVectorInfo {
        public:
            static int debug;//flag to print debugging info
        protected:
            int mnIndx; //min index for devs vector
            int mxIndx; //max index for devs vector
            int flgResample;     //flag to resample values
            int flgInitVals;     //flag to read init vals
            int phase;           //phase to turn on devs
            double priorWgt;     //likelihood weight for prior
            dvector initVals;    //initial values on the "natural" scale
            dvector priorParams; //parameters for priors
            dvector priorConsts; //constants for priors
        public:
            adstring name;             //external id
            ModelPDFInfo* pMPI;        //ptr to prior info
            ModelTransformInfo* pMTI;  //ptr to transform info
        public:
            DevsVectorInfo(adstring & name){this->name=name;pMPI=0;pMTI=0;}
            ~DevsVectorInfo();
            int getPhase(){return phase;}
            double getPriorWgt(){return priorWgt;}
            int getMinIndex(){return mnIndx;}
            int getMaxIndex(){return mxIndx;}
            void setPriorType(adstring prior);
            void setTransformType(adstring type);
            dvector getInitParamVals();
            dvector drawInitParamVals(random_number_generator& rng);
            void updateInitVals(dvector vals){initVals=vals;}
            dvar_vector calcLogPriors(dvar_vector& xv);
            void read(cifstream & is);
            void write(ostream & os);
            void writeToR(ostream& os);
            friend cifstream& operator >>(cifstream & is, DevsVectorInfo & obj){obj.read(is); return is;}
            friend ostream& operator <<(ostream & os, DevsVectorInfo & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          ModelDevsVector
//  Encapsulates a param_init_vector that functions as a devs vector
//--------------------------------------------------------------------------------
    class ModelDevsVector {
        public:
            static int debug;//flag to print debugging info
        public:
            DevsVectorInfo* pDVI;   //ptr to dev vector info
            dvar_vector* pASDvs;    //natural-scale devs
            param_init_vector* pDvs;//ptr to param_init_vector
        public:
            ModelDevsVector(adstring & name) {pDVI=new DevsVectorInfo(name); pASDvs=0;pDvs=0;}
            ~ModelDevsVector(){delete pDVI; if (pASDvs) delete pASDvs;pASDvs=0;pDvs=0;}
            adstring getName(){return pDVI->name;}
            int isActive(){if (pDvs) return active((*pDvs)); return 0;}
            int getPhase(){return pDVI->getPhase();}     //get phase for devs vector
            double getPriorWgt(){return pDVI->getPriorWgt();}
            int getMinIndex(){return pDVI->getMinIndex();}//get min index for devs vector
            int getMaxIndex(){return pDVI->getMaxIndex();}//get max index for devs vector
            dvector getInitParamVals(){return pDVI->getInitParamVals();}//returned vector has same indices as dvi.initVals
            dvector drawInitParamVals(random_number_generator& rng){return pDVI->drawInitParamVals(rng);}//returned vector has same indices as dvi.initVals
            void updateInitVals(){if (pASDvs) pDVI->updateInitVals(value(*pASDvs));}
            //call once in PRELIMINARY_CALCS_SECTION
            void setDevParamsVector(param_init_vector& pars);//set the parameter vector for the devs
            //call once in PROCEDURE_SECTION (can also call in PRELIMINARY_CALCS)
            void   calcNaturalScaleDevs(void);//calculate the natural-scale devs
            int checkIndex(int i);//check that index into dev vector is ok (1=true, 0=false)
            dvariable getNaturalScaleDev(int i);//get ith dev; call calcNaturalScaleDevs() first
            dvar_vector getNaturalScaleDevsVector();//call calcNaturalScaleDevs() first
            dvar_vector calcLogPriors();//call calcNaturalScaleDevs() first
            void writeInfoToR(ostream& os, adstring nm, int indent=0);
            void writeValuesToR(ostream& os, adstring nm, int indent=0);
            friend cifstream& operator >>(cifstream & is, ModelDevsVector & obj){is>>(*obj.pDVI);return is;}
            friend ostream& operator <<(ostream & os, ModelDevsVector & obj){os<<(*obj.pDVI);return os;}
        protected:
            void createNaturalScaleDevs(void);
    };

//--------------------------------------------------------------------------------
//          ModelDevsVectorVector
//  Encapsulates a param_init_vector_vector
//  that functions as a 1-d array of devs vectors
//--------------------------------------------------------------------------------
    class ModelDevsVectorVector {
        public:
            static int debug;//flag to print debugging info
        public:
            adstring name;
            adstring_array idxNames;
            int nDvs;               //number of dev vectors
            ModelDevsVector** ppDvs;//ptr to vector of ptrs to ModelDevsVector instances
        public:
            ModelDevsVectorVector(adstring & name){this->name=name;nDvs=0;ppDvs=0;}
            ~ModelDevsVectorVector(){deallocate();}
            ModelDevsVector* operator[](int i){if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]; else return 0;}
            void deallocate();
            int getNP(){return nDvs;}
            int getIndexForName(adstring & name);
            int isActive(int i){if (ppDvs&&(i<=nDvs)) return (ppDvs[i-1])->isActive(); return 0;}
            int getPhase(int i){if (ppDvs&&(i<=nDvs)) return (ppDvs[i-1])->getPhase(); return -999;}
            ivector getPhases();            //get phases for all devs vectors
            dvector getPriorWgts();         //get prior weights for all devs vectors
            ivector getMinIndices();        //get min indices for all devs vectors
            ivector getMaxIndices();        //get max indices for all devs vectors
            dmatrix getInitParamVals();     //get initial values for all devs vectors
            dvector drawInitParamVals(int i,random_number_generator& rng);//draw random sample for ith devs vector
            dmatrix drawInitParamVals(random_number_generator& rng);      //draw random sample for all devs vectors
            void    updateInitVals(){if (ppDvs) for (int i=0;i<nDvs;i++) ppDvs[i]->updateInitVals();}//update natural-scale initial values based on current values
//            void setNaturalScaleDevsMatrix(dvar_matrix& devs);
            //call once in PRELIMINARY_CALCS_SECTION
            void setDevParamsVectorVector(param_init_vector_vector & pv);
            //call once in PROCEDURE_SECTION
            void calcNaturalScaleDevs(void){if (ppDvs) for (int i=0;i<nDvs;i++) ppDvs[i]->calcNaturalScaleDevs();}
            //for following four functions, call calcNaturalScaleDevs() first
            int checkIndices(int i,int j);
            dvariable getNaturalScaleDev(int i,int j){if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]->getNaturalScaleDev(j);RETURN_ARRAYS_INCREMENT();dvariable res;RETURN_ARRAYS_DECREMENT();return res;}//call calcNaturalScaleDevs() first
            dvar_vector getNaturalScaleDevsVector(int i){if ((ppDvs)&&(i<=nDvs)) return ppDvs[i-1]->getNaturalScaleDevsVector();RETURN_ARRAYS_INCREMENT();dvar_vector res;RETURN_ARRAYS_DECREMENT();return res;}//call calcNaturalScaleDevs() first
            dvar_vector calcLogPriors(int i);
//            dvar_matrix calcLogPriors();
            void read(cifstream & is);
            void write(ostream & os);
            void writeInfoToR(ostream& os, adstring nm, int indent=0);
            void writeValuesToR(ostream& os, adstring nm, int indent=0);
            friend cifstream& operator >>(cifstream & is, ModelDevsVectorVector & obj){obj.read(is);return is;}
            friend ostream& operator <<(ostream & os, ModelDevsVectorVector & obj){obj.write(os);return os;}
            friend class ModelDevsVectorMatrix;
    };

//--------------------------------------------------------------------------------
//          ModelDevsVectorMatrix
//  Encapsulates a param_init_vector_vector
//  that functions as a 2-d array of devs vectors
//--------------------------------------------------------------------------------
    class ModelDevsVectorMatrix {
        public:
            static int debug;//flag to print debugging info
        public:
            adstring name;
            adstring_array rowNames;
            int nDvs;                      //total number of dev vectors
            int nDVVs;                     //number of dev vector vectors
            ModelDevsVectorVector** ppDVVs;//ptr to vector of ptrs to ModelDevsVectorVector instances
        public:
            ModelDevsVectorMatrix(adstring & name){this->name=name;nDvs=0;nDVVs=0;ppDVVs=0;}
            ~ModelDevsVectorMatrix(){deallocate();}
            ModelDevsVectorVector* operator[](int i){if ((ppDVVs)&&(i<=nDVVs)) return ppDVVs[i-1]; else return 0;}
            void deallocate();
            int getNP(){return nDvs;}
            int getNumRows(){return nDVVs;}
            int getRowIndexForName(adstring & name);
            int getIndexForNames(adstring & rowName, adstring & colName);
            ivector getPhases();            //get phases for all devs vectors
            dvector getPriorWgts();         //get prior weights for all devs vectors
            ivector getMinIndices();        //get min indices for all devs vectors
            ivector getMaxIndices();        //get max indices for all devs vectors
            dmatrix getInitParamVals();     //get initial values for all devs vectors
//             dvector drawInitParamVals(int i,random_number_generator& rng);//draw random sample for ith devs vector
            dmatrix drawInitParamVals(random_number_generator& rng);      //draw random sample for all devs vectors
            void    updateInitVals(){if (ppDVVs) for (int i=0;i<nDVVs;i++) ppDVVs[i]->updateInitVals();}//update natural-scale initial values based on current values
            //call once in PRELIMINARY_CALCS_SECTION
            void setDevParamsVectorVector(param_init_vector_vector & pv);
            //call once in PROCEDURE_SECTION
            void calcNaturalScaleDevs(void){if (ppDVVs) for (int i=0;i<nDVVs;i++) ppDVVs[i]->calcNaturalScaleDevs();}
            //for following four functions, call calcNaturalScaleDevs() first
            //Note that 'i' in the following is the row index into the param_init_vector_vector,
            //NOT the index into the array of ModelParameterVectorVector objects
            int checkIndices(int i,int j);
            dvariable getNaturalScaleDev(int i,int j);
            dvar_vector getNaturalScaleDevsVector(int i);
            dvar_vector calcLogPriors(int i);
//            dvar_matrix calcLogPriors();
            void read(cifstream & is);
            void write(ostream & os);
            void writeInfoToR(ostream& os, adstring nm, int indent=0);
            void writeValuesToR(ostream& os, adstring nm, int indent=0);
            friend cifstream& operator >>(cifstream & is, ModelDevsVectorMatrix & obj){obj.read(is);return is;}
            friend ostream& operator <<(ostream & os, ModelDevsVectorMatrix & obj){obj.write(os);return os;}
    };


#endif  //MODELDEVS_HPP
