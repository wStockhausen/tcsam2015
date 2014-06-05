//ModelParameterFunctions.cpp

#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelParameterInfoTypes.hpp"
#include "ModelParameterFunctions.hpp"



using namespace tcsam;

//******************************************************************************
//* Function: void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& cout)
//* 
//* Description: Sets initial values for a parameter vector.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (BoundedNumberVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_number_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();
        for (int i=1;i<=np;i++) p(i)=vls(i);
        rpt::echo<<"InitVals for "<<p(1).get_name()<<": "<<p<<std::endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

//******************************************************************************
//* Function: void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout)
//* 
//* Description: Sets initial values for a vector of parameter vectors.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (BoundedVectorVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug>=tcsam::dbgAll) cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<std::endl;
            for (int j=vls.indexmin();j<=vls.indexmax();j++) p(i,j)=vls(j);
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals "<<p(i).get_name()<<":"<<tb<<p(i)<<std::endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

//******************************************************************************
//* Function: void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout)
//* 
//* Description: Sets initial values for a vector of parameter vectors.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (DevsVectorVectorInfo*) 
//*     pointer to DevsNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug) cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<std::endl;
            for (int j=vls.indexmin();j<=(vls.indexmax()-1);j++) p(i,j)=vls(j);
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals "<<p(i).get_name()<<":"<<tb<<p(i)<<std::endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

//-------------------------------------------------------------------------------------
void setDevs(dvar_matrix& devs, param_init_bounded_vector_vector& pDevs, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"starting setDevs(devs,pDevs)"<<std::endl;
    int nv = pDevs.indexmax();//number of devs vectors defined
    int mni; int mxi;
    for (int v=1;v<=nv;v++){
        mni = pDevs(v).indexmin();
        mxi = pDevs(v).indexmax();
        devs(v)(mni,mxi) = pDevs(v);
        devs(v,mxi+1) = -sum(devs(v)(mni,mxi));
    }
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>nv;
        if (nv<0) exit(-1);
        cout<<"finished setDevs(devs,pDevs)"<<std::endl;
    }
}

//-------------------------------------------------------------------------------------
//Calculate priors                                      
void calcPriors(objective_function_value& objFun, NumberVectorInfo* ptrVI,param_init_number_vector& pv, int debug, std::ostream& cout){
    if (ptrVI->getSize()){
        dvar_vector tmp = pv;
        dvar_vector pri = ptrVI->calcLogPriors(tmp);//ln-scale prior (NOT NLL!)
        if (debug>=dbgPriors) cout<<"priors("<<ptrVI->name<<") = "<<pri<<std::endl;
        objFun += -ptrVI->getPriorWgts()*pri;
    }
}

//-------------------------------------------------------------------------------------
//Calculate priors                                      
void calcPriors(objective_function_value& objFun, BoundedNumberVectorInfo* ptrVI,param_init_bounded_number_vector& pv, int debug, std::ostream& cout){
    if (ptrVI->getSize()){
        dvar_vector tmp = pv;
        dvar_vector pri = ptrVI->calcLogPriors(tmp);//ln-scale prior (NOT NLL!)
        if (debug>=dbgPriors) cout<<"priors("<<ptrVI->name<<") = "<<pri<<std::endl;
        objFun += -ptrVI->getPriorWgts()*pri;
    }
}

//-------------------------------------------------------------------------------------
//Calculate contributions to objective function from priors                                      
void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, param_init_bounded_vector_vector& pm, int debug, std::ostream& cout){
    if (ptrVVI->getSize()){
        dvector wts = ptrVVI->getPriorWgts();
        for (int i=1;i<=pm.indexmax();i++) {
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (debug>=dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
        }
    }    
}

//-------------------------------------------------------------------------------------
//Calculate contributions to objective function from priors                                      
void calcPriors(objective_function_value& objFun, DevsVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout){
    if (ptrVVI->getSize()){
        dvector wts = ptrVVI->getPriorWgts();
        for (int i=1;i<=pm.indexmax();i++) {
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (debug>=dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
        }
    }    
}