/* 
 * File:   ModelParameterFunctions.hpp
 * Author: WilliamStockhausen
 *
 * Created on June 5, 2014, 9:30 AM
 */

#ifndef MODELPARAMETERFUNCTIONS_HPP
#define	MODELPARAMETERFUNCTIONS_HPP
namespace tcsam {
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
    void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& cout);

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
    void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout);

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
    void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout);

    //-------------------------------------------------------------------------------------
    void setDevs(dvar_matrix& devs, param_init_bounded_vector_vector& pDevs, int debug, std::ostream& cout);

    //-------------------------------------------------------------------------------------
    //Calculate priors                                      
    void calcPriors(objective_function_value& objFun, NumberVectorInfo* ptrVI,param_init_number_vector& pv, int debug, std::ostream& cout);

    //-------------------------------------------------------------------------------------
    //Calculate priors                                      
    void calcPriors(objective_function_value& objFun, BoundedNumberVectorInfo* ptrVI,param_init_bounded_number_vector& pv, int debug, std::ostream& cout);

    //-------------------------------------------------------------------------------------
    //Calculate contributions to objective function from priors                                      
    void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, param_init_bounded_vector_vector& pm, int debug, std::ostream& cout);

    //-------------------------------------------------------------------------------------
    //Calculate contributions to objective function from priors                                      
    void calcPriors(objective_function_value& objFun, DevsVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout);

} //namespace tcsam

#endif	/* MODELPARAMETERFUNCTIONS_HPP */

