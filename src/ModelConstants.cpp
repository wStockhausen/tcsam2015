#include <admodel.h>
#include "ModelConstants.hpp"

using namespace tcsam;

std::ofstream rpt::echo("EchoData.dat",std::ios::trunc);

int tcsam::getMaturityType(adstring s){
    s.to_upper();
    if (s==STR_IMMATURE)     return IMMATURE;     else
    if (s==STR_MATURE)       return MATURE;       else
    if (s==STR_ANY_MATURITY) return ANY_MATURITY; else
    std::cout<<"Unrecognized MaturityType '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getMaturityType(int i){
    if (i==IMMATURE)     return STR_IMMATURE;     else
    if (i==MATURE)       return STR_MATURE;       else
    if (i==ANY_MATURITY) return STR_ANY_MATURITY; else
    std::cout<<"Unrecognized MaturityType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getSexType(adstring s){
    s.to_upper();
    if (s==STR_FEMALE)  return FEMALE;  else
    if (s==STR_MALE)    return MALE;    else
    if (s==STR_ANY_SEX) return ANY_SEX; else
    std::cout<<"Unrecognized SexType '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getSexType(int i){
    if (i==FEMALE)  return STR_FEMALE;  else
    if (i==MALE)    return STR_MALE;    else
    if (i==ANY_SEX) return STR_ANY_SEX; else
    std::cout<<"Unrecognized SexType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getShellType(adstring s){
    s.to_upper();
    if (s==STR_NEW_SHELL) return NEW_SHELL; else
    if (s==STR_OLD_SHELL) return OLD_SHELL; else
    if (s==STR_ANY_SHELL) return ANY_SHELL; else
    std::cout<<"Unrecognized ShellType '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getShellType(int i){
    if (i==NEW_SHELL) return STR_NEW_SHELL; else
    if (i==OLD_SHELL) return STR_OLD_SHELL; else
    if (i==ANY_SHELL) return STR_ANY_SHELL; else
    std::cout<<"Unrecognized ShellType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getSRType(adstring s){
    s.to_lower();
    if (s==STR_CONSTANT) return SRTYPE_CONSTANT;
    if (s==STR_BEVHOLT)  return SRTYPE_BEVHOLT;
    if (s==STR_RICKER)   return SRTYPE_RICKER;
    std::cout<<"Unrecognized Stock-Recruit type (SRType) '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getSRType(int i){
    if (i==SRTYPE_CONSTANT) return STR_CONSTANT;
    if (i==SRTYPE_BEVHOLT)  return STR_BEVHOLT;
    if (i==SRTYPE_RICKER)   return STR_RICKER;
    std::cout<<"Unrecognized Stock-Recruit type (SRType) '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getScaleType(adstring s){
    if (s==STR_VAR) return SCLTYPE_VAR;
    if (s==STR_STD) return SCLTYPE_STD;
    if (s==STR_CV)  return SCLTYPE_CV;
    std::cout<<"Unrecognized scale type '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getScaleType(int i){
    if (i==SCLTYPE_VAR) return STR_VAR;
    if (i==SCLTYPE_STD) return STR_STD;
    if (i==SCLTYPE_CV)  return STR_CV;
    std::cout<<"Unrecognized scale type '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getFitType(adstring s){
    if (s==FIT_NONE_STR)     return FIT_NONE;
    if (s==FIT_BY_TOTAL_STR) return FIT_BY_TOTAL;
    if (s==FIT_BY_SEX_STR)   return FIT_BY_SEX;
    if (s==FIT_BY_SEX_EXTENDED_STR)   return FIT_BY_SEX_EXTENDED;
    if (s==FIT_BY_SEX_MAT_EXTENDED_STR) return FIT_BY_SEX_MAT_EXTENDED;
    std::cout<<"Unrecognized fit type '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return -1;
}

adstring tcsam::getFitType(int i){
    if (i==FIT_NONE)     return FIT_NONE_STR;
    if (i==FIT_BY_TOTAL) return FIT_BY_TOTAL_STR;
    if (i==FIT_BY_SEX)   return FIT_BY_SEX_STR;
    if (i==FIT_BY_SEX_EXTENDED)   return FIT_BY_SEX_EXTENDED_STR;
    if (i==FIT_BY_SEX_MAT_EXTENDED) return FIT_BY_SEX_MAT_EXTENDED_STR;
    std::cout<<"Unrecognized fit type '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
double tcsam::convertToStdDev(double sclVal, double mnVal, int sclFlg){
    double sdv = 0.0;
    if (sclFlg==SCLTYPE_CV) {
        sdv = sclVal*mnVal;
    } else if (sclFlg==SCLTYPE_STD) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
dvector tcsam::convertToStdDev(dvector sclVal, dvector mnVal, int sclFlg){
    dvector sdv = 0.0*sclVal;
    if (sclFlg==SCLTYPE_CV) {
        sdv = elem_prod(sclVal,mnVal);
    } else if (sclFlg==SCLTYPE_STD) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}
    
/**
 * Gets multiplicative conversion factor from "from" units to "to" units. 
 * @param from - adstring UNITS_ keyword
 * @param to   - adstring UNITS_ keyword
 * @return - multiplicative factor: to_units  = factor*from_units
 */
double tcsam::getConversionMultiplier(adstring from,adstring to){
    double cnv = 1.0;
    if (from==UNITS_ONES){
        if (to==UNITS_ONES)      return 1.0E0;
        if (to==UNITS_THOUSANDS) return 1.0E-3;
        if (to==UNITS_MILLIONS)  return 1.0E-6;
        if (to==UNITS_BILLIONS)  return 1.0E-9;
    }
    if (from==UNITS_THOUSANDS){
        if (to==UNITS_ONES)      return 1.0E+3;
        if (to==UNITS_THOUSANDS) return 1.0E+0;
        if (to==UNITS_MILLIONS)  return 1.0E-3;
        if (to==UNITS_BILLIONS)  return 1.0E-6;
    }
    if (from==UNITS_MILLIONS){
        if (to==UNITS_ONES)      return 1.0E6;
        if (to==UNITS_THOUSANDS) return 1.0E3;
        if (to==UNITS_MILLIONS)  return 1.0E0;
        if (to==UNITS_BILLIONS)  return 1.0E-3;
    }
    if (from==UNITS_BILLIONS){
        if (to==UNITS_ONES)      return 1.0E-9;
        if (to==UNITS_THOUSANDS) return 1.0E-6;
        if (to==UNITS_MILLIONS)  return 1.0E-3;
        if (to==UNITS_BILLIONS)  return 1.0E0;
    }
    if (from==UNITS_GM){
        if (to==UNITS_GM)   return 1.0E+0;
        if (to==UNITS_KG)   return 1.0E-3;
        if (to==UNITS_MT)   return 1.0E-6;
        if (to==UNITS_KMT)  return 1.0E-9;
        if (to==UNITS_LBS)  return 1.0E-3 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-9 * CONV_KGtoLBS;
    }
    if (from==UNITS_KG){
        if (to==UNITS_GM)   return 1.0E+3;
        if (to==UNITS_KG)   return 1.0E+0;
        if (to==UNITS_MT)   return 1.0E-3;
        if (to==UNITS_KMT)  return 1.0E-6;
        if (to==UNITS_LBS)  return 1.0E+0 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-6 * CONV_KGtoLBS;
    }
    if (from==UNITS_MT){
        if (to==UNITS_GM)   return 1.0E+6;
        if (to==UNITS_KG)   return 1.0E+3;
        if (to==UNITS_MT)   return 1.0E+0;
        if (to==UNITS_KMT)  return 1.0E-3;
        if (to==UNITS_LBS)  return 1.0E+3 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-3 * CONV_KGtoLBS;
    }
    if (from==UNITS_LBS){
        if (to==UNITS_GM)   return 1.0E+3/CONV_KGtoLBS;
        if (to==UNITS_KG)   return 1.0E+0/CONV_KGtoLBS;
        if (to==UNITS_MT)   return 1.0E-3/CONV_KGtoLBS;
        if (to==UNITS_KMT)  return 1.0E-6/CONV_KGtoLBS;
        if (to==UNITS_LBS)  return 1.0E+0;
        if (to==UNITS_MLBS) return 1.0E-6;
    }
    if (from==UNITS_MLBS){
        if (to==UNITS_GM)   return 1.0E+9/CONV_KGtoLBS;
        if (to==UNITS_KG)   return 1.0E+6/CONV_KGtoLBS;
        if (to==UNITS_MT)   return 1.0E+3/CONV_KGtoLBS;
        if (to==UNITS_KMT)  return 1.0E+0/CONV_KGtoLBS;
        if (to==UNITS_LBS)  return 1.0E+6;
        if (to==UNITS_MLBS) return 1.0E+0;
    }
    std::cout<<"Error converting units!"<<std::endl;
    std::cout<<"No conversion defined from '"<<from<<"' to '"<<to<<"'."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return cnv;
}

int tcsam::getLikelihoodType(adstring llType){
    int type = 0;
    if (llType==LL_LOGNORMAL_STR)   return LL_LOGNORMAL;
    if (llType==LL_MULTINOMIAL_STR) return LL_MULTINOMIAL;
    if (llType==LL_NONE_STR)        return LL_NONE;
    if (llType==LL_NORM2_STR)       return LL_NORM2;
    if (llType==LL_NORMAL_STR)      return LL_NORMAL;
    std::cout<<"Likelihood type keyword '"<<llType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return type;
}

adstring tcsam::getLikelihoodType(int llType){
    adstring type = LL_NONE_STR;
    if (llType==LL_LOGNORMAL)   return LL_LOGNORMAL_STR;
    if (llType==LL_MULTINOMIAL) return LL_MULTINOMIAL_STR;
    if (llType==LL_NONE)        return LL_NONE_STR;
    if (llType==LL_NORM2)       return LL_NORM2_STR;
    if (llType==LL_NORMAL)      return LL_NORMAL_STR;
    std::cout<<"Likelihood type integer '"<<llType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return type;
}
