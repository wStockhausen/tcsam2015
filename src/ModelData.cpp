#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelData.hpp"
#include "SummaryFunctions.hpp"

using namespace tcsam;

//**********************************************************************
//  Includes
//      SizeFrequencyData
//      BioData
//      ModelDatasets
//**********************************************************************
int AggregateCatchData::debug = 0;
int SizeFrequencyData::debug  = 0;
int BioData::debug            = 0;
int ModelDatasets::debug      = 0;
//----------------------------------------------------------------------
//          AggregateCatchData
//----------------------------------------------------------------------
const adstring AggregateCatchData::KW_ABUNDANCE_DATA = "AGGREGATE_ABUNDANCE";
const adstring AggregateCatchData::KW_BIOMASS_DATA   = "AGGREGATE_BIOMASS";

/*****************************************************************\n
 * Save the negative log-likelihoods from a model fit (values only).\n
 * \n
 * @param nlls - vector of negative log-likelihood components \n
 */
void AggregateCatchData::saveNLLs(dvar_vector& nlls){
    this->nlls = value(nlls);
}
/****************************************************************
 * Replace catch data C_xy with new data. Units are MILLIONS for 
 * abundance data and 1000's mt for biomass data.
 * Also modifies inpC_yc to reflect new data, but keeps original units.
 * Error-related quantities remain the same.
 * 
 * @param dmatrix newC_yx
 */
void AggregateCatchData::replaceCatchData(dmatrix& newC_yx){
    if (debug) std::cout<<"starting AggregateCatchData::replaceCatchData(dmatrix& newC_yx)"<<std::endl;
    int mnY = newC_yx.indexmin();
    int mxY = newC_yx.indexmax();
    C_xy.initialize();
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
 //        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
//        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
    for (int y=1;y<=ny;y++){
        int yr = yrs(y);
        if ((mnY<=yr)&&(yr<=mxY)) {
            for (int x=1;x<=nSXs;x++) C_xy(x,y) = newC_yx(yr,x);
            C_xy(ALL_SXs,y) = sum(newC_yx(yr));
            for (int x=1;x<=ALL_SXs;x++){
                inpC_yc(y,2*x)=C_xy(x,y)/convFac;            
            }
        }
    }
    if (debug) std::cout<<"finished AggregateCatchData::replaceCatchData(dmatrix& newC_yx)"<<std::endl;
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void AggregateCatchData::read(cifstream & is){
    if (debug) {
        std::cout<<"start AggregateCatchData::read(...) "<<this<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        std::cout<<"Apparent error reading AggregateCatchData."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<std::endl;
    if (str==KW_ABUNDANCE_DATA){type = KW_ABUNDANCE_DATA;} else
    if (str==KW_BIOMASS_DATA)  {type = KW_BIOMASS_DATA;}   else
    {   std::cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyword '"<<KW_ABUNDANCE_DATA<<"' or '"<<KW_BIOMASS_DATA<<"' but got '"<<str<<"'."<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<std::endl;
    is>>ny;//number of years of catch data
    rpt::echo<<ny<<tb<<"#number of years"<<std::endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units"<<std::endl;
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
    inpC_yc.allocate(1,ny,1,7);
    is>>inpC_yc;
    rpt::echo<<"#year males  cv_m  females  cv_f  total  cv_t"<<std::endl<<inpC_yc<<std::endl;
    
    yrs.allocate(1,ny);
    C_xy.allocate(1,ALL_SXs);
    cvs_xy.allocate(1,ALL_SXs);
    stdv_xy.allocate(1,ALL_SXs);
    yrs = (ivector) column(inpC_yc,1);
    for (int x=1;x<=ALL_SXs;x++) {
        C_xy(x)    = convFac*column(inpC_yc,2*x);
        cvs_xy(x)  = column(inpC_yc,2*x+1);
        stdv_xy(x) = sqrt(log(1.0+elem_prod(cvs_xy(x),cvs_xy(x))));
    }
    if (debug) std::cout<<"end AggregateCatchData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void AggregateCatchData::write(ostream & os){
    if (debug) std::cout<<"start AggregateCatchData::write(...) "<<this<<std::endl;
    os<<type<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<std::endl;
    os<<ny<<tb<<"#number of years of catch data"<<std::endl;
    os<<units<<tb<<"#units for catch data"<<std::endl;
    os<<"#year   males  cv_m  females  cv_f  total  cv_t"<<std::endl<<inpC_yc;
    if (debug) std::cout<<"end AggregateCatchData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void AggregateCatchData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) std::cout<<"AggregateCatchData::writing to R"<<std::endl;
    
    adstring x = qt+STR_MALE+qt+cc+qt+STR_FEMALE+qt+cc+qt+STR_ALL_SXs+qt;
    adstring y  = wts::to_qcsv(yrs);
    adstring str; adstring unitsp;
    if (type==KW_ABUNDANCE_DATA) {
        str="abundance"; 
        unitsp = tcsam::UNITS_MILLIONS;
    } else {
        str = "biomass";
        unitsp = tcsam::UNITS_KMT;
    }
    
    for (int n=0;n<indent;n++) os<<tb;
        os<<str<<"=list("<<std::endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc<<std::endl;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<std::endl;
            os<<"units="<<qt<<unitsp<<qt<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrs); os<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="; wts::writeToR(os,C_xy,x,y); os<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cvs="; wts::writeToR(os,cvs_xy,x,y); os<<std::endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")";
    if (debug) std::cout<<"AggregateCatchData::done writing to R"<<std::endl;
}
/////////////////////////////////end AggregateCatchData/////////////////////////
//----------------------------------------------------------------------
//          SizeFrequencyData
//----------------------------------------------------------------------
const adstring SizeFrequencyData::KW_SIZEFREQUENCY_DATA = "SIZE_FREQUENCY_DATA";
/**
 * Save the negative log-likelihoods from a model fit (values only).
 * 
 * @param nlls - matrix of negativie log-likelihood components
 */
void SizeFrequencyData::saveNLLs(dvar_matrix& nlls){
    this->nlls = value(nlls);
}
/**
 * Normalize the size frequency data to sum to 1 over x,m,s,z.
 */
void SizeFrequencyData::normalize(void){
    dvector nT(1,ny); nT.initialize();
    for (int y=1;y<=ny;y++){
        //calculate total numbers
        for (int x=1;x<=ALL_SXs;x++){
            for (int m=1;m<=ALL_MSs;m++){
                for (int s=1;s<=ALL_SCs;s++) nT(y) += sum(NatZ_xmsyz(x,m,s,y));//sum over size bins
            }
        }
        //normalize numbers-at-size by total numbers
        for (int x=1;x<=ALL_SXs;x++){
            for (int m=1;m<=ALL_MSs;m++){
                for (int s=1;s<=ALL_SCs;s++) PatZ_xmsyz(x,m,s,y) = NatZ_xmsyz(x,m,s,y)/nT(y);//norm over size bins
            }
        }
    }
}

/*******************************************************\n
 * Replace catch-at-size data NatZ_xmsyz with new data. 
 * Also modifies inpNatZ_xmsyc to reflect new data.
 * Error-related quantities remain the same.
 * 
 * @param d5_array newNatZ_yxmsz
 *****************************************************/
void SizeFrequencyData::replaceSizeFrequencyData(d5_array& newNatZ_yxmsz){
    if (debug) std::cout<<"starting SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<std::endl;
    int mnY = newNatZ_yxmsz.indexmin();
    int mxY = newNatZ_yxmsz.indexmax();
    //zero out old data
    for (int y=1;y<=ny;y++){
        int yr = yrs(y);
        if ((mnY<=yr)&&(yr<=mxY)){
            for (int x=1;x<=ALL_SXs;x++){
                for (int m=1;m<=ALL_MSs;m++){
                    for (int s=1;s<=ALL_SCs;s++) {
                        NatZ_xmsyz(x,m,s,y) = 0.0;
                        inpNatZ_xmsyc(x,m,s,y)(3,2+nZCs-1) = 0.0;
                    }
                }
            }
        }
    }
    //put in new data only where it existed before
    int nc = factors.indexmax();
    for (int y=1;y<=ny;y++){
        int yr = yrs(y);
        if ((mnY<=yr)&&(yr<=mxY)){
            for (int c=1;c<=nc;c++){
                int x = tcsam::getSexType(factors(c,1));
                int m = tcsam::getMaturityType(factors(c,2));
                int s = tcsam::getShellType(factors(c,3));
                dvector n_z = extractFromYXMSZ(yr,x,m,s,newNatZ_yxmsz);
                NatZ_xmsyz(x,m,s,y) = n_z;
                inpNatZ_xmsyc(x,m,s,y)(3,2+nZCs-1).shift(1) = n_z;
//                            std::cout<<"Bounds NatZ_xmsyz(x,m,s,y)   : "<<wts::getBounds(NatZ_xmsyz(x,m,s,y))<<std::endl;
//                            std::cout<<"Bounds inpNatZ_xmsyc(x,m,s,y): "<<wts::getBounds(inpNatZ_xmsyc(x,m,s,y))<<std::endl;
            }
        }
    }
    normalize();
    if (debug) std::cout<<"end SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<std::endl;
}

/*******************************************************\n
*   read from input stream.\n
******************************************************/
void SizeFrequencyData::read(cifstream & is){
    if (debug) {
        std::cout<<"start SizeFrequencyData::read(...) "<<this<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        std::cout<<"Apparent error reading SizeFrequencyData."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<std::endl;
    if (!(str==KW_SIZEFREQUENCY_DATA)){
        std::cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyowrd '"<<KW_SIZEFREQUENCY_DATA<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    
    //NUMBERS-AT-SIZE 
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>ny;//number of years of numbers-at-size data
    rpt::echo<<ny<<tb<<"#number of years for size frequency data"<<std::endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units for numbers at size data"<<std::endl;
    is>>nZCs;
    rpt::echo<<nZCs<<tb<<"#number of size bin cut points (nZCs)"<<std::endl;
    zCs.allocate(1,nZCs);
    is>>zCs;
    rpt::echo<<zCs<<tb<<"#zCutPts (mm CW)"<<std::endl;
    zBs.allocate(1,nZCs-1);
    zBs = 0.5*(zCs(1,nZCs-1)+(--zCs(2,nZCs)));
    if (debug) std::cout<<zBs<<tb<<"#zBins"<<std::endl;
    
    yrs.allocate(1,ny);
    ss_xmsy.allocate(1,ALL_SXs,1,ALL_MSs,1,ALL_SCs,1,ny);
    NatZ_xmsyz.allocate(1,ALL_SXs,1,ALL_MSs,1,ALL_SCs,1,ny,1,nZCs-1);
    PatZ_xmsyz.allocate(1,ALL_SXs,1,ALL_MSs,1,ALL_SCs,1,ny,1,nZCs-1);
    inpNatZ_xmsyc.allocate(1,ALL_SXs,1,ALL_MSs,1,ALL_SCs,1,ny,1,2+(nZCs-1));
    
    yrs.initialize();
    ss_xmsy.initialize();
    NatZ_xmsyz.initialize();
    PatZ_xmsyz.initialize();
    inpNatZ_xmsyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    rpt::echo<<nc<<tb<<"#number of factor combinations to read"<<std::endl;
    factors.allocate(1,nc,1,3);
    for (int i=0;i<nc;i++){
        is>>factors(i+1);
        int x = tcsam::getSexType(factors(i+1,1));
        int m = tcsam::getMaturityType(factors(i+1,2));
        int s = tcsam::getShellType(factors(i+1,3));
        rpt::echo<<factors(i+1,1)<<tb<<factors(i+1,2)<<tb<<factors(i+1,3)<<tb<<"#factors"<<std::endl;
        if (x&&m&&s){
            is>>inpNatZ_xmsyc(x,m,s);
            rpt::echo<<"#year    ss     "<<zBs<<std::endl;
            rpt::echo<<inpNatZ_xmsyc(x,m,s)<<std::endl;
            yrs = (ivector) column(inpNatZ_xmsyc(x,m,s),1);
            ss_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
            for (int y=1;y<=ny;y++){
                NatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(3,2+(nZCs-1))).shift(1);
//                std::cout<<"Bounds NatZ_xmsyz(x,m,s,y)   : "<<wts::getBounds(NatZ_xmsyz(x,m,s,y))<<std::endl;
//                std::cout<<"Bounds inpNatZ_xmsyc(x,m,s,y): "<<wts::getBounds(inpNatZ_xmsyc(x,m,s,y))<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"At least one factor for NatZ not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(i+1,1)<<tb<<factors(i+1,2)<<tb<<factors(i+1,3)<<std::endl;
            std::cout<<"Aborting..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }
    normalize();
    if (debug) std::cout<<"end SizeFrequencyData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void SizeFrequencyData::write(ostream & os){
    if (debug) std::cout<<"start SizeFrequencyData::write(...) "<<this<<std::endl;
    os<<KW_SIZEFREQUENCY_DATA<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<ny<<tb<<"#number of years of size data"<<std::endl;
    os<<units<<tb<<"#units for numbers at size"<<std::endl;
    os<<nZCs<<tb<<"#number of size bin cutpoints"<<std::endl;
    os<<"#size bin cutpoints (mm CW)"<<std::endl<<zCs<<std::endl;
    
//    int nCs = 0;
//    for (int x=1;x<=ALL_SXs;x++){
//        for (int m=1;m<=ALL_MSs;m++){
//            for (int s=1;s<=ALL_SCs;s++){
//                if (sum(NatZ_xmsyz(x,m,s))>0) nCs++;
//            }
//        }
//    }
    
    int nCs = factors.indexmax();
    os<<nCs<<tb<<"#number of sex x shell x maturity factor combinations to read in"<<std::endl;
//    adstring_array factors(1,3);
//    for (int x=1;x<=ALL_SXs;x++){
//        factors(1) = tcsam::getSexType(x);
//        for (int m=1;m<=ALL_MSs;m++){
//            factors(2) = tcsam::getMaturityType(m);
//            for (int s=1;s<=ALL_SCs;s++){
//                factors(3) = tcsam::getShellType(s);
//                if (sum(NatZ_xmsyz(x,m,s))>0){ //only print out non-zero matrices
//                    os<<"#-------"<<factors(1)<<cc<<factors(2)<<cc<<factors(3)<<std::endl;
//                    os<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<std::endl;
//                    os<<"#year    ss     "<<zBs<<std::endl;
//                    os<<inpNatZ_xmsyc(x,m,s)<<std::endl;
//                }
//            }
//        }
    for (int c=1;c<=nCs;c++){
        int x = tcsam::getSexType(factors(c,1));
        int m = tcsam::getMaturityType(factors(c,2));
        int s = tcsam::getShellType(factors(c,3));
        os<<"#-------"<<factors(c,1)<<cc<<factors(c,2)<<cc<<factors(c,3)<<std::endl;
        os<<factors(c,1)<<tb<<factors(c,2)<<tb<<factors(c,3)<<std::endl;
        os<<"#year    ss     "<<zBs<<std::endl;
        os<<inpNatZ_xmsyc(x,m,s)<<std::endl;
    }
    if (debug) std::cout<<"end SizeFrequencyData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void SizeFrequencyData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) std::cout<<"SizeFrequencyData::writing to R"<<std::endl;
    adstring x = qt+tcsam::STR_MALE   +qt+cc+qt+tcsam::STR_FEMALE     +qt+cc+qt+tcsam::STR_ALL_SXs+qt;
    adstring m = qt+tcsam::STR_IMMATURE +qt+cc+qt+tcsam::STR_MATURE   +qt+cc+qt+tcsam::STR_ALL_MSs+qt;
    adstring s = qt+tcsam::STR_NEW_SHELL+qt+cc+qt+tcsam::STR_OLD_SHELL+qt+cc+qt+tcsam::STR_ALL_SCs+qt;
    adstring y = wts::to_qcsv(yrs);
    adstring z = wts::to_qcsv(zBs);        
    for (int n=0;n<indent;n++) os<<tb;
        os<<"nAtZ=list(units="<<qt<<units<<qt<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"years="; wts::writeToR(os,yrs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"cutpts="; wts::writeToR(os,zCs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"sample.sizes="; wts::writeToR(os,ss_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<std::endl;
        wts::writeToR(os,NatZ_xmsyz,x,m,s,y,z); os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) std::cout<<"SizeFrequencyData::done writing to R"<<std::endl;
}
/////////////////////////////////end SizeFrequencyData/////////////////////////
//----------------------------------------------------------------------
//          BioData
//----------------------------------------------------------------------
const adstring BioData::KW_BIO_DATA = "BIO_DATA";
/***************************************************************
*   read.                                                      *
***************************************************************/
void BioData::read(cifstream & is){
    if (debug) std::cout<<"start BioData::read(...) "<<this<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    rpt::echo<<"#BioData"<<std::endl;
    rpt::echo<<"#file name is "<<is.get_file_name()<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    if (!is) {
        std::cout<<"Apparent error reading Bio Data."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_BIO_DATA)){
        std::cout<<"#Error reading bio data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyword '"<<KW_BIO_DATA<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    rpt::echo<<str<<tb<<"#required keyword"<<std::endl;
    
    //RECRUITMENT LAG
    is>>recLag;
    rpt::echo<<recLag<<tb<<"#recLag"<<std::endl;
    
    //SIZE BINS
    is>>nZBins;
    rpt::echo<<nZBins<<tb<<"#nZBins"<<std::endl;
    zBins.allocate(1,nZBins);
    is>>zBins;
    rpt::echo<<zBins<<tb<<"#zBins (mm CW)"<<std::endl;
    
    //WEIGHT-AT-SIZE
    {is>>unitsWatZ;
    rpt::echo<<unitsWatZ<<tb<<"#unitsWatZ"<<std::endl;
    double convToKG = tcsam::getConversionMultiplier(unitsWatZ,tcsam::UNITS_KG);
    rpt::echo<<"#using conversion factor from "<<unitsWatZ<<" to kg: "<<convToKG<<std::endl;
    wAtZ_xmz.allocate(1,nSXs,1,nMSs,1,nZBins);
    wAtZ_xmz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<std::endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = tcsam::getSexType(factors(1));
        int mat   = tcsam::getMaturityType(factors(2));
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<std::endl;
        if (sex||mat){
            is>>wAtZ_xmz(sex,mat);
            rpt::echo<<wAtZ_xmz(sex,mat)<<std::endl;
            wAtZ_xmz(sex,mat) *= convToKG;//convert to kg
            if (debug) {
                std::cout<<"#wAtZ_xmz("<<factors(1)<<cc<<factors(2)<<") [kg] ="<<std::endl;
                std::cout<<"#"<<zBins<<std::endl;
                std::cout<<wAtZ_xmz(sex,mat)<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"Factors for wAtZ not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<std::endl;
            std::cout<<"Terminating..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }}
    
    //PROBABILITY OF MATURING AT SIZE
    {prMature_xz.allocate(1,nSXs,1,nZBins);
    prMature_xz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<std::endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<tb<<"#factors"<<std::endl;
        int sex   = tcsam::getSexType(factors(1));
        if (sex){
            is>>prMature_xz(sex);
            rpt::echo<<prMature_xz(sex)<<std::endl;
            if (debug) {
                std::cout<<"#prMature_xz("<<factors(1)<<")="<<std::endl;
                std::cout<<"#"<<zBins<<std::endl;
                std::cout<<prMature_xz(sex)<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"Factors for prMature_xz not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(1)<<std::endl;
            std::cout<<"Terminating..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }}
    
    //FRACTION MATURE BY SEX, SHELL CONDITION
    {frMature_xsz.allocate(1,nSXs,1,nMSs,1,nZBins);
    frMature_xsz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<std::endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<std::endl;
        int sex   = tcsam::getSexType(factors(1));
        int shl   = tcsam::getShellType(factors(2));
        if (sex||shl){
            is>>frMature_xsz(sex,shl);
            rpt::echo<<frMature_xsz(sex,shl)<<std::endl;
            if (debug) {
                std::cout<<"#frMature_xsz("<<factors(1)<<cc<<factors(2)<<")="<<std::endl;
                std::cout<<"#"<<zBins<<std::endl;
                std::cout<<frMature_xsz(sex,shl)<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"Factors for frMature_xsz not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<std::endl;
            std::cout<<"Terminating..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }}
    
    //CV for mean size at min, max size bins
    {rpt::echo<<"#CV for min, max sizes"<<std::endl;
    cvMnMxZ_xc.allocate(1,nSXs,1,2);
    cvMnMxZ_xc.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<std::endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<"#factors"<<std::endl;
        int sex   = tcsam::getSexType(factors(1));
        if (sex){
            is>>cvMnMxZ_xc(sex);
            rpt::echo<<cvMnMxZ_xc(sex)<<std::endl;
            if (debug) {
                std::cout<<"#cvMnMxZ_xc("<<factors(1)<<")="<<std::endl;
                std::cout<<cvMnMxZ_xc(sex)<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"Factors for cvMnMxZ_xc not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(1)<<std::endl;
            std::cout<<"Terminating..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }}
    
    //MIDPOINTS OF FISHERY SEASONS BY YEAR
    rpt::echo<<"#Timing of fisheries and mating"<<std::endl;
    is>>nyTiming;
    rpt::echo<<nyTiming<<tb<<"#number of years"<<std::endl;
    timing_yc.allocate(1,nyTiming,1,3);
    is>>timing_yc;
    rpt::echo<<"#timing"<<std::endl<<timing_yc<<std::endl;
    fshTiming_y.allocate(timing_yc(1)(1),timing_yc(nyTiming)(1));
    matTiming_y.allocate(timing_yc(1)(1),timing_yc(nyTiming)(1));
    for (int i=1;i<=nyTiming;i++){
        fshTiming_y(timing_yc(i,1)) = timing_yc(i,2);
        matTiming_y(timing_yc(i,1)) = timing_yc(i,3);
    }
    
    if (debug) std::cout<<"end BioData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void BioData::write(ostream & os){
    if (debug) std::cout<<"start BioData::write(...) "<<this<<std::endl;
    
    os<<KW_BIO_DATA<<std::endl;
    
    os<<"#-----------RECRUITMENT LAG---------------------------#"<<std::endl;
    os<<nZBins<<tb<<"#recLag (recruitment lag in years)"<<std::endl;
    
    os<<"#-----------SIZE BINS---------------------------------#"<<std::endl;
    os<<nZBins<<tb<<"#number of size bins"<<std::endl;
    os<<"#size bins (mm CW)"<<std::endl<<zBins<<std::endl;
    
    {os<<"#-----------WEIGHT-AT-SIZE----------------------------#"<<std::endl;
    os<<unitsWatZ<<tb<<"#units for weight-at-size"<<std::endl;
    os<<nSXs*nMSs<<tb<<"#number of factor combinations (sex x maturity state)"<<std::endl;
    adstring_array factors(1,2);
    for (int sex=1;sex<=nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int mat=1;mat<=nMSs;mat++){
            factors(2) = tcsam::getMaturityType(mat);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<std::endl;
            os<<factors(1)<<tb<<factors(2)<<std::endl;
            os<<wAtZ_xmz(sex,mat)<<std::endl;
        }
    }}
    
    {os<<"#-----------PROBABILITY OF MATURING-------------------#"<<std::endl;
    os<<nSXs<<tb<<"#number of factor combinations (sex)"<<std::endl;
    adstring_array factors(1,1);
    for (int sex=1;sex<=nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        os<<"#-------"<<factors(1)<<std::endl;
        os<<factors(1)<<std::endl;
        os<<prMature_xz(sex)<<std::endl;
    }}
    
    {os<<"#-----------FRACTION MATURE BY SEX, SHELL CONDITION---#"<<std::endl;
    os<<nSXs*nSCs<<tb<<"#number of factor combinations (sex x shell condition)"<<std::endl;
    adstring_array factors(1,2);
    for (int sex=1;sex<=nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int shl=1;shl<=nSCs;shl++){
            factors(2) = tcsam::getShellType(shl);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<std::endl;
            os<<factors(1)<<tb<<factors(2)<<std::endl;
            os<<frMature_xsz(sex,shl)<<std::endl;
        }
    }}
    
    {os<<"#-----------CV FOR MIN, MAX SIZES---------------------#"<<std::endl;
    os<<nSXs<<tb<<"#number of factor combinations (sex)"<<std::endl;
    adstring_array factors(1,1);
    for (int sex=1;sex<=nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        os<<"#-------"<<factors(1)<<std::endl;
        os<<factors(1)<<std::endl;
        os<<cvMnMxZ_xc(sex)<<std::endl;
    }}
    
    {os<<"#-----------TIMING------------------------------------#"<<std::endl;
    os<<nyTiming<<tb<<"#number of years"<<std::endl;
    os<<"#year   midptFisheries  matingTime"<<std::endl;
    os<<timing_yc<<std::endl;}
    
    if (debug) std::cout<<"end BioData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void BioData::writeToR(ostream& os, string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<std::endl;
    indent++;
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"recLag="<<recLag<<","<<std::endl;
        
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"zBins=";wts::writeToR(os,zBins); os<<","<<std::endl;
        
        //weight-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"wAtZ=list(units="<<qt<<unitsWatZ<<qt<<cc<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<std::endl;
            adstring x = qt+tcsam::STR_MALE+qt    +cc+ qt+tcsam::STR_FEMALE+qt;
            adstring s = qt+tcsam::STR_IMMATURE+qt+cc+ qt+tcsam::STR_MATURE+qt;
            adstring c = wts::to_qcsv(zBins);
            wts::writeToR(os,wAtZ_xmz,x,s,c); os<<std::endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<std::endl;
        
        //probability of maturing-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"prMat="<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_MALE+qt  +cc+ qt+tcsam::STR_FEMALE+qt;
            adstring z = wts::to_qcsv(zBins);
            wts::writeToR(os,prMature_xz,x,z); os<<","<<std::endl;
        indent--;}
        
        //fraction mature-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"frMat="<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_MALE+qt     +cc+ qt+tcsam::STR_FEMALE+qt;
            adstring s = qt+tcsam::STR_NEW_SHELL+qt+cc+ qt+tcsam::STR_OLD_SHELL+qt;
            adstring z = wts::to_qcsv(zBins);
            wts::writeToR(os,frMature_xsz,x,s,z); os<<","<<std::endl;
        indent--;}
        
        //cv for min, max sizes
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"cvZs="<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_MALE+qt  +cc+ qt+tcsam::STR_FEMALE+qt;
            adstring c = "'minZ','maxZ'";
            wts::writeToR(os,cvMnMxZ_xc,x,c); os<<","<<std::endl;
        indent--;}
        
        //timing of fishery season midpoints and mating
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"timing="<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring cols = "'midptFisheries','matingTime'";
            ivector yrs = wts::to_ivector(column(timing_yc,1));
            dmatrix tmp(1,2,1,nyTiming);
            tmp(1) = column(timing_yc,2);
            tmp(2) = column(timing_yc,3);
            wts::writeToR(os,trans(tmp),yrs,cols); os<<std::endl;
        indent--;}
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end BioData/////////////////////////
//----------------------------------------------------------------------
//          ModelDatasets
//----------------------------------------------------------------------
/***************************************************************
*   Instantiation                                              *
***************************************************************/
ModelDatasets::ModelDatasets(ModelConfiguration* ptrMC){
    pMC=ptrMC;
    ptrBio=0;
//    ppFsh=0;
    ppSrv=0;
}
/***************************************************************
*   Destruction                                                *
***************************************************************/
ModelDatasets::~ModelDatasets(){
    pMC=0;
    delete ptrBio;  ptrBio=0;
//    if (ppFsh) {                                                        TODO: uncomment
//        for (int f=0;f<nFsh;f++) delete ppFsh[f];
//        delete ppFsh; ppFsh = 0;
//    } 
    if (ppSrv) {
        for (int s=0;s<nSrv;s++) delete ppSrv[s];
        delete ppSrv; ppSrv = 0;
    } 
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::read(cifstream & is){
    rpt::echo<<"#------------------------------------------"<<std::endl;
    rpt::echo<<"reading Model Datasets file."<<std::endl;
    rpt::echo<<"#file name is '"<<is.get_file_name()<<"'"<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    is>>fnBioData;
    rpt::echo<<fnBioData<<tb<<"#bio data filename"<<std::endl;
    is>>nFsh;
    rpt::echo<<nFsh<<tb<<"#number of fishery datasets to read in"<<std::endl;
    if (nFsh){
        fnsFisheryData.allocate(1,nFsh);
        is>>fnsFisheryData; 
        for (int i=1;i<=nFsh;i++) rpt::echo<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<std::endl;
    }
    is>>nSrv;
    rpt::echo<<nSrv<<tb<<"#number of survey datasets to read in"<<std::endl;
    if (nSrv){
        fnsSurveyData.allocate(1,nSrv);
        is>>fnsSurveyData; 
        for (int i=1;i<=nSrv;i++) rpt::echo<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<std::endl;
    }
    
    //          Bio data
    {ptrBio = new BioData();
    cifstream strm(fnBioData,ios::in);
    rpt::echo<<"#------Bio data ------------"<<std::endl;
    strm>>(*ptrBio);
    }  
    
    //          Fishery data 
    if (nFsh) {
        rpt::echo<<"#-------Fishery Datasets---------"<<std::endl;
        ppFsh = new FisheryData*[nFsh];
        for (int i=0;i<nFsh;i++) {
            ppFsh[i] = new FisheryData();
            cifstream strm(fnsFisheryData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Fishery Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppFsh[i]);
        }
    }
    
    //          Survey data
    if (nSrv) {
        rpt::echo<<"#-------Survey Datasets---------"<<std::endl;
        ppSrv = new SurveyData*[nSrv];
        for (int i=0;i<nSrv;i++) {
            ppSrv[i] = new SurveyData();
            cifstream strm(fnsSurveyData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Survey Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppSrv[i]);
        }
    }
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::write(ostream & os){
    if (debug) std::cout<<"start ModelDatasets::write(...) "<<this<<std::endl;
    os<<fnBioData<<tb<<"#tanner crab biological data file"<<std::endl;
    os<<"#-------fishery data files---------"<<std::endl;
    os<<nFsh<<tb<<"#number of fishery data files"<<std::endl;
    for (int i=1;i<+nFsh;i++) os<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<std::endl;
    os<<"#-------survey data files---------"<<std::endl;
    os<<nSrv<<tb<<"#number of survey data files"<<std::endl;
    for (int i=1;i<+nSrv;i++) os<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<std::endl;
    os<<"#-----biological data---------"<<std::endl;
    os<<(*ptrBio)<<std::endl;
    os<<"#-------fishery data ---------"<<std::endl;
    if (nFsh){
        for (int i=1;i<nFsh;i++) os<<"#----fishery dataset "<<i<<std::endl<<(*ppFsh[i-1])<<std::endl;
        os<<"#----fishery dataset "<<nFsh<<std::endl<<(*ppFsh[nFsh-1])<<std::endl;
    }
    os<<"#-------survey data ---------"<<std::endl;
    if (nSrv){
        for (int i=1;i<nSrv;i++) os<<"#----survey dataset "<<i<<std::endl<<(*ppSrv[i-1])<<std::endl;
        os<<"#----survey dataset "<<nSrv<<std::endl<<(*ppSrv[nSrv-1]);
    }
   if (debug) std::cout<<"end ModelDatasets::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelDatasets::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list("<<std::endl;
    indent++;
        //bio data as list
            ptrBio->writeToR(os,"bio",indent); os<<cc<<std::endl;
        
        //survey data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"surveys=list("<<std::endl;
        indent++;
            if (ppSrv) {
                for (int i=0;i<(nSrv-1);i++) {
                    ppSrv[i]->writeToR(os,(char*)ppSrv[i]->name,indent); os<<cc<<std::endl;
                }
                ppSrv[nSrv-1]->writeToR(os,(char*)ppSrv[nSrv-1]->name,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<std::endl;            
        
        //fishery data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"fisheries=list("<<std::endl;
        indent++;
            if (ppFsh) {
                for (int i=0;i<(nFsh-1);i++) {
                    ppFsh[i]->writeToR(os,(char*)ppFsh[i]->name,indent); os<<cc<<std::endl;
                }
                ppFsh[nFsh-1]->writeToR(os,(char*)ppFsh[nFsh-1]->name,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<")"<<std::endl;            
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end ModelDatasets/////////////////////////
