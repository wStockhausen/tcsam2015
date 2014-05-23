#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelData.hpp"

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
    if (debug) cout<<"starting AggregateCatchData::replaceCatchData(dmatrix& newC_yx)"<<endl;
    int mnY = newC_yx.indexmin();
    int mxY = newC_yx.indexmax();
    C_xy.initialize();
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
 //        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
//        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<endl;
    }
    for (int y=1;y<=ny;y++){
        int yr = yrs(y);
        if ((mnY<=yr)&&(yr<=mxY)) {
            for (int x=1;x<=nSXs;x++) C_xy(x,y) = newC_yx(yr,x);
            C_xy(ANY_SEX,y) = sum(newC_yx(yr));
            for (int x=1;x<=(nSXs+1);x++){
                inpC_yc(y,2*x)=C_xy(x,y)/convFac;            
            }
        }
    }
    if (debug) cout<<"finished AggregateCatchData::replaceCatchData(dmatrix& newC_yx)"<<endl;
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void AggregateCatchData::read(cifstream & is){
    if (debug) {
        cout<<"start AggregateCatchData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading AggregateCatchData."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<endl;
    if (str==KW_ABUNDANCE_DATA){type = KW_ABUNDANCE_DATA;} else
    if (str==KW_BIOMASS_DATA)  {type = KW_BIOMASS_DATA;}   else
    {   cout<<"#Error reading effort data from "<<is.get_file_name()<<endl;
        cout<<"Expected keyowrd '"<<KW_ABUNDANCE_DATA<<"' or '"<<KW_BIOMASS_DATA<<"' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<endl;
    is>>ny;//number of years of catch data
    rpt::echo<<ny<<tb<<"#number of years"<<endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units"<<endl;
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<endl;
    }
    inpC_yc.allocate(1,ny,1,7);
    is>>inpC_yc;
    rpt::echo<<"#year females  cv_f  males  cv_m  total  cv_t"<<endl<<inpC_yc<<endl;
    
    yrs.allocate(1,ny);
    C_xy.allocate(1,nSXs+1);
    cvs_xy.allocate(1,nSXs+1);
    stdv_xy.allocate(1,nSXs+1);
    yrs = (ivector) column(inpC_yc,1);
    for (int x=1;x<=(nSXs+1);x++) {
        C_xy(x)    = convFac*column(inpC_yc,2*x);
        cvs_xy(x)  = column(inpC_yc,2*x+1);
        stdv_xy(x) = sqrt(log(1.0+elem_prod(cvs_xy(x),cvs_xy(x))));
    }
    if (debug) cout<<"end AggregateCatchData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void AggregateCatchData::write(ostream & os){
    if (debug) cout<<"start AggregateCatchData::write(...) "<<this<<endl;
    os<<type<<tb<<"#required keyword"<<endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<endl;
    os<<tcsam::getLikelihoodType(optFit)<<tb<<"#likelihood type"<<endl;
    os<<ny<<tb<<"#number of years of catch data"<<endl;
    os<<units<<tb<<"#units for catch data"<<endl;
    os<<"#year   females  cv_f  males  cv_m  total  cv_t"<<endl<<inpC_yc;
    if (debug) cout<<"end AggregateCatchData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void AggregateCatchData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"AggregateCatchData::writing to R"<<endl;
    
    adstring x = qt+STR_FEMALE+qt+cc+qt+STR_MALE+qt+cc+qt+STR_ANY_SEX+qt;
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
        os<<str<<"=list("<<endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc<<endl;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<endl;
            os<<"units="<<qt<<unitsp<<qt<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrs); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="; wts::writeToR(os,C_xy,x,y); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cvs="; wts::writeToR(os,cvs_xy,x,y); os<<endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")";
    if (debug) cout<<"AggregateCatchData::done writing to R"<<endl;
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
        for (int x=1;x<=(nSXs+1);x++){
            for (int m=1;m<=(nMSs+1);m++){
                for (int s=1;s<=(nSCs+1);s++) nT(y) += sum(NatZ_xmsyz(x,m,s,y));//sum over size bins
            }
        }
        //normalize numbers-at-size by total numbers
        for (int x=1;x<=(nSXs+1);x++){
            for (int m=1;m<=(nMSs+1);m++){
                for (int s=1;s<=(nSCs+1);s++) PatZ_xmsyz(x,m,s,y) = NatZ_xmsyz(x,m,s,y)/nT(y);//norm over size bins
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
    if (debug) cout<<"starting SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<endl;
    int mnY = newNatZ_yxmsz.indexmin();
    int mxY = newNatZ_yxmsz.indexmax();
    for (int y=1;y<=ny;y++){
        int yr = yrs(y);
        if ((mnY<=yr)&&(yr<=mxY)){
            for (int x=1;x<=(nSXs+1);x++){
                for (int m=1;m<=(nMSs+1);m++){
                    for (int s=1;s<=(nSCs+1);s++) {
                        if ((x<=nSXs)&&(m<=nMSs)&&(s<=nSCs)){
                            NatZ_xmsyz(x,m,s,y) = newNatZ_yxmsz(yr,x,m,s);
                            inpNatZ_xmsyc(x,m,s,y)(3,2+nZCs-1).shift(1) = newNatZ_yxmsz(yr,x,m,s);
//                            cout<<"Bounds NatZ_xmsyz(x,m,s,y)   : "<<wts::getBounds(NatZ_xmsyz(x,m,s,y))<<endl;
//                            cout<<"Bounds inpNatZ_xmsyc(x,m,s,y): "<<wts::getBounds(inpNatZ_xmsyc(x,m,s,y))<<endl;
                        } else{
                            NatZ_xmsyz(x,m,s,y) = 0.0;//summary of some sort
                            inpNatZ_xmsyc(x,m,s,y)(3,2+nZCs-1) = 0.0;
                        }
                    }
                }
            }
        } else {
            for (int x=1;x<=(nSXs+1);x++){
                for (int m=1;m<=(nMSs+1);m++){
                    for (int s=1;s<=(nSCs+1);s++) {
                        NatZ_xmsyz(x,m,s,y) = 0.0;
                        inpNatZ_xmsyc(x,m,s,y)(3,2+nZCs-1) = 0.0;
                    }
                }
            }
        }
    }
    normalize();
    if (debug) cout<<"end SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<endl;
}

/*******************************************************\n
*   read from input stream.\n
******************************************************/
void SizeFrequencyData::read(cifstream & is){
    if (debug) {
        cout<<"start SizeFrequencyData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading SizeFrequencyData."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<endl;
    if (!(str==KW_SIZEFREQUENCY_DATA)){
        cout<<"#Error reading effort data from "<<is.get_file_name()<<endl;
        cout<<"Expected keyowrd '"<<KW_SIZEFREQUENCY_DATA<<"' but got '"<<str<<"'"<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    //NUMBERS-AT-SIZE 
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<endl;
    is>>ny;//number of years of numbers-at-size data
    rpt::echo<<ny<<tb<<"#number of years for size frequency data"<<endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units for numbers at size data"<<endl;
    is>>nZCs;
    rpt::echo<<nZCs<<tb<<"#number of size bin cut points (nZCs)"<<endl;
    zCs.allocate(1,nZCs);
    is>>zCs;
    rpt::echo<<zCs<<tb<<"#zCutPts (mm CW)"<<endl;
    zBs.allocate(1,nZCs-1);
    zBs = 0.5*(zCs(1,nZCs-1)+(--zCs(2,nZCs)));
    if (debug) cout<<zBs<<tb<<"#zBins"<<endl;
    
    yrs.allocate(1,ny);
    ss_xmsy.allocate(1,nSXs+1,1,nMSs+1,1,nSCs+1,1,ny);
    NatZ_xmsyz.allocate(1,nSXs+1,1,nMSs+1,1,nSCs+1,1,ny,1,nZCs-1);
    PatZ_xmsyz.allocate(1,nSXs+1,1,nMSs+1,1,nSCs+1,1,ny,1,nZCs-1);
    inpNatZ_xmsyc.allocate(1,nSXs+1,1,nMSs+1,1,nSCs+1,1,ny,1,2+(nZCs-1));
    
    yrs.initialize();
    ss_xmsy.initialize();
    NatZ_xmsyz.initialize();
    PatZ_xmsyz.initialize();
    inpNatZ_xmsyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    rpt::echo<<nc<<tb<<"#number of factor combinations to read"<<endl;
    adstring_array factors(1,3);
    for (int i=0;i<nc;i++){
        is>>factors;
        int x  = tcsam::getSexType(factors(1));
        int m  = tcsam::getMaturityType(factors(2));
        int s = tcsam::getShellType(factors(3));
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<tb<<"#factors"<<endl;
        if (x&&m&&s){
            is>>inpNatZ_xmsyc(x,m,s);
            rpt::echo<<"#year    ss     "<<zBs<<endl;
            rpt::echo<<inpNatZ_xmsyc(x,m,s)<<endl;
            yrs = (ivector) column(inpNatZ_xmsyc(x,m,s),1);
            ss_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
            for (int y=1;y<=ny;y++){
                NatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(3,2+(nZCs-1))).shift(1);
//                cout<<"Bounds NatZ_xmsyz(x,m,s,y)   : "<<wts::getBounds(NatZ_xmsyz(x,m,s,y))<<endl;
//                cout<<"Bounds inpNatZ_xmsyc(x,m,s,y): "<<wts::getBounds(inpNatZ_xmsyc(x,m,s,y))<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"At least one factor for NatZ not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<endl;
            cout<<"Aborting..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }
    normalize();
    if (debug) cout<<"end SizeFrequencyData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void SizeFrequencyData::write(ostream & os){
    if (debug) cout<<"start SizeFrequencyData::write(...) "<<this<<endl;
    os<<KW_SIZEFREQUENCY_DATA<<tb<<"#required keyword"<<endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<endl;
    os<<ny<<tb<<"#number of years of size data"<<endl;
    os<<units<<tb<<"#units for numbers at size"<<endl;
    os<<nZCs<<tb<<"#number of size bin cutpoints"<<endl;
    os<<"#size bin cutpoints (mm CW)"<<endl<<zCs<<endl;
    
    int nCs = 0;
    for (int x=1;x<=ANY_SEX;x++){
        for (int m=1;m<=ANY_MATURITY;m++){
            for (int s=1;s<=ANY_SHELL;s++){
                if (sum(NatZ_xmsyz(x,m,s))>0) nCs++;
            }
        }
    }
    
    os<<nCs<<tb<<"#number of sex x shell x maturity factor combinations to read in"<<endl;
    adstring_array factors(1,3);
    for (int x=1;x<=ANY_SEX;x++){
        factors(1) = tcsam::getSexType(x);
        for (int m=1;m<=ANY_MATURITY;m++){
            factors(2) = tcsam::getMaturityType(m);
            for (int s=1;s<=ANY_SHELL;s++){
                factors(3) = tcsam::getShellType(s);
                if (sum(NatZ_xmsyz(x,m,s))>0){ //only print out non-zero matrices
                    os<<"#-------"<<factors(1)<<cc<<factors(2)<<cc<<factors(3)<<endl;
                    os<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<endl;
                    os<<"#year    ss     "<<zBs<<endl;
                    os<<inpNatZ_xmsyc(x,m,s)<<endl;
                }
            }
        }
    }
    if (debug) cout<<"end SizeFrequencyData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void SizeFrequencyData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"SizeFrequencyData::writing to R"<<endl;
    adstring x = qt+tcsam::STR_FEMALE   +qt+cc+qt+tcsam::STR_MALE     +qt+cc+qt+tcsam::STR_ANY_SEX+qt;
    adstring m = qt+tcsam::STR_IMMATURE +qt+cc+qt+tcsam::STR_MATURE   +qt+cc+qt+tcsam::STR_ANY_MATURITY+qt;
    adstring s = qt+tcsam::STR_NEW_SHELL+qt+cc+qt+tcsam::STR_OLD_SHELL+qt+cc+qt+tcsam::STR_ANY_SHELL+qt;
    adstring y = wts::to_qcsv(yrs);
    adstring z = wts::to_qcsv(zBs);        
    for (int n=0;n<indent;n++) os<<tb;
        os<<"nAtZ=list(units="<<qt<<units<<qt<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc<<endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"years="; wts::writeToR(os,yrs); os<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"cutpts="; wts::writeToR(os,zCs); os<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"sample.sizes="; wts::writeToR(os,ss_xmsy,x,m,s,y); os<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<endl;
        wts::writeToR(os,NatZ_xmsyz,x,m,s,y,z); os<<endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) cout<<"SizeFrequencyData::done writing to R"<<endl;
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
    if (debug) cout<<"start BioData::read(...) "<<this<<endl;
    rpt::echo<<"#------------------------------------------"<<endl;
    rpt::echo<<"#BioData"<<endl;
    rpt::echo<<"#file name is "<<is.get_file_name()<<endl;
    rpt::echo<<"#------------------------------------------"<<endl;
    if (!is) {
        cout<<"Apparent error reading Bio Data."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_BIO_DATA)){
        cout<<"#Error reading bio data from "<<is.get_file_name()<<endl;
        cout<<"Expected keyowrd '"<<KW_BIO_DATA<<"' but got '"<<str<<"'"<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    rpt::echo<<str<<tb<<"#required keyword"<<endl;
    //SIZE BINS
    is>>nZBins;
    rpt::echo<<nZBins<<tb<<"#nZBins"<<endl;
    zBins.allocate(1,nZBins);
    is>>zBins;
    rpt::echo<<zBins<<tb<<"#zBins (mm CW)"<<endl;
    
    //WEIGHT-AT-SIZE
    {is>>unitsWatZ;
    rpt::echo<<unitsWatZ<<tb<<"#unitsWatZ"<<endl;
    double convToKG = tcsam::getConversionMultiplier(unitsWatZ,tcsam::UNITS_KG);
    rpt::echo<<"#using conversion factor from "<<unitsWatZ<<" to kg: "<<convToKG<<endl;
    wAtZ_xmz.allocate(FEMALE,MALE,IMMATURE,MATURE,1,nZBins);
    wAtZ_xmz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = tcsam::getSexType(factors(1));
        int mat   = tcsam::getMaturityType(factors(2));
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        if (sex||mat){
            is>>wAtZ_xmz(sex,mat);
            rpt::echo<<wAtZ_xmz(sex,mat)<<endl;
            wAtZ_xmz(sex,mat) *= convToKG;//convert to kg
            if (debug) {
                cout<<"#wAtZ_xmz("<<factors(1)<<cc<<factors(2)<<") [kg] ="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<wAtZ_xmz(sex,mat)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for wAtZ not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //PROBABILITY OF MATURING AT SIZE
    {prMature_xz.allocate(FEMALE,MALE,1,nZBins);
    prMature_xz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<tb<<"#factors"<<endl;
        int sex   = tcsam::getSexType(factors(1));
        if (sex){
            is>>prMature_xz(sex);
            rpt::echo<<prMature_xz(sex)<<endl;
            if (debug) {
                cout<<"#prMature_xz("<<factors(1)<<")="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<prMature_xz(sex)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for prMature_xz not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //FRACTION MATURE BY SEX, SHELL CONDITION
    {frMature_xsz.allocate(FEMALE,MALE,IMMATURE,MATURE,1,nZBins);
    frMature_xsz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        int sex   = tcsam::getSexType(factors(1));
        int shl   = tcsam::getShellType(factors(2));
        if (sex||shl){
            is>>frMature_xsz(sex,shl);
            rpt::echo<<frMature_xsz(sex,shl)<<endl;
            if (debug) {
                cout<<"#frMature_xsz("<<factors(1)<<cc<<factors(2)<<")="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<frMature_xsz(sex,shl)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for frMature_xsz not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //CV for mean size at min, max size bins
    {rpt::echo<<"#CV for min, max sizes"<<endl;
    cvMnMxZ_xc.allocate(FEMALE,MALE,1,2);
    cvMnMxZ_xc.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        rpt::echo<<factors(1)<<"#factors"<<endl;
        int sex   = tcsam::getSexType(factors(1));
        if (sex){
            is>>cvMnMxZ_xc(sex);
            rpt::echo<<cvMnMxZ_xc(sex)<<endl;
            if (debug) {
                cout<<"#cvMnMxZ_xc("<<factors(1)<<")="<<endl;
                cout<<cvMnMxZ_xc(sex)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for cvMnMxZ_xc not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //MIDPOINTS OF FISHERY SEASONS BY YEAR
    rpt::echo<<"#Timing of fisheries and mating"<<endl;
    is>>nyTiming;
    rpt::echo<<nyTiming<<tb<<"#number of years"<<endl;
    timing_yc.allocate(1,nyTiming,1,3);
    is>>timing_yc;
    rpt::echo<<"#timing"<<endl<<timing_yc<<endl;
    fshTiming_y.allocate(timing_yc(1)(1),timing_yc(nyTiming)(1));
    matTiming_y.allocate(timing_yc(1)(1),timing_yc(nyTiming)(1));
    for (int i=1;i<=nyTiming;i++){
        fshTiming_y(timing_yc(i,1)) = timing_yc(i,2);
        matTiming_y(timing_yc(i,1)) = timing_yc(i,3);
    }
    
    if (debug) cout<<"end BioData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void BioData::write(ostream & os){
    if (debug) cout<<"start BioData::write(...) "<<this<<endl;
    
    os<<KW_BIO_DATA<<endl;
    
    os<<"#-----------SIZE BINS---------------------------------#"<<endl;
    os<<nZBins<<tb<<"#number of size bins"<<endl;
    os<<"#size bins (mm CW)"<<endl<<zBins<<endl;
    
    {os<<"#-----------WEIGHT-AT-SIZE----------------------------#"<<endl;
    os<<unitsWatZ<<tb<<"#units for weight-at-size"<<endl;
    os<<MALE*MATURE<<tb<<"#number of factor combinations (sex x maturity state)"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int mat=IMMATURE;mat<=MATURE;mat++){
            factors(2) = tcsam::getMaturityType(mat);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<wAtZ_xmz(sex,mat)<<endl;
        }
    }}
    
    {os<<"#-----------PROBABILITY OF MATURING-------------------#"<<endl;
    os<<MALE<<tb<<"#number of factor combinations (sex)"<<endl;
    adstring_array factors(1,1);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = tcsam::getSexType(sex);
        os<<"#-------"<<factors(1)<<endl;
        os<<factors(1)<<endl;
        os<<prMature_xz(sex)<<endl;
    }}
    
    {os<<"#-----------FRACTION MATURE BY SEX, SHELL CONDITION---#"<<endl;
    os<<MALE*OLD_SHELL<<tb<<"#number of factor combinations (sex x shell condition)"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int shl=NEW_SHELL;shl<=OLD_SHELL;shl++){
            factors(2) = tcsam::getShellType(shl);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<frMature_xsz(sex,shl)<<endl;
        }
    }}
    
    {os<<"#-----------CV FOR MIN, MAX SIZES---------------------#"<<endl;
    os<<MALE<<tb<<"#number of factor combinations (sex)"<<endl;
    adstring_array factors(1,1);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = tcsam::getSexType(sex);
        os<<"#-------"<<factors(1)<<endl;
        os<<factors(1)<<endl;
        os<<cvMnMxZ_xc(sex)<<endl;
    }}
    
    {os<<"#-----------TIMING------------------------------------#"<<endl;
    os<<nyTiming<<tb<<"#number of years"<<endl;
    os<<"#year   midptFisheries  matingTime"<<endl;
    os<<timing_yc<<endl;}
    
    if (debug) cout<<"end BioData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void BioData::writeToR(ostream& os, string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<endl;
    indent++;
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"zBins=";wts::writeToR(os,zBins); os<<","<<endl;
        
        //weight-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"wAtZ=list(units="<<qt<<unitsWatZ<<qt<<cc<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            adstring x = qt+tcsam::STR_FEMALE+qt    +cc+ qt+tcsam::STR_MALE+qt;
            adstring s = qt+tcsam::STR_IMMATURE+qt  +cc+ qt+tcsam::STR_MATURE+qt;
            adstring c = wts::to_qcsv(zBins);
            wts::writeToR(os,wAtZ_xmz,x,s,c); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        
        //probability of maturing-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"prMat="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_FEMALE+qt    +cc+ qt+tcsam::STR_MALE+qt;
            adstring z = wts::to_qcsv(zBins);
            wts::writeToR(os,prMature_xz,x,z); os<<","<<endl;
        indent--;}
        
        //fraction mature-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"frMat="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_FEMALE+qt     +cc+ qt+tcsam::STR_MALE+qt;
            adstring s = qt+tcsam::STR_NEW_SHELL+qt  +cc+ qt+tcsam::STR_OLD_SHELL+qt;
            adstring z = wts::to_qcsv(zBins);
            wts::writeToR(os,frMature_xsz,x,s,z); os<<","<<endl;
        indent--;}
        
        //cv for min, max sizes
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"cvZs="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = qt+tcsam::STR_FEMALE+qt    +cc+ qt+tcsam::STR_MALE+qt;
            adstring c = "'minZ','maxZ'";
            wts::writeToR(os,cvMnMxZ_xc,x,c); os<<","<<endl;
        indent--;}
        
        //timing of fishery season midpoints and mating
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"timing="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring cols = "'midptFisheries','matingTime'";
            ivector yrs = wts::to_ivector(column(timing_yc,1));
            dmatrix tmp(1,2,1,nyTiming);
            tmp(1) = column(timing_yc,2);
            tmp(2) = column(timing_yc,3);
            wts::writeToR(os,trans(tmp),yrs,cols); os<<endl;
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
    rpt::echo<<"#------------------------------------------"<<endl;
    rpt::echo<<"reading Model Datasets file."<<endl;
    rpt::echo<<"#file name is '"<<is.get_file_name()<<"'"<<endl;
    rpt::echo<<"#------------------------------------------"<<endl;
    is>>fnBioData;
    rpt::echo<<fnBioData<<tb<<"#bio data filename"<<endl;
    is>>nFsh;
    rpt::echo<<nFsh<<tb<<"#number of fishery datasets to read in"<<endl;
    if (nFsh){
        fnsFisheryData.allocate(1,nFsh);
        is>>fnsFisheryData; 
        for (int i=1;i<=nFsh;i++) rpt::echo<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<endl;
    }
    is>>nSrv;
    rpt::echo<<nSrv<<tb<<"#number of survey datasets to read in"<<endl;
    if (nSrv){
        fnsSurveyData.allocate(1,nSrv);
        is>>fnsSurveyData; 
        for (int i=1;i<=nSrv;i++) rpt::echo<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<endl;
    }
    
    //          Bio data
    {ptrBio = new BioData();
    cifstream strm(fnBioData,ios::in);
    rpt::echo<<"#------Bio data ------------"<<endl;
    strm>>(*ptrBio);
    }  
    
    //          Fishery data 
    if (nFsh) {
        rpt::echo<<"#-------Fishery Datasets---------"<<endl;
        ppFsh = new FisheryData*[nFsh];
        for (int i=0;i<nFsh;i++) {
            ppFsh[i] = new FisheryData();
            cifstream strm(fnsFisheryData(i+1),ios::in);
            rpt::echo<<endl<<"#----------Fishery Data "<<i+1<<"-----"<<endl;
            strm>>(*ppFsh[i]);
        }
    }
    
    //          Survey data
    if (nSrv) {
        rpt::echo<<"#-------Survey Datasets---------"<<endl;
        ppSrv = new SurveyData*[nSrv];
        for (int i=0;i<nSrv;i++) {
            ppSrv[i] = new SurveyData();
            cifstream strm(fnsSurveyData(i+1),ios::in);
            rpt::echo<<endl<<"#----------Survey Data "<<i+1<<"-----"<<endl;
            strm>>(*ppSrv[i]);
        }
    }
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::write(ostream & os){
    if (debug) cout<<"start ModelDatasets::write(...) "<<this<<endl;
    os<<fnBioData<<tb<<"#tanner crab biological data file"<<endl;
    os<<"#-------fishery data files---------"<<endl;
    os<<nFsh<<tb<<"#number of fishery data files"<<endl;
    for (int i=1;i<+nFsh;i++) os<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<endl;
    os<<"#-------survey data files---------"<<endl;
    os<<nSrv<<tb<<"#number of survey data files"<<endl;
    for (int i=1;i<+nSrv;i++) os<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<endl;
    os<<"#-----biological data---------"<<endl;
    os<<(*ptrBio)<<endl;
    os<<"#-------fishery data ---------"<<endl;
    if (nFsh){
        for (int i=1;i<nFsh;i++) os<<"#----fishery dataset "<<i<<endl<<(*ppFsh[i-1])<<endl;
        os<<"#----fishery dataset "<<nFsh<<endl<<(*ppFsh[nFsh-1])<<endl;
    }
    os<<"#-------survey data ---------"<<endl;
    if (nSrv){
        for (int i=1;i<nSrv;i++) os<<"#----survey dataset "<<i<<endl<<(*ppSrv[i-1])<<endl;
        os<<"#----survey dataset "<<nSrv<<endl<<(*ppSrv[nSrv-1]);
    }
   if (debug) cout<<"end ModelDatasets::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelDatasets::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list("<<endl;
    indent++;
        //bio data as list
            ptrBio->writeToR(os,"bio",indent); os<<cc<<endl;
        
        //survey data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"surveys=list("<<endl;
        indent++;
            if (ppSrv) {
                for (int i=0;i<(nSrv-1);i++) {
                    ppSrv[i]->writeToR(os,(char*)ppSrv[i]->name,indent); os<<cc<<endl;
                }
                ppSrv[nSrv-1]->writeToR(os,(char*)ppSrv[nSrv-1]->name,indent); os<<endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<endl;            
        
        //fishery data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"fisheries=list("<<endl;
        indent++;
            if (ppFsh) {
                for (int i=0;i<(nFsh-1);i++) {
                    ppFsh[i]->writeToR(os,(char*)ppFsh[i]->name,indent); os<<cc<<endl;
                }
                ppFsh[nFsh-1]->writeToR(os,(char*)ppFsh[nFsh-1]->name,indent); os<<endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<")"<<endl;            
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end ModelDatasets/////////////////////////
