/* 
 * File:   ModelData.hpp
 * Author: william.stockhausen
 *
 * Created on 2014-02-11.
 * 
 * 20140930: 1. Added recLag to BioData.
 */

#ifndef MODELDATA_HPP
#define MODELDATA_HPP

//**********************************************************************
//  Includes
//      AggregateCatchData
//      EffortData
//      SizeFrequencyData
//      BioData
//      SurveyData
//      RetainedCatchData
//      ObservedCatchData
//      FisehryData
//      ModelDatasets
//**********************************************************************
class ModelConfiguration; //forward definition
class IndexRange;

//--------------------------------------------------------------------------------
//          AggregateCatchData
//--------------------------------------------------------------------------------
    class AggregateCatchData {
    public:
        static int debug;
        const static adstring KW_ABUNDANCE_DATA;//keyword indicating abundance (numbers) data
        const static adstring KW_BIOMASS_DATA;  //keyword indicating biomass (weight) data
    protected:
        dmatrix inpC_yc; //input aggregate catch data (year, female abundance, female cv, male abundance, male cv)
    public:
        adstring type;  //type (abundance, biomass) of data (matches one of the keywords above)
        int optFit;     //objective function fitting option
        int llType;     //likelihood function type
        dvector nlls;   //negative log-likelihood components
        int ny;         //number of years of aggregate catch
        adstring units; //units for aggregate catch data
        ivector yrs;    //years for aggregate catch data
        dmatrix C_xy;   //aggregate catch by sex, year (converted from units to THOUSANDS of crab or MT))
        dmatrix cv_xy;  //aggregate catch cv's by sex, year
        dmatrix sd_xy;  //aggregate catch stdv's by sex, year
        
    public:
        AggregateCatchData(){}
        ~AggregateCatchData(){}
        /**
         * Replace catch data C_xy with new data. 
         * Also modifies inpC_yc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param dmatrix newC_yx (note index order)
         */
        void replaceCatchData(int iSeed,random_number_generator& rng,dmatrix& newC_yx);
        /**
         * Save the negative log-likelihoods from a model fit (values only).
         * 
         * @param nlls
         */
        void saveNLLs(dvar_vector& nlls);
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, AggregateCatchData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   AggregateCatchData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          EffortData
//--------------------------------------------------------------------------------
    class EffortData {
    public:
        static int debug;
        const static adstring KW_EFFORT_DATA;//keyword indicating effort data
    public:
        int ny;               //number of years of effort data
        IndexRange* ptrAvgIR; //interval to average effort/fishing mortality
        adstring units;       //units for potlifts
        dmatrix  inpEff_yc;   //input effort data (year)x(year,potlifts)
        dvector yrs;          //years w/ effort data
        dvector eff_y;        //effort data vector
    public:
        EffortData(){ptrAvgIR=0;}
        ~EffortData();
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, EffortData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   EffortData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          SizeFrequencyData
//--------------------------------------------------------------------------------
    class SizeFrequencyData {
    public:
        static int debug;
        const static adstring KW_SIZEFREQUENCY_DATA;//keyword indicating size frequency data
    private:
        wts::adstring_matrix factors;//factor combinations for input numbers-at-size
        dmatrix inpSS_yc;             //input sample size matrix TODO: based on num tows??
        d5_array inpNatZ_xmsyc;       //input numbers-at-size data (sex,maturity state,shell condition,year,year+sample_size+nAtZ)
    public:
        int optFit; //objective function fitting option
        int llType; //likelihood function type
        
        int nZCs;    //number of size bin cut pts
        dvector zCs; //cut points for size bins
        dvector zBs; //size bins
        
        int ny;             //number of years of size frequency data
        adstring units;     //units for numbers-at-size data
        ivector  yrs;       //years of size frequency data
        d4_array ss_xmsy;   //sample sizes for size frequency data
        d5_array NatZ_xmsyz;//raw size frequency data
        d5_array PatZ_xmsyz;//normalized size frequency data (sums to 1 over xmsz for each y)
    public:
        SizeFrequencyData(){}
        ~SizeFrequencyData(){}
        /**
         * Replace catch-at-size data NatZ_xmsyz with new data. 
         * Also modifies inpNatZ_xmsyc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param d5_array newNatZ_yxmsz
         */
        void replaceSizeFrequencyData(int iSeed,random_number_generator& rng,d5_array& newNatZ_xmsyz);
        /**
         * Save the negative log-likelihoods from a model fit (values only).
         * 
         * @param nlls
         */
        void saveNLLs(dvar_matrix& nlls);
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, SizeFrequencyData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   SizeFrequencyData & obj){obj.write(os); return os;}
    public:
        void normalize(void);
    };

//--------------------------------------------------------------------------------
//          BioData
//--------------------------------------------------------------------------------
    class BioData {
    public:
        static int debug;
        const static adstring KW_BIO_DATA;
    public:
        int nZBins;           //number of size bin cut pts
        dvector zBins;        //size bins
        adstring unitsWatZ;   //units for weight-at-size
        d3_array wAtZ_xmz;    //weight at size (in kg)
        dmatrix prMature_xz;  //probability of maturity at size by sex
        d3_array frMature_xsz;//fraction mature at size by sex, shell condition
        dmatrix cvMnMxZ_xc;   //cv of min size, max size by sex
        
        int recLag;         //recruitment lag (in years)
        double fshTimingTypical;//typical midpoint of fishery seasons
        double matTimingTypical;//typical timing of mating
        int nyAtypicalT;        //number of years for atypical fishery season midpoints, time-at-mating
        dvector fshTiming_y;    //timing of midpoint of fisheries seasons
        dvector matTiming_y;    //timing of mating
    protected:
        dmatrix timing_yc;//midpoint of fishery season, time-at-mating by year
    public:
        BioData(){}
        ~BioData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, BioData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   BioData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          CatchData
//--------------------------------------------------------------------------------
    class CatchData {
    public:
        static int debug;
        const static adstring KW_CATCH_DATA;
    protected:
        dmatrix inpN_yc;       //input catch abundance data (year,female abundance,female cv,male abundance, male cv, total abundance, cv_total)
        dmatrix inpB_yc;       //input catch biomass   data (year,female biomass,female cv,male biomass, male cv, total biomass, cv_total)
        //TODO: need sample sizes (tows/potlifts, non-zero tows/potlifts, num individ.s by xms, calculated ss)
        d5_array inpNatZ_xmsyc;//input numbers-at-size data (sex,maturity,shell,year) x (year,sample_size,nAtZ)
    public:
        adstring type;    //data type (e.g. survey or fishery)
        adstring name;    //data source name
        int hasN;         //abundance data flag
        AggregateCatchData* ptrN;//pointer to aggregate abundance data
        int hasB;         //biomass data flag
        AggregateCatchData* ptrB;//pointer to aggregate biomass data
        int hasZFD;       //numbers-at-size data flag
        SizeFrequencyData* ptrZFD;//pointer to numbers-at-size data
        
    public:
        CatchData();
        ~CatchData();
        /**
         * Replaces catch data based on newNatZ_yxmsz.
         * 
         * @param newNatZ_yxmsz - catch-at-size by sex/maturity/shell condition/year
         * @param wAtZ_xmz - weight-at-size by sex/maturity
         */
        virtual void replaceCatchData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz, d3_array& wAtZ_xmz);
        virtual void read(cifstream & is);//read file in ADMB format
        virtual void write(ostream & os); //write object to file in ADMB format
        virtual void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, CatchData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   CatchData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          SurveyData
//--------------------------------------------------------------------------------
class SurveyData: public CatchData {
    public:
        static int debug;
        const static adstring KW_SURVEY_DATA;
    public:
        SurveyData():CatchData(){type="survey";}
        ~SurveyData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
    };

//--------------------------------------------------------------------------------
//          FisheryData
//--------------------------------------------------------------------------------
    class FisheryData {
    public:
        static int debug;
        const static adstring KW_FISHERY_DATA;
    public:
        adstring name;     //fishery name
        int hasEff;        //flag indicating effort data
        EffortData* ptrEff;//pointer to effort data
        int hasRCD;        //flag indicating retained catch data
        CatchData* ptrRCD; //pointer to CatchData object for retained catch
        int hasTCD;        //flag indicating observed total catch data
        CatchData* ptrTCD; //pointer to CatchData object for observed total catch
        int hasDCD;        //flag indicating observed discard catch data
        CatchData* ptrDCD; //pointer to CatchData object for observed discards
    public:
        FisheryData();
        ~FisheryData();
        /**
         * Replace existing catch data with new values.
         * 
         * @param newCatZ_yxmsz - new total catch-at-size
         * @param newRatZ_yxmsz - new retained catch-at-size
         * @param wAtZ_xmz - weight-at-size
         */
        void replaceCatchData(int iSeed,random_number_generator& rng,d5_array& newCatZ_yxmsz,d5_array& newRatZ_yxmsz,d3_array& wAtZ_xmz);
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, FisheryData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   FisheryData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//         ModelDatasets
//--------------------------------------------------------------------------------
    class ModelDatasets {
    public:
        static int debug;
    public:
        ModelConfiguration* pMC;//pointer to ModelConfiguration object
        adstring fnBioData;//bio data file name
        BioData* ptrBio;   //pointer to bio dataset object
        
        int nFsh;//number of fishery datasets to read
        adstring_array fnsFisheryData;//fishery data file names       
        FisheryData**  ppFsh;         //pointer to array of pointers to fishery dataset objects
        
        int nSrv;//number of survey datasets to read
        adstring_array fnsSurveyData;//survey data files names
        SurveyData**   ppSrv;        //pointer to array of pointers to survey dataset objects
    public:
        ModelDatasets(ModelConfiguration* ptrMC);
        ~ModelDatasets();
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, ModelDatasets & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   ModelDatasets & obj){obj.write(os); return os;}
    };
#endif  /* MODELDATA_HPP */

