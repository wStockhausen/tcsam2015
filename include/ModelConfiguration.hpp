/* 
 * File:   ModelConfiguration.hpp
 * Author: william.stockhausen
 *
 * Created on March 12, 2013, 8:59 AM
 * 
 * Includes:
 *  class ModelConfiguration
 *  class ModelOptions
 * 
 * History:
 * 2014-10-30: 1. Removed nSXs, nSCs, nMSs as static variables that were set through
 *                the model configuration file. These should now be changed in ModelConstants.hpp
 *                and the model recompiled for different configurations.
 * 2014-12-03: 1. Added asYr (assessment year) as input, replacing mxYr. Now, mxYr = asYr-1.
 * 2015-03-01: 1. Added "dims" adstring variables to facilitate output to R.
 *             2. Added optsGrowth and optsInitNatZ options to ModelOptions.
 * 2015-05-12: 1. Added cvFDevsPen, phsDecrFDevsPen,phsZeroFDevsPen for F-devs penalties to ModelOptions
 * 2015-05-13: 1. Added phsLastDevPen, wgtLastDevPen to ModelOptions
 * 2015-05-26: 1. Added penWgtSmthLgtPrMat, penWgtNonDecLgtPrMat to ModelOptions
 */

#ifndef MODELCONFIGURATION_HPP
    #define	MODELCONFIGURATION_HPP

//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
    class ModelConfiguration {
    public:
        static int debug;  //flag to print debug info
        
        static int nZBs;//number of size bins 
        static int mnYr;//min model year
        static int asYr;//assessment year (final pop numbers given for July 1, asYr)
        static int mxYr;//max model year (mxYr = asYr-1)
        static int nFsh;//number of fisheries 
        static int nSrv;//number of surveys 
        
        static int    jitter;  //flag to jitter initial parameter values
        static double jitFrac; //fraction to jitter bounded parameter values
        static int    resample;//flag to resample initial parameter values
        static double vif;     //variance inflation factor for resampling parameter values
    public:
        adstring cfgName;//model configuration name

        int runOpMod;    //flag to run operating model only
        int fitToPriors; //flag to fit to priors
        
        dvector zMidPts;     //size bin midpoints (CW in mm)
        dvector zCutPts;     //size bin cutpoints (CW in mm)
        dvector onesZMidPts; //vector of 1's same size as zMidPts

        adstring fnMDS; //model datasets file name
        adstring fnMPI; //model parameters info file name
        adstring fnMOs; //model options file name
        
        adstring_array lblsFsh;//labels for fisheries
        adstring_array lblsSrv;//labels for surveys
        
        adstring csvYrs;  //csv string of model years (mnYr:mxYr)
        adstring csvYrsP1;//csv string of model years+1 (mnYr:mxYr+1)
        adstring csvSXs;//csv string of sexes
        adstring csvMSs;//csv string of maturity states
        adstring csvSCs;//csv string of shell conditions
        adstring csvZCs;//csv string of size bin cutpoints
        adstring csvZBs;//csv string of size bin midpoints
        adstring csvFsh;//csv string of fishery labels
        adstring csvSrv;//csv string of survey labels

        adstring dimYrsToR;  //R dim string of model years (mnYr:mxYr)
        adstring dimYrsP1ToR;//R dim string of model years+1 (mnYr:mxYr+1)
        adstring dimSXsToR;//R dim string of sexes
        adstring dimMSsToR;//R dim string of maturity states
        adstring dimSCsToR;//R dim string of shell conditions
        adstring dimZCsToR;//R dim string of size bin cutpoints
        adstring dimZBsToR;//R dim string of size bin midpoints
        adstring dimZPsToR;//R dim string of size bin midpoints (alternative)
        adstring dimFshToR;//R dim string of fishery labels
        adstring dimSrvToR;//R dim string of survey labels

    public:
        ModelConfiguration();
        ~ModelConfiguration();
        ModelConfiguration& operator =(const ModelConfiguration & rhs);
        
        /**
         * Set max model year (for retrospective model runs).
         * 
         * @param yr - new max model year
         */
        void setMaxModelYear(int yr);
        /**
         * Tests if mnYr <= yr <= mxYr. 
         * 
         * @param yr
         * @return 1 if true, 0 if false
         */
        int isModelYear(int yr){if ((mnYr<=yr)&&(yr<=mxYr)) return 1; return 0;}
        void read(const adstring & fn);   //read file in ADMB format
        void write(const adstring & fn);  //write object to file in ADMB format
        void read(cifstream & is);        //read file in ADMB format
        void write(std::ostream & os);         //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream&    operator >>(cifstream & is, ModelConfiguration & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os,   ModelConfiguration & obj){obj.write(os);;return os;}
    };

//--------------------------------------------------------------------------------
//          ModelOptions
//--------------------------------------------------------------------------------
    class ModelOptions {
    public:
        static int debug;  //flag to print debug info        
    public:
        ModelConfiguration* ptrMC;      //pointer to model configuration object
        adstring_array lblsFcAvgOpts;   //labels for capture rate averaging options
        ivector optsFcAvg;              //selected options for averaging capture rate
        adstring_array lblsGrowthOpts;  //labels for growth options
        int optsGrowth;                 //selected option for growth calculations
        adstring_array lblsInitNatZOpts;//labels for initial n-at-z options
        int optsInitNatZ;               //selected option for initial n-at-z calculations
        double cvFDevsPen;              //penalty for F-devs
        int phsDecrFDevsPen;            //phase to start decreasing fpenCV
        int phsZeroFDevsPen;            //phase at which to turn off penalties on F-devs
        double wgtLastDevsPen;          //penalty for last dev in each devs vector
        int phsLastDevsPen;             //phase to start the penalty on the last devs
        double wgtSmthLgtPrMat;      //penalty on maturity ogive smoothness
        double wgtNonDecLgtPrMat;    //penalty on maturity ogive non-decreasing

    public:
        ModelOptions(ModelConfiguration& mc);
        ~ModelOptions(){}
        
        void read(cifstream & is);        //read file in ADMB format
        void write(std::ostream & os);         //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream&    operator >>(cifstream & is, ModelOptions & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os,   ModelOptions & obj){obj.write(os);;return os;}
    };
    
#endif	/* MODELCONFIGURATION_HPP */

