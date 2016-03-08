/* 
 * File:   OFLCalcs.hpp
 * Author: WilliamStockhausen
 *
 * Created on March 4, 2016, 10:16 AM
 */

#ifndef OFLCALCS_HPP
#define	OFLCALCS_HPP

/**
 * Class to facilitate population dynamics calculations for one sex.
 */
class PopDyInfo {
    public:
        static int debug;
        static ostream& cout;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        
    public:
        dvector  R_z;      //recruitment size distribution
        dmatrix  w_mz;     //weight-at-size
        d3_array M_msz;    //natural mortality
        dmatrix  Th_sz;    //pr(molt to maturity|pre-molt size, molt)
        d3_array T_szz;    //growth matrices (indep. of molt to maturity)
        
    public:
        PopDyInfo(int npZBs);
        ~PopDyInfo(){}
        /**
         * Calculate single-sex mature biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * 
         * @return mature biomass (units??)
         */
        double calcMatureBiomass(d3_array& n_msz);
        /**
         * Calculate total single-sex biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * 
         * @return total biomass (units??)
         */
        double calcTotalBiomass(d3_array& n_msz);
        /**
         * Calculate single-sex survival probabilities, given natural mortality rates,
         * over a given period of time (dt)
         * 
         * @param dt - period of time (in years)
         * 
         * @return d3_array S_msz
         */
        d3_array calcSurvival(double dt);
        d3_array applyNM(double dt, d3_array& n_msz);
        d3_array applyMG(d3_array& n_msz);
};//PopDyInfo

/**
 * Class to facilitate fishery catch calculations for one sex.
 */
class CatchInfo {
    public:
        static int debug;
        static ostream& cout;
        int nMSs;
        int nSCs;
        int nZBs;
        int nFsh;   //number of fisheries
        
    public:
        double maxF; //maximum capture rate in target fishery
        d4_array capF_fmsz;//fishery capture rates
        d4_array retF_fmsz;//retention function
        dvector  hm_f;     //discard mortality rates
        
    public:
        d3_array ct_msz; //catch mortality (abundance)
        d4_array rm_fmsz;//retained mortality (abundance)
        d4_array dm_fmsz;//discard mortality (abundance)
        
    public:
        CatchInfo(int npZBs, int npFsh);
        ~CatchInfo(){}
        /**
         * Calculate maximum fishing capture rate for target fishery (f=1)
         * 
         * Modifies: maxF
         * 
         * @return maxF
         */
        double findMaxTargetCaptureRate(void);
        /**
         * Calculate catch abundance (ct_msz, rm_fmsz, dm_fmsz) and post-fisheries
         * abundance based on target fishery capture rate 'dirF' and initial population 
         * abundance n_msz.
         * 
         * Modifies: ct_msz, rm_fmsz, dm_fmsz
         * 
         * @param dirF - fully-selected directed fishery capture rate
         * @param n_msz - pre-fisheries population size
         * 
         * @return np_msz - post-fisheries population size (d3_array)
         * 
         */
        d3_array applyFM(double dirF, d3_array& n_msz);
        /**
         * Calculates probabilities of surviving fisheries, given directed 
         * fishing capture rate 'dirF'.
         * 
         * @param dirF - fully-selected capture rate for target fishery (f=1)
         * 
         * @return d3_array of survival probabilities S_msz.
         */
        d3_array calcSurvival(double dirF);
        /**
         * Set sex-specific fishery capture rates.
         * 
         * @param x - sex to select
         * @param capF_fxmsz - input capture rates
         */
        void setCaptureRates(int x, d5_array& capF_fxmsz);
        /**
         * Set sex-specific retention functions.
         * 
         * @param x - sex to select
         * @param retF_fxmsz - input retention functions
         */
        void setRetentionFcns(int x, d5_array& retF_fxmsz);
        /**
         * Set handling mortality rates, by fishery.
         * 
         * @param pHM_f - dvector of handling mortality rates
         */
        void setHandlingMortality(dvector& pHM_f);
};//CatchInfo

class Tier3_Calculator {
    public:
        static int debug;
        static ostream& cout;
        int nMSs;
        int nSCs;
        int nZBs;
        int nFsh;   //number of fisheries
        double dtF; //time to fisheries
        double dtM; //time to molting/growth
        
    public:
        PopDyInfo* pPI;//pointer to PopDyInfo object for males
        CatchInfo* pCI;//pointer to CatchInfo object for males
    public:
        double XX;
        double avR;
        double B100;
        double Fmsy;
        double Bmsy;
        
    public:
        Tier3_Calculator(int npZBs, int npFsh);
        ~Tier3_Calculator(){delete pPI; delete pCI;}
        double calcB100(double R);
        double calcBmsy(double R);
        double calcFmsy(double R);
        
    public:
        d3_array calcEqNatZ(dvector& R_z,d3_array& S1_msz, 
                            dmatrix& Th_sz, d3_array& T_szz, d3_array& S2_msz);
        d3_array calcEqNatZF0(double R);
        d3_array calcEqNatZFM(double R, double dirF);
        double   calcEqMMBatF(double R, double dirF);
};

/**
 * Class to enable sex-specific population projections.
 */
class PopProjector{
    public:
        static int debug;
        static ostream& cout;
        int nMSs;
        int nSCs;
        int nZBs;
        int nFsh;   //number of fisheries
        double dtF; //time to fisheries
        double dtM; //time to molting/growth
        
    public:
        double spB;       //mature biomass at time of mating
    
    public:
        PopDyInfo* pPI;//pointer to single sex PopDyInfo object
        CatchInfo* pCI;//pointer to single sex CatchInfo object
    
    public:
        PopProjector(int npZBs, int npFsh);
        ~PopProjector(){delete pPI; delete pCI;}
        /**
         * Project sex-specific population abundance forward one year,
         * based on sex-specific population abundance on July 1.
         * Also calculates:
         *      spB - spawning biomass at mating time
         *      ct_msz - total fishing mortality (abundance)
         *      rm_fmsz - retained mortality, by fishery
         *      dm_fmsz - discard mortality, by fishery
         * 
         * @param n_msz - initial abundance
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * 
         * @return final sex-specific abundance
         */
        d3_array project(double dirF, d3_array& n_msz);
        /**
         * Calculate mature biomass based on single-sex population
         * abundance on July 1.
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * 
         * @return mature biomass (units??)
         */
        double calcMatureBiomass(double dirF, d3_array& n_msz);
        /**
         * Project sex-specific population abundance forward by time interval "dt", 
         * applying only natural mortality.
         * 
         * @param dt - time interval
         * @param n_msz - initial abundance
         * @return - final abundance
         */
        d3_array applyNM(double dt, d3_array& n_msz);
        /**
         * Project sex-specific population abundance forward through 
         * pulse fisheries.
         * 
         * @param n_msz - initial abundance
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @return - final abundance after fisheries
         */
        d3_array applyFM(double dirF, d3_array& n_msz);
        /**
         * Project sex-specific population abundance forward through 
         * molting/growth.
         * 
         * @param n_msz - initial abundance
         * @return - final abundance after molting/growth
         */
        d3_array applyMG(d3_array& n_msz);
};

class OFL_Calculator{
    public:
        static int debug;
        static ostream& cout;
        int nSXs;
        int nMSs;
        int nSCs;
        int nZBs;
        int nFsh;   //number of fisheries

    public:
        int tier;
        double alpha; 
        double beta;
        dmatrix ofl_fx;
        
    public:
        Tier3_Calculator* pT3C;//pointer to a Tier3_Calculator
        PopProjector*     pPrjM;//pointer to PopProjector for males
        PopProjector*     pPrjF;//pointer to PopProjector for females
    public:
        OFL_Calculator(int Tier, int npZBs);
        ~OFL_Calculator(){delete pT3C;delete pPrjM;delete pPrjF;}
        double calcBmsy(double R);
        double calcFmsy(double R);
        double calcFofl(double Bmsy, double Fmsy, d3_array& n_msz);
        double calcOFL(double Fofl, d4_array& n_xmsz);
        double calcHCR(double currMMB, double Bmsy, double Fmsy);
};

#endif	/* OFLCALCS_HPP */

