/* 
 * File:   OFLCalcs.hpp
 * Author: WilliamStockhausen
 *
 * Created on March 4, 2016, 10:16 AM
 */

#ifndef OFLCALCS_HPP
#define	OFLCALCS_HPP

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
        dvector  R_z;      //recruitment size distribution
        dmatrix  w_mz;     //weight-at-size
        d3_array M_msz;    //natural mortality
        dmatrix  Th_sz;    //pr(molt to maturity|pre-molt size, molt)
        d3_array T_szz;    //growth matrices (indep. of molt to maturity)
        d4_array capF_fmsz;//fishery capture rates
        d4_array retF_fmsz;//retention function
        dvector  hm_f;     //discard mortality
        
    public:
        double XX;
        double avR;
        double B100;
        double Fmsy;
        double Bmsy;
        d3_array ct_msz; //catch mortality
        d4_array rm_fmsz;//retained mortality
        d4_array dm_fmsz;//discard mortality
        
    public:
        Tier3_Calculator(int npZBs);
        ~Tier3_Calculator(){}
        double calcB100(double R);
        double calcBmsy(double R);
        double calcFmsy(double R);
        
    public:
        d3_array calcEqNatZ(dvector& R_z,d3_array& S1_msz, 
                            dmatrix& Th_sz, d3_array& T_szz, d3_array& S2_msz);
        d3_array calcEqNatZF0(double R);
        d3_array calcEqNatZFM(double R);
        double calcEqMMBatF(double R, double dirF);
        double calcMMB(d3_array& n_msz);
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
        dmatrix  w_mz;     //weight-at-size
        d3_array M_msz;    //natural mortality
        dmatrix  Th_sz;    //pr(molt to maturity|pre-molt size, molt)
        d3_array T_szz;    //growth matrices (indep. of molt to maturity)
        d4_array capF_fmsz;//fishery capture rates
        d4_array retF_fmsz;//retention function
        dvector  hm_f;     //discard mortality
        
    public:
        double spB;       //mature biomass at time of mating
        d3_array ct_msz;  //total catch
        d4_array rm_fmsz; //retained catch mortality
        d4_array dm_fmsz; //discard catch mortality
    
    public:
        PopProjector(int npZBs);
        ~PopProjector(){}
        /**
         * Project sex-specific population abundance forward one year,
         * based on sex-specific population abundance on July 1.
         * 
         * @param n_msz - initial abundance
         * 
         * @return final sex-specific abundance
         */
        d3_array project(d3_array& n_msz);
        /**
         * Calculate mature biomass-at-mating based on sex-specific population
         * abundance on July 1.
         * 
         * @param n_msz - initial abundance
         * 
         * @return mature biomass (units??)
         */
        double calcMatBio(d3_array& n_msz);
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
         * @return - final abundance after fisheries
         */
        d3_array applyFM(d3_array& n_msz);
        /**
         * Project sex-specific population abundance forward through 
         * molting/growth.
         * 
         * @param n_msz - initial abundance
         * @return - final abundance after molting/growth
         */
        d3_array applyMG(d3_array& n_msz);
        /**
         * Calculate total biomass from population abundance
         * @param n_msz
         * @return biomass (units??)
         */
        double calcBiomass(d3_array& n_msz);
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
        double calcOFL(void);
        double calcHCR(double currMMB, double Bmsy, double Fmsy);
};

#endif	/* OFLCALCS_HPP */

