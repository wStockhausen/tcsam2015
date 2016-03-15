    #include <math.h>
    #include <time.h>
    #include <limits>
    #include <admodel.h>
    #include "TCSAM.hpp"
    adstring model  = "TCSAM2015";
    adstring modVer = "2016.03.15"; 
    
    time_t start,finish;
    
    //model objects
    ModelConfiguration*  ptrMC; //ptr to model configuration object
    ModelParametersInfo* ptrMPI;//ptr to model parameters info object
    ModelOptions*        ptrMOs;//ptr to model options object
    ModelDatasets*       ptrMDS;//ptr to model datasets object
    ModelDatasets*       ptrSimMDS;//ptr to simulated model datasets object
    OFLResults*          ptrOFL;   //pointer to OFL results object for MCMC calculations
        
    //dimensions for R output
    adstring yDms;
    adstring xDms;
    adstring mDms;
    adstring sDms;
    adstring fDms;
    adstring vDms;
    adstring ypDms;
    adstring zbDms;
    adstring zpDms;
    adstring zcDms;
    
    //file streams and filenames
    std::ofstream mcmc;        //stream for mcmc output
    
    //filenames
    adstring fnMCMC = "TCSAM2015.MCMC.R";
    adstring fnConfigFile;//configuration file
    adstring fnResultFile;//results file name 
    adstring fnPin;       //pin file
    
    //runtime flags (0=false)
    int jitter     = 0;//use jittering for initial parameter values
    int resample   = 0;//use resampling for initial parameter values
    int opModMode  = 0;//run as operating model, no fitting
    int usePin     = 0;//flag to initialize parameter values using a pin file
    int doRetro    = 0;//flag to facilitate a retrospective model run
    int fitSimData = 0;//flag to fit model to simulated data calculated in the PRELIMINARY_CALCs section
    
    int yRetro = 0; //number of years to decrement for retrospective model run
    int iSeed =  0;//default random number generator seed
    random_number_generator rng(-1);//random number generator
    int iSimDataSeed = 0;
    random_number_generator rngSimData(-1);//random number generator for data simulation
    
    //debug flags
    int debugModelConfig     = 0;
    int debugModelDatasets   = 0;
    int debugModelParamsInfo = 0;
    int debugModelParams     = 0;
    int debugModelOptions    = 0;
    
    int debugDATA_SECTION    = 0;
    int debugPARAMS_SECTION  = 0;
    int debugPRELIM_CALCS    = 0;
    int debugPROC_SECTION    = 0;
    int debugREPORT_SECTION  = 0;
    
    int showActiveParams = 0;    
    int debugRunModel    = 0;    
    int debugObjFun      = 0;
    
    int debugMCMC = 0;
    
    //note: consider using std::bitset to implement debug functionality
    int dbgOFL = 0;
    int dbgCalcProcs = 10;
    int dbgObjFun = 20;
    int dbgNLLs   = 25;
    int dbgPriors = tcsam::dbgPriors;
    int dbgPopDy  = 70;
    int dbgApply  = 80;
    int dbgDevs   = 90;
    int dbgAll    = tcsam::dbgAll;
    
    int nSXs    = tcsam::nSXs;
    int MALE    = tcsam::MALE;
    int FEMALE  = tcsam::FEMALE;
    int ALL_SXs = tcsam::ALL_SXs;
    
    int nMSs     = tcsam::nMSs;
    int IMMATURE = tcsam::IMMATURE;
    int MATURE   = tcsam::MATURE;
    int ALL_MSs  = tcsam::ALL_MSs;
    
    int nSCs      = tcsam::nSCs;
    int NEW_SHELL = tcsam::NEW_SHELL;
    int OLD_SHELL = tcsam::OLD_SHELL;
    int ALL_SCs   = tcsam::ALL_SCs;
    
    double smlVal = 0.00001;//small value to keep things > 0
    
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <TCSAM2015.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
    rpt::echo<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    rpt::echo<<"#Starting DATA_SECTION"<<endl;
    cout<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    cout<<"#Starting DATA_SECTION"<<endl;
    int on = 0;
    int flg = 0;
    rpt::echo<<"#------Reading command line options---------"<<endl;
    //configFile
    fnConfigFile = "TCSAM2015_ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //resultsFile
    fnResultFile = "TCSAM2015_ResultFile";//default results file name
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resultsFile"))>-1) {
        fnResultFile = ad_comm::argv[on+1];
        rpt::echo<<"#results file changed to '"<<fnResultFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //opModMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-opModMode"))>-1) {
        opModMode=1;
        rpt::echo<<"#operating model mode turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        rpt::echo<<"#debugModelDatasets turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParamsInfo
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParamsInfo"))>-1) {
        debugModelParamsInfo=1;
        rpt::echo<<"#debugModelParamsInfo turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelOptions
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelOptions"))>-1) {
        debugModelOptions=1;
        rpt::echo<<"#debugModelOptions turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //doRetro
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-doRetro"))>-1) {
        doRetro=1;
        cout<<"#doRetro turned ON"<<endl;
        rpt::echo<<"#doRetro turned ON"<<endl;
        if (on+1<argc) {
            yRetro=atoi(ad_comm::argv[on+1]);
            cout<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
            rpt::echo<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //fitSimData
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-fitSimData"))>-1) {
        fitSimData=1;
        if (on+1<argc) {
            iSimDataSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed (0 -> deterministic) for data simulation: ";
            cin>>iSimDataSeed;
        }
        if (iSimDataSeed) rng.reinitialize(iSimDataSeed);
        rpt::echo<<"#Simulating data to fit using "<<iSimDataSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //seed
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-seed"))>-1) {
        if (on+1<argc) {
            iSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed for jittering/resampling: ";
            cin>>iSeed;
        }
        rng.reinitialize(iSeed);
        rpt::echo<<"#Random number seed set to "<<iSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //jitter
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-jitter"))>-1) {
        jitter=1;
        rpt::echo<<"#Jittering for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        rpt::echo<<"#Resampling for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParams"))>-1) {
        debugModelParams=1;
        cout<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        rpt::echo<<"#debugDATA_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPARAMS_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPARAMS_SECTION"))>-1) {
        debugPARAMS_SECTION=1;
        rpt::echo<<"#debugPARAMS_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPRELIM_CALCS
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPRELIM_CALCS"))>-1) {
        debugPRELIM_CALCS=1;
        rpt::echo<<"debugPRELIM_CALCS turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPROC_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPROC_SECTION"))>-1) {
        debugPROC_SECTION=1;
        rpt::echo<<"#debugPROC_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugREPORT_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugREPORT_SECTION"))>-1) {
        debugREPORT_SECTION=1;
        rpt::echo<<"#debugREPORT_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugRunModel
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugRunModel")>-1) {
        debugRunModel=1;
        rpt::echo<<"#debugRunModel turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugObjFun
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugObjFun")>-1) {
        debugObjFun=1;
        rpt::echo<<"#debugObjFun turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //showActiveParams
    if (option_match(ad_comm::argc,ad_comm::argv,"-showActiveParams")>-1) {
        showActiveParams=1;
        rpt::echo<<"#showActiveParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading configuration file '"<<fnConfigFile<<"'"<<endl;
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    if (doRetro){mxYr = mxYr-yRetro; ptrMC->setMaxModelYear(mxYr);}
    if (jitter)   {ptrMC->jitter=1;}
    if (resample) {ptrMC->resample = 1;}
    
    rpt::echo<<"#------------------ModelConfiguration-----------------"<<endl;
    rpt::echo<<(*ptrMC);
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#----finished model configuration---"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelConfiguration-----------------"<<endl;
        cout<<(*ptrMC);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model configuration---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    
    mxYrp1 = mxYr+1;
    nFsh   = ptrMC->nFsh;
    nSrv   = ptrMC->nSrv;
    nZBs   = ptrMC->nZBs;
  zBs.allocate(1,nZBs);
zBs  = ptrMC->zMidPts;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading parameters info file '"<<ptrMC->fnMPI<<"'"<<endl;
    if (debugModelParamsInfo) ModelParametersInfo::debug=1;
    ptrMPI = new ModelParametersInfo(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMPI);
    ptrMPI->read(*(ad_comm::global_datafile));
    if (debugModelParamsInfo) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelParamsInfo;
        if (debugModelParamsInfo<0) exit(1);
        ModelParametersInfo::debug=debugModelParamsInfo;
    }
    rpt::echo<<"#----finished reading model parameters info---"<<endl;
    rpt::echo<<"#-----------ModelParametersInfo---------------"<<endl;
    rpt::echo<<(*ptrMPI)<<endl;
    rpt::echo<<"#----finished ModelParametersInfo---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelParametersInfo-----------------"<<endl;
        cout<<(*ptrMPI);
        cout<<"#----finished model parameters info---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file '"<<ptrMC->fnMDS<<"'"<<endl;
    if (debugModelDatasets) {
        BioData::debug=1;
        FleetData::debug=1;
    }
    ptrMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrMDS->read(*(ad_comm::global_datafile));
    if (debugModelDatasets) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelDatasets;
        if (debugModelDatasets<0) exit(1);
        ModelDatasets::debug=debugModelDatasets;
        BioData::debug=debugModelDatasets;
        FleetData::debug=debugModelDatasets;
    }
    rpt::echo<<"#----finished model datasets---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelDatasets-----------------"<<endl;
        cout<<(*ptrMDS);
        cout<<"#----finished model datasets---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file again to create SimMDS object '"<<ptrMC->fnMDS<<"'"<<endl;
    ptrSimMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrSimMDS->read(*(ad_comm::global_datafile));
    rpt::echo<<"---SimMDS object after reading datasets---"<<endl;
    rpt::echo<<(*ptrSimMDS);
    rpt::echo<<"#finished SimMDS object"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading model options file '"<<ptrMC->fnMOs<<"'"<<endl;
    if (debugModelOptions) ModelOptions::debug=1;
    ptrMOs = new ModelOptions(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMOs);
    ptrMOs->read(*(ad_comm::global_datafile));
    if (debugModelOptions) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelOptions;
        if (debugModelOptions<0) exit(1);
        ModelOptions::debug=debugModelOptions;
    }
    rpt::echo<<"#------------------ModelOptions-----------------"<<endl;
    rpt::echo<<(*ptrMOs);
    rpt::echo<<"#----finished model options---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelOptions-----------------"<<endl;
        cout<<(*ptrMOs);
        cout<<"#----finished model options---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
  mapD2MFsh.allocate(1,nFsh);
  mapM2DFsh.allocate(1,nFsh);
    {
     int idx;
     for (int f=1;f<=nFsh;f++){
         idx = wts::which(ptrMDS->ppFsh[f-1]->name,ptrMC->lblsFsh);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying fishery names and labels in data file and config file."<<endl;
             cout<<"Incorrect fishery name in data file is '"<<ptrMDS->ppFsh[f-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MFsh(f)   = idx;//map from fishery data object f to model fishery idx
         mapM2DFsh(idx) = f;  //map from model fishery idx to fishery data object f
     }
     rpt::echo<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
     cout<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
    }
  mapD2MSrv.allocate(1,nSrv);
  mapM2DSrv.allocate(1,nSrv);
    {
     int idx;
     for (int v=1;v<=nSrv;v++){
         idx = wts::which(ptrMDS->ppSrv[v-1]->name,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in data file and config file."<<endl;
             cout<<"Incorrect survey name in data file is '"<<ptrMDS->ppSrv[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MSrv(v)   = idx;//map from survey data object v to model survey idx
         mapM2DSrv(idx) = v;  //map from model survey idx to survey data object v
     }
     rpt::echo<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
     cout<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
    }
  phsLnR.allocate();
  lbLnR.allocate();
  ubLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnR,npLnR,lbLnR,ubLnR,phsLnR,rpt::echo);
  phsLnRCV.allocate();
  lbLnRCV.allocate();
  ubLnRCV.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRCV,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV,rpt::echo);
  phsLgtRX.allocate();
  lbLgtRX.allocate();
  ubLgtRX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLgtRX,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX,rpt::echo);
  phsLnRa.allocate();
  lbLnRa.allocate();
  ubLnRa.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRa,npLnRa,lbLnRa,ubLnRa,phsLnRa,rpt::echo);
  phsLnRb.allocate();
  lbLnRb.allocate();
  ubLnRb.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRb,npLnRb,lbLnRb,ubLnRb,phsLnRb,rpt::echo);
  mniDevsLnR.allocate();
  mxiDevsLnR.allocate();
  idxsDevsLnR.allocate();
  lbDevsLnR.allocate();
  ubDevsLnR.allocate();
  phsDevsLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pDevsLnR,npDevsLnR,mniDevsLnR,mxiDevsLnR,idxsDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,rpt::echo);
  phsLnM.allocate();
  lbLnM.allocate();
  ubLnM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnM,npLnM,lbLnM,ubLnM,phsLnM,rpt::echo);
  phsLnDMT.allocate();
  lbLnDMT.allocate();
  ubLnDMT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMT,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT,rpt::echo);
  phsLnDMX.allocate();
  lbLnDMX.allocate();
  ubLnDMX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMX,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX,rpt::echo);
  phsLnDMM.allocate();
  lbLnDMM.allocate();
  ubLnDMM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMM,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM,rpt::echo);
  phsLnDMXM.allocate();
  lbLnDMXM.allocate();
  ubLnDMXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMXM,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM,rpt::echo);
zMref = ptrMPI->ptrNM->zRef;
  mniLgtPrMat.allocate();
  mxiLgtPrMat.allocate();
  idxsLgtPrMat.allocate();
  lbLgtPrMat.allocate();
  ubLgtPrMat.allocate();
  phsLgtPrMat.allocate();
tcsam::setParameterInfo(ptrMPI->ptrMat->pLgtPrMat,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,idxsLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,rpt::echo);
  phsLnGrA.allocate();
  lbLnGrA.allocate();
  ubLnGrA.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrA,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA,rpt::echo);
  phsLnGrB.allocate();
  lbLnGrB.allocate();
  ubLnGrB.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrB,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB,rpt::echo);
  phsLnGrBeta.allocate();
  lbLnGrBeta.allocate();
  ubLnGrBeta.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrBeta,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta,rpt::echo);
nSel = ptrMPI->ptrSel->nPCs;//number of selectivity functions defined
  phsS1.allocate();
  lbS1.allocate();
  ubS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS1,npS1,lbS1,ubS1,phsS1,rpt::echo);
  phsS2.allocate();
  lbS2.allocate();
  ubS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS2,npS2,lbS2,ubS2,phsS2,rpt::echo);
  phsS3.allocate();
  lbS3.allocate();
  ubS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS3,npS3,lbS3,ubS3,phsS3,rpt::echo);
  phsS4.allocate();
  lbS4.allocate();
  ubS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS4,npS4,lbS4,ubS4,phsS4,rpt::echo);
  phsS5.allocate();
  lbS5.allocate();
  ubS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS5,npS5,lbS5,ubS5,phsS5,rpt::echo);
  phsS6.allocate();
  lbS6.allocate();
  ubS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS6,npS6,lbS6,ubS6,phsS6,rpt::echo);
  mniDevsS1.allocate();
  mxiDevsS1.allocate();
  idxsDevsS1.allocate();
  lbDevsS1.allocate();
  ubDevsS1.allocate();
  phsDevsS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS1,npDevsS1,mniDevsS1,mxiDevsS1,idxsDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,rpt::echo);
  mniDevsS2.allocate();
  mxiDevsS2.allocate();
  idxsDevsS2.allocate();
  lbDevsS2.allocate();
  ubDevsS2.allocate();
  phsDevsS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS2,npDevsS2,mniDevsS2,mxiDevsS2,idxsDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,rpt::echo);
  mniDevsS3.allocate();
  mxiDevsS3.allocate();
  idxsDevsS3.allocate();
  lbDevsS3.allocate();
  ubDevsS3.allocate();
  phsDevsS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS3,npDevsS3,mniDevsS3,mxiDevsS3,idxsDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,rpt::echo);
  mniDevsS4.allocate();
  mxiDevsS4.allocate();
  idxsDevsS4.allocate();
  lbDevsS4.allocate();
  ubDevsS4.allocate();
  phsDevsS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS4,npDevsS4,mniDevsS4,mxiDevsS4,idxsDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,rpt::echo);
  mniDevsS5.allocate();
  mxiDevsS5.allocate();
  idxsDevsS5.allocate();
  lbDevsS5.allocate();
  ubDevsS5.allocate();
  phsDevsS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS5,npDevsS5,mniDevsS5,mxiDevsS5,idxsDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,rpt::echo);
  mniDevsS6.allocate();
  mxiDevsS6.allocate();
  idxsDevsS6.allocate();
  lbDevsS6.allocate();
  ubDevsS6.allocate();
  phsDevsS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS6,npDevsS6,mniDevsS6,mxiDevsS6,idxsDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,rpt::echo);
  phsHM.allocate();
  lbHM.allocate();
  ubHM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pHM,npHM,lbHM,ubHM,phsHM,rpt::echo);
  phsLnC.allocate();
  lbLnC.allocate();
  ubLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnC,npLnC,lbLnC,ubLnC,phsLnC,rpt::echo);
  phsLnDCT.allocate();
  lbLnDCT.allocate();
  ubLnDCT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCT,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT,rpt::echo);
  phsLnDCX.allocate();
  lbLnDCX.allocate();
  ubLnDCX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCX,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX,rpt::echo);
  phsLnDCM.allocate();
  lbLnDCM.allocate();
  ubLnDCM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCM,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM,rpt::echo);
  phsLnDCXM.allocate();
  lbLnDCXM.allocate();
  ubLnDCXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCXM,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM,rpt::echo);
  mniDevsLnC.allocate();
  mxiDevsLnC.allocate();
  idxsDevsLnC.allocate();
  lbDevsLnC.allocate();
  ubDevsLnC.allocate();
  phsDevsLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDevsLnC,npDevsLnC,mniDevsLnC,mxiDevsLnC,idxsDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,rpt::echo);
  phsLnQ.allocate();
  lbLnQ.allocate();
  ubLnQ.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnQ,npLnQ,lbLnQ,ubLnQ,phsLnQ,rpt::echo);
  phsLnDQT.allocate();
  lbLnDQT.allocate();
  ubLnDQT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQT,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT,rpt::echo);
  phsLnDQX.allocate();
  lbLnDQX.allocate();
  ubLnDQX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQX,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX,rpt::echo);
  phsLnDQM.allocate();
  lbLnDQM.allocate();
  ubLnDQM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQM,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM,rpt::echo);
  phsLnDQXM.allocate();
  lbLnDQXM.allocate();
  ubLnDQXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQXM,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM,rpt::echo);
  dtF_y.allocate(mnYr,mxYr);
dtF_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
  dtM_y.allocate(mnYr,mxYr);
dtM_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
  optsFcAvg.allocate(1,nFsh);
optsFcAvg = ptrMOs->optsFcAvg;
  eff_fy.allocate(1,nFsh,mnYr,mxYr);
  avgEff.allocate(1,nFsh);
optsGrowth = ptrMOs->optsGrowth;
optsInitNatZ = ptrMOs->optsInitNatZ;
npcRec = ptrMPI->ptrRec->nPCs;
npcNM = ptrMPI->ptrNM->nPCs;
npcMat = ptrMPI->ptrMat->nPCs;
npcGr = ptrMPI->ptrGr->nPCs;
npcSel = ptrMPI->ptrSel->nPCs;
npcFsh = ptrMPI->ptrFsh->nPCs;
npcSrv = ptrMPI->ptrSrv->nPCs;
  idxDevsLnC_fy.allocate(1,nFsh,mnYr,mxYr);
yDms = ptrMC->dimYrsToR;//years (mny:mxy)
xDms = ptrMC->dimSXsToR;//sex
mDms = ptrMC->dimMSsToR;//maturity
sDms = ptrMC->dimSCsToR;//shell condition
fDms = ptrMC->dimFshToR;//fisheries
vDms = ptrMC->dimSrvToR;//surveys
ypDms = ptrMC->dimYrsP1ToR;//years (mny:asy)
zbDms = ptrMC->dimZBsToR;//size bin midpoints
zpDms = ptrMC->dimZPsToR;//size bin midpoints (alternative)
zcDms = ptrMC->dimZCsToR;//size bin cuptoints
  hasF_fy.allocate(1,nFsh,mnYr,mxYr);
ptrOFL = new OFLResults();
    rpt::echo<<"#finished DATA_SECTION"<<endl;
    cout<<"#finished DATA_SECTION"<<endl;
}

void model_parameters::initializationfunction(void)
{
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
rpt::echo<<"#Starting PARAMETER_SECTION"<<endl;
cout<<"#Starting PARAMETER_SECTION"<<endl;
  pLnR.allocate(1,npLnR,lbLnR,ubLnR,phsLnR,"pLnR");
  pLnRCV.allocate(1,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV,"pLnRCV");
  pLgtRX.allocate(1,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX,"pLgtRX");
  pLnRa.allocate(1,npLnRa,lbLnRa,ubLnRa,phsLnRa,"pLnRa");
  pLnRb.allocate(1,npLnRb,lbLnRb,ubLnRb,phsLnRb,"pLnRb");
  pDevsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,"pDevsLnR");
  devsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR+1,"devsLnR");
  #ifndef NO_AD_INITIALIZE
    devsLnR.initialize();
  #endif
  pLnM.allocate(1,npLnM,lbLnM,ubLnM,phsLnM,"pLnM");
  pLnDMT.allocate(1,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT,"pLnDMT");
  pLnDMX.allocate(1,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX,"pLnDMX");
  pLnDMM.allocate(1,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM,"pLnDMM");
  pLnDMXM.allocate(1,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM,"pLnDMXM");
  pLnGrA.allocate(1,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA,"pLnGrA");
  pLnGrB.allocate(1,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB,"pLnGrB");
  pLnGrBeta.allocate(1,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta,"pLnGrBeta");
  pLgtPrMat.allocate(1,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,"pLgtPrMat");
  pS1.allocate(1,npS1,lbS1,ubS1,phsS1,"pS1");
  pS2.allocate(1,npS2,lbS2,ubS2,phsS2,"pS2");
  pS3.allocate(1,npS3,lbS3,ubS3,phsS3,"pS3");
  pS4.allocate(1,npS4,lbS4,ubS4,phsS4,"pS4");
  pS5.allocate(1,npS5,lbS5,ubS5,phsS5,"pS5");
  pS6.allocate(1,npS6,lbS6,ubS6,phsS6,"pS6");
  pDevsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,"pDevsS1");
  pDevsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,"pDevsS2");
  pDevsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,"pDevsS3");
  pDevsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,"pDevsS4");
  pDevsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,"pDevsS5");
  pDevsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,"pDevsS6");
  devsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1+1,"devsS1");
  #ifndef NO_AD_INITIALIZE
    devsS1.initialize();
  #endif
  devsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2+1,"devsS2");
  #ifndef NO_AD_INITIALIZE
    devsS2.initialize();
  #endif
  devsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3+1,"devsS3");
  #ifndef NO_AD_INITIALIZE
    devsS3.initialize();
  #endif
  devsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4+1,"devsS4");
  #ifndef NO_AD_INITIALIZE
    devsS4.initialize();
  #endif
  devsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5+1,"devsS5");
  #ifndef NO_AD_INITIALIZE
    devsS5.initialize();
  #endif
  devsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6+1,"devsS6");
  #ifndef NO_AD_INITIALIZE
    devsS6.initialize();
  #endif
  pHM.allocate(1,npHM,lbHM,ubHM,phsHM,"pHM");
  pLnC.allocate(1,npLnC,lbLnC,ubLnC,phsLnC,"pLnC");
  pLnDCT.allocate(1,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT,"pLnDCT");
  pLnDCX.allocate(1,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX,"pLnDCX");
  pLnDCM.allocate(1,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM,"pLnDCM");
  pLnDCXM.allocate(1,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM,"pLnDCXM");
  pDevsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,"pDevsLnC");
  devsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC+1,"devsLnC");
  #ifndef NO_AD_INITIALIZE
    devsLnC.initialize();
  #endif
  pLnQ.allocate(1,npLnQ,lbLnQ,ubLnQ,phsLnQ,"pLnQ");
  pLnDQT.allocate(1,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT,"pLnDQT");
  pLnDQX.allocate(1,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX,"pLnDQX");
  pLnDQM.allocate(1,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM,"pLnDQM");
  pLnDQXM.allocate(1,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM,"pLnDQXM");
  objFun.allocate("objFun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  spB_yx.allocate(mnYr,mxYr,1,nSXs,"spB_yx");
  #ifndef NO_AD_INITIALIZE
    spB_yx.initialize();
  #endif
  n_yxmsz.allocate(mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_yxmsz");
  #ifndef NO_AD_INITIALIZE
    n_yxmsz.initialize();
  #endif
  nmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"nmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    nmN_yxmsz.initialize();
  #endif
  tmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmN_yxmsz.initialize();
  #endif
  initMnR.allocate("initMnR");
  #ifndef NO_AD_INITIALIZE
  initMnR.initialize();
  #endif
  R_y.allocate(mnYr,mxYr,"R_y");
  #ifndef NO_AD_INITIALIZE
    R_y.initialize();
  #endif
  Rx_c.allocate(1,npcRec,"Rx_c");
  #ifndef NO_AD_INITIALIZE
    Rx_c.initialize();
  #endif
  R_yx.allocate(mnYr,mxYr,1,nSXs,"R_yx");
  #ifndef NO_AD_INITIALIZE
    R_yx.initialize();
  #endif
  R_cz.allocate(1,npcRec,1,nZBs,"R_cz");
  #ifndef NO_AD_INITIALIZE
    R_cz.initialize();
  #endif
  R_yz.allocate(mnYr,mxYr,1,nZBs,"R_yz");
  #ifndef NO_AD_INITIALIZE
    R_yz.initialize();
  #endif
  R_yxz.allocate(mnYr,mxYr,1,nSXs,1,nZBs,"R_yxz");
  #ifndef NO_AD_INITIALIZE
    R_yxz.initialize();
  #endif
  stdvDevsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"stdvDevsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    stdvDevsLnR_cy.initialize();
  #endif
  zscrDevsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"zscrDevsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    zscrDevsLnR_cy.initialize();
  #endif
  M_cxm.allocate(1,npcNM,1,nSXs,1,nMSs,"M_cxm");
  #ifndef NO_AD_INITIALIZE
    M_cxm.initialize();
  #endif
  M_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"M_yxmsz");
  #ifndef NO_AD_INITIALIZE
    M_yxmsz.initialize();
  #endif
  prMat_cz.allocate(1,npcMat,1,nZBs,"prMat_cz");
  #ifndef NO_AD_INITIALIZE
    prMat_cz.initialize();
  #endif
  prMat_yxz.allocate(mnYr,mxYr,1,nSXs,1,nZBs,"prMat_yxz");
  #ifndef NO_AD_INITIALIZE
    prMat_yxz.initialize();
  #endif
  mnGrZ_cz.allocate(1,npcGr,1,nZBs,"mnGrZ_cz");
  #ifndef NO_AD_INITIALIZE
    mnGrZ_cz.initialize();
  #endif
  prGr_czz.allocate(1,npcGr,1,nZBs,1,nZBs,"prGr_czz");
  #ifndef NO_AD_INITIALIZE
    prGr_czz.initialize();
  #endif
  mnGrZ_yxsz.allocate(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs,"mnGrZ_yxsz");
  #ifndef NO_AD_INITIALIZE
    mnGrZ_yxsz.initialize();
  #endif
  prGr_yxszz.allocate(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs,1,nZBs,"prGr_yxszz");
  #ifndef NO_AD_INITIALIZE
    prGr_yxszz.initialize();
  #endif
  sel_cz.allocate(1,npcSel,1,nZBs,"sel_cz");
  #ifndef NO_AD_INITIALIZE
    sel_cz.initialize();
  #endif
  sel_iyz.allocate(1,nSel,mnYr,mxYr+1,1,nZBs,"sel_iyz");
  #ifndef NO_AD_INITIALIZE
    sel_iyz.initialize();
  #endif
  dvsLnC_fy.allocate(1,nFsh,mnYr,mxYr,"dvsLnC_fy");
  #ifndef NO_AD_INITIALIZE
    dvsLnC_fy.initialize();
  #endif
  cpF_fxmsy.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr,"cpF_fxmsy");
  #ifndef NO_AD_INITIALIZE
    cpF_fxmsy.initialize();
  #endif
  avgFc_fxms.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,"avgFc_fxms");
  #ifndef NO_AD_INITIALIZE
    avgFc_fxms.initialize();
  #endif
  avgRatioFc2Eff.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,"avgRatioFc2Eff");
  #ifndef NO_AD_INITIALIZE
    avgRatioFc2Eff.initialize();
  #endif
  hmF_fy.allocate(1,nFsh,mnYr,mxYr,"hmF_fy");
  #ifndef NO_AD_INITIALIZE
    hmF_fy.initialize();
  #endif
  cpF_fyxms.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,"cpF_fyxms");
  #ifndef NO_AD_INITIALIZE
    cpF_fyxms.initialize();
  #endif
  sel_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"sel_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    sel_fyxmsz.initialize();
  #endif
  ret_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"ret_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    ret_fyxmsz.initialize();
  #endif
  cpF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cpF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cpF_fyxmsz.initialize();
  #endif
  rmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmF_fyxmsz.initialize();
  #endif
  dmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmF_fyxmsz.initialize();
  #endif
  tmF_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmF_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmF_yxmsz.initialize();
  #endif
  cpN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cpN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cpN_fyxmsz.initialize();
  #endif
  dsN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dsN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dsN_fyxmsz.initialize();
  #endif
  rmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmN_fyxmsz.initialize();
  #endif
  dmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmN_fyxmsz.initialize();
  #endif
  mb_vyx.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,"mb_vyx");
  #ifndef NO_AD_INITIALIZE
    mb_vyx.initialize();
  #endif
  q_vyxms.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,"q_vyxms");
  #ifndef NO_AD_INITIALIZE
    q_vyxms.initialize();
  #endif
  s_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"s_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    s_vyxmsz.initialize();
  #endif
  q_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"q_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    q_vyxmsz.initialize();
  #endif
  n_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    n_vyxmsz.initialize();
  #endif
  fPenRecDevs.allocate(1,npDevsLnR,"fPenRecDevs");
  #ifndef NO_AD_INITIALIZE
    fPenRecDevs.initialize();
  #endif
  fPenSmoothLgtPrMat.allocate(1,npLgtPrMat,"fPenSmoothLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenSmoothLgtPrMat.initialize();
  #endif
  fPenNonDecLgtPrMat.allocate(1,npLgtPrMat,"fPenNonDecLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenNonDecLgtPrMat.initialize();
  #endif
  fPenDevsS1.allocate(1,npDevsS1,"fPenDevsS1");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS1.initialize();
  #endif
  fPenDevsS2.allocate(1,npDevsS2,"fPenDevsS2");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS2.initialize();
  #endif
  fPenDevsS3.allocate(1,npDevsS3,"fPenDevsS3");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS3.initialize();
  #endif
  fPenDevsS4.allocate(1,npDevsS4,"fPenDevsS4");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS4.initialize();
  #endif
  fPenDevsS5.allocate(1,npDevsS5,"fPenDevsS5");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS5.initialize();
  #endif
  fPenDevsS6.allocate(1,npDevsS6,"fPenDevsS6");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS6.initialize();
  #endif
  fPenDevsLnC.allocate(1,npDevsLnC,"fPenDevsLnC");
  #ifndef NO_AD_INITIALIZE
    fPenDevsLnC.initialize();
  #endif
  R_z.allocate(1,nZBs,"R_z");
  #ifndef NO_AD_INITIALIZE
    R_z.initialize();
  #endif
  S1_msz.allocate(1,nMSs,1,nSCs,1,nZBs,"S1_msz");
  #ifndef NO_AD_INITIALIZE
    S1_msz.initialize();
  #endif
  Th_sz.allocate(1,nSCs,1,nZBs,"Th_sz");
  #ifndef NO_AD_INITIALIZE
    Th_sz.initialize();
  #endif
  T_szz.allocate(1,nSCs,1,nZBs,1,nZBs,"T_szz");
  #ifndef NO_AD_INITIALIZE
    T_szz.initialize();
  #endif
  S2_msz.allocate(1,nMSs,1,nSCs,1,nZBs,"S2_msz");
  #ifndef NO_AD_INITIALIZE
    S2_msz.initialize();
  #endif
  n_xmsz.allocate(1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_xmsz");
  #ifndef NO_AD_INITIALIZE
    n_xmsz.initialize();
  #endif
  nllRecDevs.allocate(1,npcRec,"nllRecDevs");
  #ifndef NO_AD_INITIALIZE
    nllRecDevs.initialize();
  #endif
  sdrLnR_y.allocate(mnYr,mxYr,"sdrLnR_y");
  sdrSpB_xy.allocate(1,nSXs,mnYr+5,mxYr,"sdrSpB_xy");
cout<<"#finished PARAMETER_SECTION"<<endl;
rpt::echo<<"#finished PARAMETER_SECTION"<<endl;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
    rpt::echo<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    cout<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    int debug=1;
    
    //set initial values for all parameters
    if (usePin) {
        rpt::echo<<"NOTE: setting initial values for parameters using pin file"<<endl;
    } else {
        rpt::echo<<"NOTE: setting initial values for parameters using setInitVals(...)"<<endl;
        setInitVals();
    }
    cout<<"testing setAllDevs()"<<endl;
    setAllDevs(0,rpt::echo);
        
    {cout<<"writing data to R"<<endl;
     ofstream echo1; echo1.open("ModelData.R", ios::trunc);
     ReportToR_Data(echo1,0,cout);
    }
    
    {cout<<"writing parameters info to R"<<endl;
     ofstream echo1; echo1.open("ModelParametersInfo.R", ios::trunc);
     ptrMPI->writeToR(echo1);
    }
    
    //calculate average effort for fisheries over specified time periods
    avgEff = 0.0;
    for (int fd=1;fd<=nFsh;fd++){//fishery data object
        if (ptrMDS->ppFsh[fd-1]->ptrEff){
            IndexRange* pir = ptrMDS->ppFsh[fd-1]->ptrEff->ptrAvgIR;
            int fm = mapD2MFsh(fd);//index of corresponding model fishery
            int mny = max(mnYr,pir->getMin());
            int mxy = min(mxYr,pir->getMax());
            eff_fy(fm).deallocate();
            eff_fy(fm).allocate(mny,mxy);
            eff_fy(fm) = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(mny,mxy);
            avgEff(fm) = mean(ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(mny,mxy));
        }
    }
    if (debug) {
        rpt::echo<<"eff_fy = "<<endl<<eff_fy<<endl;
        rpt::echo<<"avgEff = "<<avgEff<<endl;
    }
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")<0) {
        cout<<"testing calcRecruitment():"<<endl;
        calcRecruitment(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcNatMort():"<<endl;
        calcNatMort(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcGrowth():"<<endl;
        calcGrowth(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcMaturity():"<<endl;
        calcMaturity(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcSelectivities():"<<endl;
        calcSelectivities(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcFisheryFs():"<<endl;
        calcFisheryFs(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcSurveyQs():"<<endl;
        calcSurveyQs(dbgCalcProcs+1,cout);
        rpt::echo<<"testing runPopDyMod():"<<endl;
        runPopDyMod(dbgCalcProcs+1,cout);
        rpt::echo<<"n_yxm:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                   rpt::echo<<y<<cc;
                   rpt::echo<<tcsam::getSexType(x)<<cc;
                   rpt::echo<<tcsam::getMaturityType(m)<<cc;
                   rpt::echo<<sum(n_yxmsz(y,x,m))<<endl;
                }
            }
        }
        rpt::echo<<"n_yxmsz:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                       rpt::echo<<y<<cc;
                       rpt::echo<<tcsam::getSexType(x)<<cc;
                       rpt::echo<<tcsam::getMaturityType(m)<<cc;
                       rpt::echo<<tcsam::getShellType(s)<<cc;
                       rpt::echo<<n_yxmsz(y,x,m,s)<<endl;
                    }
                }
            }
        }
        
        cout<<"Test OFL calculations"<<endl;
        ofstream echoOFL; echoOFL.open("calcOFL.txt", ios::trunc);
        echoOFL<<"----Testing calcOFL()"<<endl;
        calcOFL(mxYr,100,echoOFL);//updates ptrOFL
        echoOFL<<"----Finished testing calcOFL()!"<<endl;
        echoOFL.close();
        cout<<"Finished testing OFL calculations!"<<endl;
        if (fitSimData){
            cout<<"creating sim data to fit in model"<<endl;
            rpt::echo<<"creating sim data to fit in model"<<endl;
            createSimData(1,rpt::echo,iSimDataSeed,ptrMDS);//stochastic if iSimDataSeed<>0
            {cout<<"re-writing data to R"<<endl;
             rpt::echo<<"re-writing data to R"<<endl;
             ofstream echo1; echo1.open("ModelData.R", ios::trunc);
             ReportToR_Data(echo1,0,cout);
            }
        }
        
        cout<<"Testing calcObjFun()"<<endl;
        rpt::echo<<"Testing calcObjFun()"<<endl;
        calcObjFun(-1,rpt::echo);
        rpt::echo<<"Testing calcObjFun() again"<<endl;
        calcObjFun(dbgAll,rpt::echo);
        
        {cout<<"writing model results to R"<<endl;
            rpt::echo<<"writing model results to R"<<endl;
            ofstream echo1; echo1.open("ModelRes0.R", ios::trunc);
            ReportToR(echo1,1,cout);
        }
        
        {cout<<"writing model sim data to file"<<endl;
            rpt::echo<<"writing model sim data to file"<<endl;
            createSimData(1,rpt::echo,0,ptrSimMDS);//deterministic
            ofstream echo1; echo1.open("ModelSimData0.dat", ios::trunc);
            writeSimData(echo1,0,rpt::echo,ptrSimMDS);
        }
        cout<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
        rpt::echo<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
        rpt::echo<<"#----------------------------------"<<endl<<endl;
    } else {
        writeMCMCHeader();
        cout<<"MCEVAL is on"<<endl;
        rpt::echo<<"MCEVAL is on"<<endl;
    }
    
    
}

void model_parameters::userfunction(void)
{
  objFun =0.0;
    objFun.initialize();
    runPopDyMod(0,rpt::echo);
    calcObjFun(0,rpt::echo);
    
    if (sd_phase()){
        sdrLnR_y = log(R_y);
        for (int x=1;x<=nSXs;x++){
            for (int y=mnYr+ptrMDS->ptrBio->recLag; y<=mxYr; y++){
                sdrSpB_xy(x,y) = spB_yx(y,x);
            }
        }
    }
    
    if (mceval_phase()){
        updateMPI(0, cout);
        calcOFL(mxYr,0,cout);//update ptrOFL
        writeMCMCtoR(mcmc);
    }
}

void model_parameters::setInitVals(void)
{
    //recruitment parameters
    setInitVals(ptrMPI->ptrRec->pLnR,    pLnR,    0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRCV,  pLnRCV,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLgtRX,  pLgtRX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRa,   pLnRa,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRb,   pLnRb,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,0,rpt::echo);
    //natural mortality parameters
    setInitVals(ptrMPI->ptrNM->pLnM,   pLnM,   0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMT, pLnDMT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMX, pLnDMX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMM, pLnDMM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMXM,pLnDMXM,0,rpt::echo);
    //growth parameters
    setInitVals(ptrMPI->ptrGr->pLnGrA,   pLnGrA,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGr->pLnGrB,   pLnGrB,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGr->pLnGrBeta,pLnGrBeta,0,rpt::echo);
    //maturity parameters
    setInitVals(ptrMPI->ptrMat->pLgtPrMat,pLgtPrMat,0,rpt::echo);
    //selectivity parameters
    setInitVals(ptrMPI->ptrSel->pS1, pS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS2, pS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS3, pS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS4, pS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS5, pS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS6, pS6,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS1, pDevsS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS2, pDevsS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS3, pDevsS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS4, pDevsS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS5, pDevsS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS6, pDevsS6,0,rpt::echo);
    //fully-selected fishing capture rate parameters
    setInitVals(ptrMPI->ptrFsh->pHM,     pHM,     0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnC,    pLnC,    0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCT,  pLnDCT,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCX,  pLnDCX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCM,  pLnDCM,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCXM, pLnDCXM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,0,rpt::echo);
    //survey catchability parameters
    setInitVals(ptrMPI->ptrSrv->pLnQ,   pLnQ,   0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQT, pLnDQT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQX, pLnDQX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQM, pLnDQM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQXM,pLnDQXM,0,rpt::echo);
}

void model_parameters::writeMCMCHeader(void)
{
    mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::trunc);
    mcmc<<"mcmc=list("<<endl;
    mcmc.close();
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,NumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedNumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,DevsVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ofstream& mcmc)
{
    mcmc.open((char *) fnMCMC, ofstream::out|ofstream::app);
    mcmc<<"list(objFun="<<objFun<<cc<<endl;
    //write parameter values
        //recruitment values
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnR);   mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRCV); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLgtRX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRa);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRb);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pDevsLnR);  mcmc<<cc<<endl;
        //natural mortality parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMXM); mcmc<<cc<<endl;
        //growth parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrA); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrB); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrBeta); mcmc<<cc<<endl;
        //maturity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrMat->pLgtPrMat); mcmc<<cc<<endl;
        //selectivity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS6); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS6); mcmc<<cc<<endl;
        //fully-selected fishing capture rate parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pHM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCXM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDevsLnC); mcmc<<cc<<endl;
        //survey catchability parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnQ);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQT);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQX);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQM);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQXM); mcmc<<cc<<endl;
    
        //write other quantities
        mcmc<<"R_y="; wts::writeToR(mcmc,value(R_y)); mcmc<<cc<<endl;
        ivector bnds = wts::getBounds(spB_yx);
        mcmc<<"MB_xy="; wts::writeToR(mcmc,trans(value(spB_yx)),xDms,yDms); mcmc<<cc<<endl;
        ptrOFL->writeToR(mcmc,"oflResults",0);//mcm<<cc<<endl;
        
    mcmc<<")"<<cc<<endl;
    mcmc.close();
    
}

void model_parameters::createSimData(int debug, ostream& cout, int iSimDataSeed, ModelDatasets* ptrSim)
{
    if (debug)cout<<"simulating model results as data"<<endl;
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vcN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    for (int f=1;f<=nFsh;f++) {
        if (debug) cout<<"fishery f: "<<f<<endl;
        (ptrSim->ppFsh[f-1])->replaceFisheryCatchData(iSimDataSeed,rngSimData,vcN_fyxmsz(f),vrmN_fyxmsz(f),ptrSim->ptrBio->wAtZ_xmz);
    }
    for (int v=1;v<=nSrv;v++) {
        if (debug) cout<<"survey "<<v<<endl;
        (ptrSim->ppSrv[v-1])->replaceIndexCatchData(iSimDataSeed,rngSimData,vn_vyxmsz(v),ptrSim->ptrBio->wAtZ_xmz);
    }
    if (debug) cout<<"finished simulating model results as data"<<endl;
     
}

void model_parameters::writeSimData(ostream& os, int debug, ostream& cout, ModelDatasets* ptrSim)
{
    if (debug)cout<<"writing model results as data"<<endl;
    for (int v=1;v<=nSrv;v++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppSrv[v-1]))<<endl;
    }
    //     cout<<4<<endl;
    for (int f=1;f<=nFsh;f++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppFsh[f-1]))<<endl;
    }
    if (debug) cout<<"finished writing model results as data"<<endl;
     
}

void model_parameters::setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = pI->getInitVals();//defaults are initial values
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            NumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);//assign initial value based on resampling prior pdf
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
     
}

void model_parameters::setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = 0.5*(pI->getUpperBounds()+pI->getLowerBounds());//defaults are midpoints of ranges
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            BoundedNumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<"jittering "<<p(i).label()<<endl;
                p(i) = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//only done if parameter phase > 0
            } else 
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            BoundedVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            DevsVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setAllDevs(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"starting setAllDevs()"<<endl;
    tcsam::setDevs(devsLnR, pDevsLnR,debug,cout);
    tcsam::setDevs(devsS1, pDevsS1,debug,cout);
    tcsam::setDevs(devsS2, pDevsS2,debug,cout);
    tcsam::setDevs(devsS3, pDevsS3,debug,cout);
    tcsam::setDevs(devsS4, pDevsS4,debug,cout);
    tcsam::setDevs(devsS5, pDevsS5,debug,cout);
    tcsam::setDevs(devsS6, pDevsS6,debug,cout);
    
    tcsam::setDevs(devsLnC, pDevsLnC,debug,cout);
    if (debug>=dbgAll) cout<<"finished setAllDevs()"<<endl;
}

void model_parameters::calcEqNatZF100(dvariable& R, int yr, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting dvar calcEqNatZF100()"<<endl;
    n_xmsz.initialize();//equilibrium n-at-z
    for (int x=1;x<=nSXs;x++){
        S1_msz.initialize(); //survival until molting/mating
        Th_sz.initialize();  //pr(molt to maturity|pre-molt size, molt)
        T_szz.initialize();  //growth matrices (indep. of molt to maturity)
        S2_msz.initialize(); //survival after molting/mating
        R_z.initialize();    //recruitment size distribution
        R_z = R*R_yx(yr,x)*R_yz(yr);//initial mean recruitment by size
        for (int s=1;s<=nSCs;s++){
            Th_sz(s) = prMat_yxz(yr,x); //pr(molt to maturity|pre-molt size, molt)
            for (int z=1;z<=nZBs;z++) T_szz(s,z) = prGr_yxszz(yr,x,s,z);//growth matrices
            for (int m=1;m<=nMSs;m++){ 
                S1_msz(m,s) = exp(-M_yxmsz(yr,x,m,s)*dtM_y(yr));      //survival until molting/growth/mating
                S2_msz(m,s) = exp(-M_yxmsz(yr,x,m,s)*(1.0-dtM_y(yr)));//survival after molting/growth/mating
            }//m
        }//s
        n_xmsz(x) = calcEqNatZ(R_z, S1_msz, Th_sz, T_szz, S2_msz, debug, cout);
    }
    
    if (debug>=dbgPopDy) cout<<"finished dvar calcEqNatZF100()"<<endl;
}

dvar3_array model_parameters::calcEqNatZ(dvar_vector& R_z,dvar3_array& S1_msz, dvar_matrix& Th_sz, dvar3_array& T_szz, dvar3_array& S2_msz, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgPopDy) cout<<"starting calcEqNatZ()"<<endl;
    //the equilibrium solution
    dvar3_array n_msz(1,nMSs,1,nSCs,1,nZBs); n_msz.initialize();
    
    //create an identity matrix
    dmatrix I = identity_matrix(1,nZBs);
    
    //--calc the state transition matrices
    int i = IMMATURE; 
    int m =   MATURE;
    int n = NEW_SHELL;
    int o = OLD_SHELL;
    //immature new shell crab
    dvar_matrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dvar_matrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dvar_matrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|pre-molt size,new shell, molting)
    dvar_matrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dvar_matrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dvar_matrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dvar_matrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|pre-molt size,old shell, molting)
    dvar_matrix Ph_io = identity_matrix(1,nZBs);//pr(molt|pre-molt size, old shell) [assumed all immature crab molt]
    dvar_matrix S1_io = wts::diag(S1_msz(i,o)); //pr(survival|size) for immature old shell crab before molting occurs
    //mature new shell crab
    dvar_matrix Tr_mn = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix Tr_mo = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix S2_mn = wts::diag(S2_msz(m,n)); //pr(survival|size) for mature new shell crab after molting occurs
    dvar_matrix S1_mn = wts::diag(S1_msz(m,n)); //pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    //mature old shell crab
    dvar_matrix S2_mo = wts::diag(S2_msz(m,o)); //pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix S1_mo = wts::diag(S1_msz(m,o)); //pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    
    //full state transition matrices
    dvar_matrix lA = S2_in * Tr_in * (I-Th_in) * Ph_in * S1_in;//imm, new -> imm, new
    dvar_matrix lB = S2_in * Tr_io * (I-Th_io) * Ph_io * S1_io;//imm, old -> imm, new
    dvar_matrix lC = S2_io * (I-Ph_in) * S1_in;                //imm, new -> imm, old
    dvar_matrix lD = S2_io * (I-Ph_io) * S1_io;                //imm, old -> imm, old
    dvar_matrix lE = S2_mn * Tr_mn * Th_in * Ph_in * S1_in;    //imm, new -> mat, new (terminal molt)
    dvar_matrix lF = S2_mn * Tr_mo * Th_io * Ph_io * S1_io;    //imm, old -> mat, new (terminal molt)
    dvar_matrix lG = S2_mo * S1_mn;                            //mat, new -> mat, old
    dvar_matrix lH = S2_mo * S1_mo;                            //mat, old -> mat, old
    //--done calculating transition matrices
    
    //calculate inverses of matrix quantities
    dvar_matrix iM1 = inv(I - lD);
    dvar_matrix iM2 = inv(I - lA - lB * iM1 * lC);
    dvar_matrix iM3 = inv(I - lH);
    
    //the equilibrium solution is
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
        
    if (debug>=dbgPopDy) cout<<"finished calcEqNatZ()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return(n_msz);
}

void model_parameters::calcOFL(int yr, int debug, ostream& cout)
{
    if (debug>=dbgOFL) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL(yr,debug,cout)"<<endl;
    }
    //get initial population for "upcoming" year, yr
    d4_array n_xmsz = wts::value(n_yxmsz(yr));
    
    //NOW set yr back one year to get population rates, etc.
    yr = yr - 1;
    
    //1. Determine population rates for next year, using yr
    double dtF = dtF_y(yr);
    double dtM = dtM_y(yr);
    
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z   = value(R_yz(yr));
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = value(M_yxmsz(yr,MALE));
    pPIM->T_szz = value(prGr_yxszz(yr,MALE));
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = value(prMat_yxz(yr,MALE));
    
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z   = value(R_yz(yr));
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = value(M_yxmsz(yr,FEMALE));
    pPIF->T_szz = value(prGr_yxszz(yr,FEMALE));
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = value(prMat_yxz(yr,FEMALE));
    
    //2. Determine fishery conditions for next year based on averages for recent years
        int oflAvgPeriodYrs = 5;
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvector avgHM_f(1,nFsh);
        avgHM_f.initialize();
        for (int f=1;f<=nFsh;f++){
            ny = 0;
            for (int y=yr-oflAvgPeriodYrs+1;y<=yr;y++){
                ny         += hasF_fy(f,y);
                avgHM_f(f) += value(hmF_fy(f,y));
            }
            avgHM_f(f) /= 1.0*ny;
        }
        if (debug>=dbgOFL) cout<<"avgHm_f = "<<avgHM_f<<endl;
        d5_array avgRFcn_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//averaged retention function
        d5_array avgCapF_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//averaged capture mortality
        avgRFcn_fxmsz.initialize();
        avgCapF_fxmsz.initialize();
        for (int f=1;f<=nFsh;f++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) {
                        for (int z=1;z<=nZBs;z++){
                            ny = 0;
                            for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                                ny += hasF_fy(f,y);
                                avgRFcn_fxmsz(f,x,m,s,z) += value(ret_fyxmsz(f,y,x,m,s,z));
                                avgCapF_fxmsz(f,x,m,s,z) += value(cpF_fyxmsz(f,y,x,m,s,z));
                            }
                            avgCapF_fxmsz(f,x,m,s,z) /= 1.0*ny;
                            avgRFcn_fxmsz(f,x,m,s,z) /= 1.0*ny;
                        }
                    }
                }
            }
        }
        cout<<"avgCapF_fxmsz(1,MALE,MATURE,NEW_SHELL) = "<<avgCapF_fxmsz(1,MALE,MATURE,NEW_SHELL)<<endl;
        cout<<"avgRFcn_fxmsz(1,MALE,MATURE,NEW_SHELL) = "<<avgRFcn_fxmsz(1,MALE,MATURE,NEW_SHELL)<<endl;
        
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
        pCIM->setCaptureRates(MALE, avgCapF_fxmsz);
        pCIM->setRetentionFcns(MALE, avgRFcn_fxmsz);
        pCIM->setHandlingMortality(avgHM_f);
        double maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug>=dbgOFL) cout<<"maxCapF = "<<maxCapF<<endl;
        
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
        pCIF->setCaptureRates(FEMALE, avgCapF_fxmsz);
        pCIF->setRetentionFcns(FEMALE, avgRFcn_fxmsz);
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
        
    //3. Determine TIER LEVEL
        int tier = 3;
        
    //4. Determine mean recruitment
        dvector avgRec_x(1,nSXs);
        for (int x=1;x<=nSXs;x++) 
            avgRec_x(x)= value(mean(elem_prod(R_y(1982,mxYr),column(R_yx,x)(1982,mxYr))));
        if (debug>=dbgOFL) {
            cout<<"exp(mnLnR) = "<<exp(pLnR(1))<<endl;
            cout<<"R_y(1982,mxYr) = "<<R_y(1982,mxYr)<<endl;
            cout<<"R_yx((1982:mxYr,MALE) = "<<column(R_yx,MALE)(1982,mxYr)<<endl;
            cout<<"Average recruitment = "<<avgRec_x<<endl;
        }
        
    //5. Determine Fmsy and Bmsy
        Tier3_Calculator* pT3C = new Tier3_Calculator(nZBs,nFsh);
        pT3C->dtF = dtF;
        pT3C->dtM = dtM;
        pT3C->pPI = pPIM;//male pop dy info
        pT3C->pCI = pCIM;//male catch info
        
        double B100 = pT3C->calcB100(avgRec_x(MALE),cout);
        double Bmsy = pT3C->calcBmsy(avgRec_x(MALE),cout);
        double Fmsy = pT3C->calcFmsy(avgRec_x(MALE),cout);
        
        if (debug>=dbgOFL){
            d3_array tmp0_msz = pT3C->calcEqNatZF0(avgRec_x(MALE),cout);
            double tmpMMB0 = pT3C->pPI->calcMatureBiomass(tmp0_msz,cout);
            d3_array tmp1_msz = pT3C->calcEqNatZFM(avgRec_x(MALE),Fmsy,cout);
            double tmpMMB = pT3C->pPI->calcMatureBiomass(tmp1_msz,cout);
            cout<<"B100 = "<<B100<<endl;
            cout<<"Bmsy = "<<Bmsy<<endl;
            cout<<"Fmsy = "<<Fmsy<<endl<<endl<<endl;
            cout<<"Unfished equilibrium size distribution:"<<endl;
            wts::print(tmp0_msz,cout,2);
            cout<<"MMB100 = "<<tmpMMB0<<endl;
            cout<<"Fmsy equilibrium size distribution:"<<endl;
            wts::print(tmp1_msz,cout,2);
            cout<<"MMBmsy = "<<tmpMMB<<endl;
            cout<<"------------------------------------------"<<endl<<endl;
        }
        
    //6. Determine Fofl and OFL    
        //population projector for males
        PopProjector* pPPM = new PopProjector(nZBs,nFsh);
        pPPM->dtF = dtF;
        pPPM->dtM = dtM;
        pPPM->pPI = pPIM;//male pop dy info
        pPPM->pCI = pCIM;//male catch info
        
        //population projector for females
        PopProjector* pPPF = new PopProjector(nZBs,nFsh);
        pPPF->dtF = dtF;
        pPPF->dtM = dtM;
        pPPF->pPI = pPIF;//female pop dy info
        pPPF->pCI = pCIF;//female catch info
        
        //OFL calculator
        OFL_Calculator* pOC = new OFL_Calculator(nFsh);
        pOC->pT3C = pT3C;
        pOC->pPrjM = pPPM;
        pOC->pPrjF = pPPF;
        
        //calculate Fofl
        double Fofl = pOC->calcFofl(Bmsy,Fmsy,n_xmsz(MALE),cout);
        if (debug>=dbgOFL) cout<<"Fofl = "<<Fofl<<endl;
        //calculate OFL
        double OFL = pOC->calcOFL(Fofl,n_xmsz,cout);
        if (debug>=dbgOFL) {
            cout<<"OFL = "<<OFL<<endl;
            cout<<"retained catch:"<<tb<<pOC->ofl_fx(0,MALE)<<endl;
            for (int f=1;f<=nFsh;f++){
                cout<<"fishery "<<f<<":"<<tb<<pOC->ofl_fx(f,MALE)<<tb<<pOC->ofl_fx(f,FEMALE)<<endl;
            }
        }
        //calculate projected ("current") MMB
        double prjMMB = pOC->calcPrjMMB(Fofl,n_xmsz(MALE),cout);
        if (debug>=dbgOFL) cout<<"prjMMB = "<<prjMMB<<endl;
        
    //encapsulate results
    ptrOFL->B100 = B100;
    ptrOFL->Bmsy = Bmsy;
    ptrOFL->Fmsy = Fmsy;
    ptrOFL->Fofl = Fofl;
    ptrOFL->OFL  = OFL;
    ptrOFL->prjB = prjMMB;
    
    if (debug>=dbgOFL) {
        cout<<"finished calcOFL(yr,debug,cout)"<<endl;
        cout<<"#------------------------"<<endl<<endl<<endl;
    }
        
}

void model_parameters::initPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting initPopDyMod()"<<endl;
    
    spB_yx.initialize();
    n_yxmsz.initialize();
    nmN_yxmsz.initialize();
    tmN_yxmsz.initialize();
       
    setAllDevs(debug,cout);//set devs vectors
    
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcMaturity(debug,cout);   //calculate maturity ogives
    
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    
    if (optsInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (optsInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<optsInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }
    
    if (debug>=dbgPopDy) cout<<"finished initPopDyMod()"<<endl;
}

void model_parameters::runPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting runPopDyMod()"<<endl;
    //initialize population model
    initPopDyMod(debug, cout);
    
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,debug,cout);
        runPopDyModOneYear(y,debug,cout);        
    }
    doSurveys(mxYr+1,debug,cout);//do final surveys
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyMod()"<<endl;
    
}

void model_parameters::doSurveys(int y,int debug,ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting doSurveys("<<y<<")"<<endl;
    for (int v=1;v<=nSrv;v++){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    n_vyxmsz(v,y,x,m,s) = elem_prod(q_vyxmsz(v,y,x,m,s),n_yxmsz(y,x,m,s));
                }
            }
        }
    }
    for (int v=1;v<=nSrv;v++){
        mb_vyx(v,y) = calcSpB(n_vyxmsz(v,y),y,debug,cout);
    }
    if (debug>=dbgPopDy) cout<<"finished doSurveys("<<y<<")"<<endl;
}

void model_parameters::runPopDyModOneYear(int yr, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"Starting runPopDyModOneYear("<<yr<<")"<<endl;
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n2_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n3_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n4_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n5_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    
    if (dtF_y(yr)<=dtM_y(yr)){//fishery occurs BEFORE molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs BEFORE molting/growth/maturity"<<endl;
        //apply natural mortality before fisheries
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtF_y(yr),debug,cout);
        //conduct fisheries
        n2_xmsz = applyFshMort(n1_xmsz,yr,debug,cout);
        //apply natural mortality from fisheries to molting/growth/maturity
        if (dtF_y(yr)==dtM_y(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtM_y(yr)-dtF_y(yr),debug,cout);
        }
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spB_yx(yr) = calcSpB(n3_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n4_xmsz = applyMGM(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtM_y(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtM_y(yr),debug,cout);
        }
    } else {              //fishery occurs AFTER molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs AFTER molting/growth/maturity"<<endl;
        //apply natural mortality before molting/growth/maturity
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtM_y(yr),debug,cout);
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spB_yx(yr) = calcSpB(n1_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n2_xmsz = applyMGM(n1_xmsz,yr,debug,cout);
        //apply natural mortality from molting/growth/maturity to fisheries
        if (dtM_y(yr)==dtF_y(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtF_y(yr)-dtM_y(yr),debug,cout);
        }
        //conduct fisheries
        n4_xmsz = applyFshMort(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtF_y(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtF_y(yr),debug,cout);
        }
    }
    
    //advance surviving individuals to next year
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n_yxmsz(yr+1,x,m,s) = n5_xmsz(x,m,s);
            }
        }
    }
    //add in recruits
    for (int x=1;x<=nSXs;x++) n_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) += R_yxz(yr,x);
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyModOneYear("<<yr<<")"<<endl;
    
}

dvar_vector model_parameters::calcSpB(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector spb(1,nSXs); spb.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int s=1;s<=nSCs;s++) spb(x) += n0_xmsz(x,MATURE,s)*ptrMDS->ptrBio->wAtZ_xmz(x,MATURE);//dot product here
    }
    if (debug>dbgApply) cout<<"finished calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return spb;
    
}

dvar4_array model_parameters::applyNatMort(dvar4_array& n0_xmsz, int y, double dt, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n1_xmsz(x,m,s) = elem_prod(exp(-M_yxmsz(y,x,m,s)*dt),n0_xmsz(x,m,s));//survivors
                nmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
                tmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

dvar4_array model_parameters::applyFshMort(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector tm_z(1,nZBs);//total mortality (numbers) by size
    dvar_vector tvF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvector     tdF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                tmF_yxmsz(y,x,m,s) = 0.0;//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) tmF_yxmsz(y,x,m,s) += rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(exp(-tmF_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheries
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);  //numbers killed by all fisheries
                tmN_yxmsz(y,x,m,s) += tm_z;            //add in numbers killed by all fisheries to total killed
                
                //calculate fishing rate components (need to ensure NOT dividing by 0)
                tdF_z = value(tmF_yxmsz(y,x,m,s));
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),tmF_yxmsz(y,x,m,s)) + 
                                  wts::isEQ(tdF_z,0.0);
                for (int f=1;f<=nFsh;f++){                   
                    cpN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div( cpF_fyxmsz(f,y,x,m,s),tvF_z),tm_z);//numbers captured in fishery f
                    rmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(rmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z); //retained mortality in fishery f (numbers)
                    dmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(dmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z); //discards mortality in fishery f (numbers)
                    dsN_fyxmsz(f,y,x,m,s) = cpN_fyxmsz(f,y,x,m,s)-rmN_fyxmsz(f,y,x,m,s);//discarded catch (NOT mortality) in fishery f (numbers)                    
                }
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

dvar4_array model_parameters::applyMGM(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        n1_xmsz(x,IMMATURE,NEW_SHELL) = prGr_yxszz(y,x,NEW_SHELL)*elem_prod(1.0-prMat_yxz(y,x),n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,IMMATURE,OLD_SHELL) = 0.0;
        n1_xmsz(x,MATURE,NEW_SHELL)   = prGr_yxszz(y,x,NEW_SHELL)*elem_prod(    prMat_yxz(y,x),n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,MATURE,OLD_SHELL)   = n0_xmsz(x,MATURE,NEW_SHELL)+n0_xmsz(x,MATURE,OLD_SHELL);
    }
    if (debug>dbgApply) cout<<"finished applyNatMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

void model_parameters::calcRecruitment(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcRecruitment()"<<endl;
    RecruitmentInfo* ptrRI = ptrMPI->ptrRec;
    
    R_y.initialize();
    Rx_c.initialize();
    R_yx.initialize();
    R_cz.initialize();
    R_yz.initialize();
    R_yxz.initialize();
    stdvDevsLnR_cy.initialize();
    zscrDevsLnR_cy.initialize();
    
    int k; int y;
    dvector dzs = zBs+(zBs[2]-zBs[1])/2.0-zBs[1];
    for (int pc=1;pc<=ptrRI->nPCs;pc++){
        ivector pids = ptrRI->getPCIDs(pc);
        k=ptrRI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvariable mnLnR    = pLnR(pids[k++]);
        dvariable lnRCV    = pLnRCV(pids[k++]);
        dvariable lgtRX    = pLgtRX(pids[k++]);
        dvariable lnRa     = pLnRa(pids[k++]);
        dvariable lnRb     = pLnRb(pids[k++]);
        if (debug>dbgCalcProcs){
            cout<<"pids  = "<<pids<<endl;
            cout<<"mnLnR = "<<mnLnR<<endl;
            cout<<"lnRCV = "<<lnRCV<<endl;
            cout<<"lgtRX = "<<lgtRX<<endl;
            cout<<"lnRa  = "<<lnRa<<endl;
            cout<<"lnRb  = "<<lnRb<<endl;
        }
        int useDevs = pids[k]; k++;
        dvariable mnR;   //mean recruitment
        dvariable varLnR;//ln-scale variance in recruitment
        dvar_vector dvsLnR;
        ivector idxDevsLnR;
        varLnR = log(1.0+exp(2.0*lnRCV));//ln-scale variance
        mnR    = exp(mnLnR+varLnR/2.0);  //mean recruitment
        if (useDevs) {
            dvsLnR     = devsLnR(useDevs);
            idxDevsLnR = idxsDevsLnR(useDevs);
            if (debug>dbgCalcProcs) {
                cout<<"lims(dvsLnR) = "<<dvsLnR.indexmin()<<cc<<dvsLnR.indexmax()<<endl;
                cout<<"idx(dvsLnR) = "<<idxDevsLnR<<endl;
                cout<<"dvsLnR = "<<dvsLnR<<endl;
            }
        }
        
        Rx_c(pc) = 1.0/(1.0+exp(-lgtRX));
        R_cz(pc) = elem_prod(pow(dzs,exp(lnRa-lnRb)-1.0),exp(-dzs/exp(lnRb)));
        R_cz(pc) /= sum(R_cz(pc));//normalize to sum to 1
        imatrix idxs = ptrRI->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if (y==mnYr) initMnR = mnR;
            if ((mnYr<=y)&&(y<=mxYr)){
                if (debug>dbgCalcProcs+10) cout<<"y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
                if (useDevs){
                    R_y(y) = exp(mnLnR+dvsLnR[idxDevsLnR[y]]);
                } else {
                    R_y(y) = mnR;
                }
                if (debug>dbgCalcProcs+10) cout<<"R_y(y)="<<R_y(y)<<tb;
                if (MALE==nSXs){
                    R_yx(y,MALE) = 1.0;//only tracking males
                } else {
                    R_yx(y,MALE)   = Rx_c(pc);
                    R_yx(y,FEMALE) = 1.0-R_yx(y,MALE);
                    if (debug>dbgCalcProcs+10) cout<<R_yx(y,MALE)<<endl;
                }
                R_yz(y) = R_cz(pc);
                if (debug>dbgCalcProcs+10) cout<<"R_yz(y)="<<R_yz(y)<<endl;
                
                for (int x=1;x<=nSXs;x++) R_yxz(y,x) = R_y(y)*R_yx(y,x)*R_yz(y);
                stdvDevsLnR_cy(pc,y) = sqrt(varLnR); //ln-scale std dev
                zscrDevsLnR_cy(pc,y) = dvsLnR(idxDevsLnR(y))/stdvDevsLnR_cy(pc,y);//standardized ln-scale rec devs
            } else {
                if (debug>dbgCalcProcs) cout<<"skipping y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            }
        }//idx
    }//pc
    
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz = "<<endl<<R_yz<<endl;
        cout<<"zscr = "<<zscrDevsLnR_cy<<endl;
        cout<<"finished calcRecruitment()"<<endl;
    }
}

void model_parameters::calcNatMort(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcNatMort()"<<endl;
    
    NaturalMortalityInfo* ptrNM = ptrMPI->ptrNM;
    
    dvar_matrix lnM(1,nSXs,1,nMSs);
    dvar3_array M_xmz(1,nSXs,1,nMSs,1,nZBs);
    
    M_cxm.initialize();
    M_yxmsz.initialize();
    int y; 
    for (int pc=1;pc<=ptrNM->nPCs;pc++){
        lnM.initialize();
        ivector pids = ptrNM->getPCIDs(pc);
        int k=ptrNM->nIVs+1;//1st parameter variable column
        //add in base (ln-scale) natural mortality (mature males)
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x) += pLnM(pids[k]);}   k++;
        //add in main temporal offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x) += pLnDMT(pids[k]);} k++;
        if (FEMALE<=nSXs){
            //add in female offset
            if (pids[k]) {lnM(FEMALE) += pLnDMX(pids[k]);}                      k++;
            //add in immature offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x,IMMATURE) += pLnDMM(pids[k]);} k++;
            //add in offset immature females for stanza
            if (pids[k]) {lnM(FEMALE,IMMATURE) += pLnDMXM(pids[k]);}            k++; //advance k to zScaling in pids
        }
        
        //convert from ln-scale to arithmetic scale
        M_cxm(pc) = exp(lnM);
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"lnM:"<<endl<<lnM<<endl;
            cout<<"pc: "<<pc<<tb<<"M_xm:"<<endl<<M_cxm(pc)<<endl;
        }
        
        //add in size-scaling, if requested
        M_xmz.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                if (pids[k]&&(current_phase()>=pids[k])) {
                    M_xmz(x,m) = M_cxm(pc,x,m)*(zMref/zBs);//factor in size dependence
                } else {
                    M_xmz(x,m) = M_cxm(pc,x,m);//no size dependence
                }
            }
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrNM->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//only model index for natural mortality is year
            if ((mnYr<=y)&&(y<=mxYr)){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) M_yxmsz(y,x,m,s) = M_xmz(x,m);
                    }
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"Finished calcNatMort()"<<endl;
    
}

void model_parameters::calcMaturity(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcMaturity()"<<endl;
    MaturityInfo* ptrMI = ptrMPI->ptrMat;
    
    prMat_cz.initialize();
    prMat_yxz.initialize();
    
    int k; int y; int x;
    for (int pc=1;pc<=ptrMI->nPCs;pc++){
        ivector pids = ptrMI->getPCIDs(pc);
        k=ptrMI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvar_vector lgtPrMat = pLgtPrMat(pids[k++]);
        int vmn = lgtPrMat.indexmin();
        int vmx = lgtPrMat.indexmax();
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<vmn<<", mx = "<<vmx<<endl;
            cout<<"lgtPrMat = "<<lgtPrMat<<endl;
        }
        prMat_cz(pc) = 1.0;//default is 1
        prMat_cz(pc)(vmn,vmx) = 1.0/(1.0+exp(-lgtPrMat));
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<vmn<<", mx = "<<vmx<<endl;
            cout<<"prMat = "<<prMat_cz(pc)<<endl;
        }
        
        imatrix idxs = ptrMI->getModelIndices(pc);
        if (debug>dbgCalcProcs) cout<<"maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2);
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(x)<<endl;
                prMat_yxz(y,x) = prMat_cz(pc);//note: this change made a difference, but not sure why!
            }
        }
    }
    
    if (debug>dbgCalcProcs) cout<<"finished calcMaturity()"<<endl;
}

void model_parameters::calcGrowth(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcGrowth()"<<endl;
    
    GrowthInfo* ptrGrI = ptrMPI->ptrGr;
    
    dvariable grA;
    dvariable grB;
    dvariable grBeta;
    
    mnGrZ_cz.initialize();
    prGr_czz.initialize();
    mnGrZ_yxsz.initialize();
    prGr_yxszz.initialize();
    
    dvar_matrix prGr_zz(1,nZBs,1,nZBs);
    int y; int x;
    for (int pc=1;pc<=ptrGrI->nPCs;pc++){
        ivector pids = ptrGrI->getPCIDs(pc);
        int k=ptrGrI->nIVs+1;//1st parameter column
        grA = exp(pLnGrA(pids[k])); k++; //"a" coefficient for mean growth
        grB = exp(pLnGrB(pids[k])); k++; //"b" coefficient for mean growth
        grBeta = exp(pLnGrBeta(pids[k])); k++; //shape factor for gamma function growth transition
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<". grB:"<<tb<<grB<<". grBeta:"<<grBeta<<endl;
        }
        
        //compute growth transition matrix for this pc
        prGr_zz.initialize();
        dvar_vector mnZ = exp(grA)*pow(zBs,grB);//mean size after growth from zBs
        mnGrZ_cz(pc) = mnZ;
        if (optsGrowth==0) {
            //old style (TCSAM2013)
            dvar_vector alZ = (mnZ-zBs)/grBeta;//scaled mean growth increment from zBs
            for (int z=1;z<nZBs;z++){//pre-molt growth bin
                dvar_vector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
                if (debug) cout<<"dZs: "<<dZs.indexmin()<<":"<<dZs.indexmax()<<endl;
                dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),exp(-dZs/grBeta)); //pr(dZ|z)
                if (debug) cout<<"prs: "<<prs.indexmin()<<":"<<prs.indexmax()<<endl;
                if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                if (debug) cout<<prs<<endl;
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else if (optsGrowth==1){
            //use cumd_gamma function like gmacs
            dvar_vector sclMnZ = mnZ/grBeta;           //scaled mean growth increments
            dvar_vector sclZCs = ptrMC->zCutPts/grBeta;//scaled size bin cut points
            for (int z=1;z<nZBs;z++){
                dvar_vector cprs(z,nZBs);
                for (int zp=z;zp<=nZBs;zp++){
                    cprs(zp) = cumd_gamma(sclZCs(zp),sclMnZ(z));//cumulative pr to sclZCs(zp)
                }
                //cout<<"cprs indices: "<<cprs.indexmin()<<"  "<<cprs.indexmax()<<endl;
                dvar_vector prs(z,nZBs);
                prs(z,nZBs-1) = first_difference(cprs);
                prs(nZBs) = 1.0 - cprs(nZBs);//treat final size bin as accumulator
                //cout<<"prs indices: "<<prs.indexmin()<<"  "<<prs.indexmax()<<endl;
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }            
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else {
            cout<<"Unrecognized growth option: "<<optsGrowth<<endl;
            cout<<"Terminating!"<<endl;
            exit(-1);
        }
        
        prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., "to") z's so n+ = prGr_zz*n
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrGrI->getModelIndices(pc);
        if (debug) cout<<"growth indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1); //year index
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2); //sex index
                for (int s=1;s<=nSCs;s++){
                    mnGrZ_yxsz(y,x,s) = mnGrZ_cz(pc);
                    for (int z=1;z<=nZBs;z++) prGr_yxszz(y,x,s,z) = prGr_czz(pc,z);
                }//s
            }
        }//idx
    }
    
    if (debug>dbgCalcProcs) cout<<"finished calcGrowth()"<<endl;
}

void model_parameters::calcSelectivities(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSelectivities()"<<endl;
    
    SelectivityInfo* ptrSel = ptrMPI->ptrSel;
    double fsZ;             //fully selected size
    int idSel;              //selectivity function id
    int idxFSZ = 1;//index for fsZ in pXDs vector below
    
    ivector mniSelDevs(1,6);//min indices of devs vectors
    ivector mxiSelDevs(1,6);//max indices of devs vectors
    dvar_vector params(1,6);//vector for number_vector params
        
    sel_cz.initialize();//selectivities w/out deviations
    sel_iyz.initialize();//selectivity array
    int y;
    for (int pc=1;pc<=ptrSel->nPCs;pc++){
        params.initialize();
        ivector pids = ptrSel->getPCIDs(pc);
        dvector pXDs = ptrSel->getPCXDs(pc);
        //extract the number parameters
        int k=ptrSel->nIVs+1;//1st parameter variable column
        if (pids[k]) {params[1] = pS1(pids[k]);}   k++;
        if (pids[k]) {params[2] = pS2(pids[k]);}   k++;
        if (pids[k]) {params[3] = pS3(pids[k]);}   k++;
        if (pids[k]) {params[4] = pS4(pids[k]);}   k++;
        if (pids[k]) {params[5] = pS5(pids[k]);}   k++;
        if (pids[k]) {params[6] = pS6(pids[k]);}   k++;
        if (debug>dbgCalcProcs) {
            cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
            cout<<tb<<"params:"<<tb<<params<<endl;
        }
        
        int useDevsS1=pids[k++];
        dvar_vector dvsS1; ivector idxDevsS1;
        if (useDevsS1){
            dvsS1 = devsS1(useDevsS1);
            idxDevsS1 = idxsDevsS1(useDevsS1);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS1) = "<<idxDevsS1<<endl;
                cout<<"dvsS1      = "<<dvsS1<<endl;
            }
        }
        int useDevsS2=pids[k++];
        dvar_vector dvsS2; ivector idxDevsS2;
        if (useDevsS2){
            dvsS2 = devsS2(useDevsS2);
            idxDevsS2 = idxsDevsS2(useDevsS2);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS2) = "<<idxDevsS2<<endl;
                cout<<"dvsS2      = "<<dvsS2<<endl;
            }
        }
        int useDevsS3=pids[k++];
        dvar_vector dvsS3; ivector idxDevsS3;
        if (useDevsS3){
            dvsS3 = devsS3(useDevsS3);
            idxDevsS3 = idxsDevsS3(useDevsS3);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS3) = "<<idxDevsS3<<endl;
                cout<<"dvsS3      = "<<dvsS3<<endl;
            }
        }
        int useDevsS4=pids[k++];
        dvar_vector dvsS4; ivector idxDevsS4;
        if (useDevsS4){
            dvsS4 = devsS4(useDevsS4);
            idxDevsS4 = idxsDevsS4(useDevsS4);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS4) = "<<idxDevsS4<<endl;
                cout<<"dvsS4      = "<<dvsS4<<endl;
            }
        }
        int useDevsS5=pids[k++];
        dvar_vector dvsS5; ivector idxDevsS5;
        if (useDevsS5){
            dvsS5 = devsS5(useDevsS5);
            idxDevsS5 = idxsDevsS5(useDevsS5);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS5) = "<<idxDevsS5<<endl;
                cout<<"dvsS5      = "<<dvsS5<<endl;
            }
        }
        int useDevsS6=pids[k++];
        dvar_vector dvsS6; ivector idxDevsS6;
        if (useDevsS6){
            dvsS6 = devsS6(useDevsS6);
            idxDevsS6 = idxsDevsS6(useDevsS6);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS6) = "<<idxDevsS6<<endl;
                cout<<"dvsS6      = "<<dvsS6<<endl;
            }
        }
        fsZ   = pXDs[idxFSZ];
        idSel = pids[ptrSel->nIVs+ptrSel->nPVs+idxFSZ+1];
        if (debug>dbgCalcProcs) cout<<tb<<"fsZ: "<<fsZ<<tb<<"idSel"<<tb<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<endl;
        sel_cz(pc) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
        if (debug>dbgCalcProcs) cout<<tb<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSel->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//year
            if ((mnYr<=y)&&(y<=mxYr+1)){
                k=ptrSel->nIVs+1+6;//1st devs vector variable column
                if (useDevsS1){
                    if (idxDevsS1[y]){
                        if (debug>dbgCalcProcs) cout<<tb<<idx<<tb<<y<<tb<<useDevsS1<<tb<<idxDevsS1[y]<<endl;
                        params[1] += devsS1(useDevsS1,idxDevsS1[y]);
                    }
                }
                if (useDevsS2) if (idxDevsS2[y]){params[2] += devsS2(useDevsS2,idxDevsS2[y]);}
                if (useDevsS3) if (idxDevsS3[y]){params[3] += devsS3(useDevsS3,idxDevsS3[y]);}
                if (useDevsS4) if (idxDevsS4[y]){params[4] += devsS4(useDevsS4,idxDevsS4[y]);}
                if (useDevsS5) if (idxDevsS5[y]){params[5] += devsS5(useDevsS5,idxDevsS5[y]);}
                if (useDevsS6) if (idxDevsS6[y]){params[6] += devsS6(useDevsS6,idxDevsS6[y]);}
                sel_iyz(pc,y) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"sel: "<<sel_iyz(pc,y)<<endl;
            } else {
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
            }
        }//idx
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished calcSelectivities()"<<endl;
}

void model_parameters::calcFisheryFs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcFisheryFs()"<<endl;
    
    FisheriesInfo* ptrFsh = ptrMPI->ptrFsh;
    
    dvariable hm;                   //handling mortality
    dvar_matrix lnC(1,nSXs,1,nMSs); //ln-scale capture rate
    dvar_matrix C_xm(1,nSXs,1,nMSs);//arithmetic-scale capture rate
    
    dvsLnC_fy.initialize();
    for (int f=1;f<=nFsh;f++) idxDevsLnC_fy(f) = -1;
    
    hasF_fy.initialize();   //flags indicating whether or not fishery occurs
    hmF_fy.initialize();    //handling mortality
    cpF_fyxms.initialize(); //fully-selected capture rate
    sel_fyxmsz.initialize();//selectivity functions
    ret_fyxmsz.initialize();//retention functions
    cpF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retention rate
    dmF_fyxmsz.initialize();//discard mortality rate
    tmF_yxmsz.initialize(); //total mortality rate
    
    /******************************************************\n
     * Fully-selected annual capture rates are calculated  \n
     * using 2 approaches:                                 \n
     *  1. from parameter values (if useER=0 below)        \n
     *  2. based on effort and the average ratio between   \n
     *     effort and capture rates over some time period  \n
     *     in the fishery        (if useER=1 below)        \n
     * Consequently, calculating all the above quantities  \n
     * requires 2 passes through parameter combinations.   \n
    ******************************************************/
    int idxER = ptrFsh->idxUseER;//index into pids below for flag to use effort ratio
    int y; int f; int x; int idSel; int idRet; int useER; int useDevs;
    //Pass 1: calculations based on parameter values
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids: "<<pids<<endl;
        useER = pids[idxER];//flag to use effort ratio
        if (!useER){//calculate capture rates from parameters
            lnC.initialize();
            C_xm.initialize();
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);}                                   k++;
            //set base (ln-scale) capture rate (mature males)
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x) += pLnC(pids[k]);}   k++;
            //add in main temporal offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x) += pLnDCT(pids[k]);} k++;
            if (FEMALE<=nSXs){
                //add in female offset
                if (pids[k]) {lnC(FEMALE) += pLnDCX(pids[k]);}                      k++;
                //add in immature offsets
                if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x,IMMATURE) += pLnDCM(pids[k]);} k++;
                //add in offset immature females for stanza
                if (pids[k]) {lnC(FEMALE,IMMATURE) += pLnDCXM(pids[k]);}            k++; 
            }
            //extract devs vector
            useDevs = pids[k]; k++;
            dvar_vector dvsLnC;             
            ivector idxDevsLnC;
            if (useDevs) {
                dvsLnC     = devsLnC(useDevs);
                idxDevsLnC = idxsDevsLnC(useDevs);
                if (debug>dbgCalcProcs){
                    cout<<"y   idx    devsLnC"<<endl;
                    for (int i=idxDevsLnC.indexmin();i<=idxDevsLnC.indexmax();i++) {
                        cout<<i<<tb<<idxDevsLnC(i)<<tb;
                        if (idxDevsLnC(i)) cout<<dvsLnC[idxDevsLnC(i)];
                        cout<<endl;
                    }
                }
            } else {
                C_xm = exp(lnC);
            }
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];//selectivity function id
            idRet = pids[k++];//retention function id
        
            //convert from ln-scale to arithmetic scale
            if (debug>dbgCalcProcs){
                cout<<"pc: "<<pc<<". idSel = "<<idSel<<". idRet = "<<idRet<<"."<<endl;
                cout<<tb<<tb<<"lnC:"<<endl<<lnC<<endl;
                if (useDevs) {
                    cout<<tb<<tb<<"dvsLnC["<<dvsLnC.indexmin()<<cc<<dvsLnC.indexmax()<<"] = "<<dvsLnC<<endl;
                } else {
                    cout<<tb<<tb<<"C_xm:"<<endl<<C_xm<<endl;
                }
            }
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    x = idxs(idx,3);//sex
                    if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                    if (useDevs) {
                        idxDevsLnC_fy(f,y) = idxDevsLnC[y];
                        dvsLnC_fy(f,y)     = dvsLnC[idxDevsLnC[y]];
                        C_xm = exp(lnC+dvsLnC[idxDevsLnC[y]]);//recalculate C_xm w/ devs
                    }
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            cpF_fyxms(f,y,x,m,s)  = C_xm(x,m);                 //fully-selected capture rate
                            sel_fyxmsz(f,y,x,m,s) = sel_iyz(idSel,y);          //selectivity
                            cpF_fyxmsz(f,y,x,m,s) = C_xm(x,m)*sel_iyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                ret_fyxmsz(f,y,x,m,s) = sel_iyz(idRet,y);      //retention curves
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         cpF_fyxmsz(f,y,x,m,s));//retention mortality
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),cpF_fyxmsz(f,y,x,m,s));//discard mortality
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality
                            }
                        }
                    }
                }
            }
        }//useER=FALSE
    }
    if (debug) cout<<"finished pass 1"<<endl;
    
    //calculate ratio of average capture rate to effort
    if (debug>dbgCalcProcs) cout<<"calculating avgRatioFc2Eff"<<endl;
    dvariable tot;
    for (int f=1;f<=nFsh;f++){//model fishery objects
        int fd = mapM2DFsh(f);//index of corresponding fishery data object
        if (ptrMDS->ppFsh[fd-1]->ptrEff){
            IndexRange* pir = ptrMDS->ppFsh[fd-1]->ptrEff->ptrAvgIR;
            int mny = max(mnYr,pir->getMin());//adjust for model min
            int mxy = min(mxYr,pir->getMax());//adjust for model max
            if (debug>dbgCalcProcs) cout<<"f,mny,mxy = "<<f<<tb<<mny<<tb<<mxy<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        tot.initialize();
                        cpF_fxmsy(f,x,m,s).deallocate();
                        cpF_fxmsy(f,x,m,s).allocate(mny,mxy);
                        switch (optsFcAvg(f)){
                            case 1:
                                for (int y=mny;y<=mxy;y++) cpF_fxmsy(f,x,m,s,y) = cpF_fyxms(f,y,x,m,s);
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/(mxy-mny+1); break;
                            case 2:
                                for (int y=mny;y<=mxy;y++) cpF_fxmsy(f,x,m,s,y) = 1.0-exp(-cpF_fyxms(f,y,x,m,s));
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/(mxy-mny+1); break;
                            case 3:
                                for (int y=mny;y<=mxy;y++) cpF_fxmsy(f,x,m,s,y) = mean(cpF_fyxmsz(f,y,x,m,s));
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/(mxy-mny+1); break;
                            default:
                                cout<<"optsFcAvg("<<f<<") = "<<optsFcAvg(f)<<" to calculate average Fc is invalid."<<endl;
                                cout<<"Aborting..."<<endl;
                                exit(-1);
                        }
                        avgRatioFc2Eff(f,x,m,s) = avgFc_fxms(f,x,m,s)/avgEff(f);
                    }
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"calculated avgRatioFc2Eff"<<endl;
    
    //Pass 2: calculations based on effort and average effort:capture rate ratios
    int fd; double eff;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        useER = pids[idxER];//flag to use effort ratio
        if (useER){//calculate capture rates from parameters
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);} k++;
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];   //selectivity function id
            idRet = pids[k++];   //retention function id
            
            if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<". hm = "<<hm<<". idSel = "<<idSel<<". idRet = "<<idRet<<". Using ER"<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    x = idxs(idx,3);//sex
                    fd = mapM2DFsh(f);//index of corresponding fishery data object
                    eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                    if (debug>dbgCalcProcs) cout<<"f,y,x,eff = "<<f<<tb<<y<<tb<<x<<tb<<eff<<endl;
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            //fully-selected capture rate
                            switch(optsFcAvg(f)) {
                                case 1:
                                    cpF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                                case 2:
                                    cpF_fyxms(f,y,x,m,s) = -log(1.0-avgRatioFc2Eff(f,x,m,s)*eff); break;
                                case 3:
                                    cpF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                            }
                            cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_iyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                            }
                        }
                    }
                }
            }
        }//useER=TRUE
    }
    if (debug>dbgCalcProcs) cout<<"finished pass 2."<<endl;
    if (debug>dbgCalcProcs) cout<<"finished calcFisheryFs()"<<endl;
}

void model_parameters::calcSurveyQs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSurveyQs()"<<endl;
    
    SurveysInfo* ptrSrv = ptrMPI->ptrSrv;
    
    dvar_matrix lnQ(1,nSXs,1,nMSs);
    dvar_matrix Q_xm(1,nSXs,1,nMSs);
    
    q_vyxms.initialize();
    q_vyxmsz.initialize();
    int y; int v; int x; int idSel;
    for (int pc=1;pc<=ptrSrv->nPCs;pc++){
        lnQ.initialize();
        Q_xm.initialize();
        ivector pids = ptrSrv->getPCIDs(pc);
        int k=ptrSrv->nIVs+1;//1st parameter variable column
        //add in base (ln-scale) catchability (mature males)
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x) += pLnQ(pids[k]);}   k++;
        //add in main temporal offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x) += pLnDQT(pids[k]);} k++;
        if (FEMALE<=nSXs){
            //add in female offset
            if (pids[k]) {lnQ(FEMALE) += pLnDQX(pids[k]);}                      k++;
            //add in immature offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x,IMMATURE) += pLnDQM(pids[k]);} k++;
            //add in offset immature females for stanza
            if (pids[k]) {lnQ(FEMALE,IMMATURE) += pLnDQXM(pids[k]);}            k++; 
        }
        
        idSel = pids[k];//selectivity function id
        
        //convert from ln-scale to arithmetic scale
        Q_xm = exp(lnQ);
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"lnQ:"<<endl<<lnQ<<endl;
            cout<<"pc: "<<pc<<tb<<"Q_xm:"<<endl<<Q_xm<<endl;
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSrv->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            v = idxs(idx,1);//survey
            y = idxs(idx,2);//year
            if ((mnYr<=y)&&(y<=mxYrp1)){
                x = idxs(idx,3);//sex
                if (y <= mxYr+1){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            q_vyxms(v,y,x,m,s)  = Q_xm(x,m);
                            s_vyxmsz(v,y,x,m,s) = sel_iyz(idSel,y);
                            q_vyxmsz(v,y,x,m,s) = Q_xm(x,m)*sel_iyz(idSel,y);
                        }
                    }
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
    
}

void model_parameters::calcPenalties(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    if (debug<0) cout<<"list("<<endl;//start list of penalties by category
    if (debug<0) cout<<tb<<"maturity=list("<<endl;//start of maturity penalties list
    //smoothness penalties on maturity parameters (NOT maturity ogives)
    double penWgtSmthLgtPrMat = ptrMOs->wgtSmthLgtPrMat;
    fPenSmoothLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = 1.0*pLgtPrMat(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i)<<"),"<<endl;
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = 1.0*pLgtPrMat(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i)<<")"<<endl;
    }
    if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
    //non-decreasing penalties on maturity parameters (NOT maturity ogives)
    double penWgtNonDecLgtPrMat = ptrMOs->wgtNonDecLgtPrMat;
    fPenNonDecLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"nondecreasing=list(";//start of non-decreasing penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = calc1stDiffs(pLgtPrMat(i));
        for (int iv=v.indexmin();iv<=v.indexmax();iv++){
            posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<"),";
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = calc1stDiffs(pLgtPrMat(i));
        for (int iv=v.indexmin();iv<=v.indexmax();iv++){
            posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<")";
    }
    if (debug<0) cout<<tb<<tb<<")"<<endl;//end of non-decreasing penalties list    
    if (debug<0) cout<<tb<<"),";//end of maturity penalties list
    
    //penalties on final value of dev vectors
    double penWgt = 0.0;
    if (current_phase()>=ptrMOs->phsLastDevsPen) penWgt = ptrMOs->wgtLastDevsPen;
    if (debug<0) cout<<tb<<"final.devs=list("<<endl;//start of devs penalties list
    //recruitment devs
    if (ptrMPI->ptrRec->pDevsLnR->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsLnR=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnR,devsLnR);
        if (debug<0) cout<<cc<<endl;
    }
    //S1 devs
    if (ptrMPI->ptrSel->pDevsS1->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS1=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS1,devsS1);
        if (debug<0) cout<<cc<<endl;
    }
    //S2 devs
    if (ptrMPI->ptrSel->pDevsS2->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS2=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS2,devsS2);
        if (debug<0) cout<<cc<<endl;
    }
    //S3 devs
    if (ptrMPI->ptrSel->pDevsS3->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS3=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS3,devsS3);
        if (debug<0) cout<<cc<<endl;
    }
    //S4 devs
    if (ptrMPI->ptrSel->pDevsS4->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS4=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS4,devsS4);
        if (debug<0) cout<<cc<<endl;
    }
    //S5 devs
    if (ptrMPI->ptrSel->pDevsS5->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS5=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS5,devsS5);
        if (debug<0) cout<<cc<<endl;
    }
    //S6 devs
    if (ptrMPI->ptrSel->pDevsS6->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS6=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS6,devsS6);
        if (debug<0) cout<<cc<<endl;
    }
    //capture rate devs
    if (ptrMPI->ptrFsh->pDevsLnC->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsLnC=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnC,devsLnC);
        if (debug<0) cout<<cc<<endl;
    }
    
    if (debug<0) cout<<tb<<"NULL)"<<endl;//end of devs penalties list
    if (debug<0) cout<<")"<<cc<<endl;//end of penalties list
    
    //Apply penalties to F-devs
    {
        if (debug<0) cout<<"penFDevs=list("<<endl;
        double penWgt = 1.0/log(1.0+square(ptrMOs->cvFDevsPen));
        if (ptrMOs->phsDecrFDevsPen<=current_phase()){
            double scl = max(1.0,(double)(ptrMOs->phsZeroFDevsPen-ptrMOs->phsDecrFDevsPen));
            penWgt *= max(0.0,(double)(ptrMOs->phsZeroFDevsPen-current_phase()))/scl;
        }
        double effCV = std::numeric_limits<double>::infinity();
        if (penWgt>0) {
            effCV = sqrt(exp(1.0/penWgt)-1.0);
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = "<<effCV<<endl;
        } else {
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = Inf"<<endl;
        }
        for (int i=pDevsLnC.indexmin();i<=pDevsLnC.indexmax();i++){
            dvariable fpen = 0.5*norm2(pDevsLnC(i));
            objFun += penWgt*fpen;
            if (debug<0) {
                rpt::echo<<tb<<i<<": pen="<<fpen<<cc<<" objfun="<<penWgt*fpen<<endl;
                if (penWgt>0) {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV="<<effCV<<cc<<"pen="<<fpen<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                } else {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV=Inf"    <<cc<<"pen="<<fpen<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                }
            }
        }
        if (debug<0) cout<<tb<<tb<<"NULL)"<<endl;//end of penFDevs lists
    }
    
    if (debug>=dbgObjFun) cout<<"Finished calcPenalties()"<<endl;
}

void model_parameters::calcDevsPenalties(int debug, ostream& cout, double penWgt, param_init_bounded_vector_vector& pDevs, dvar_matrix devs)
{
    double scale = 1.0e-2;//scale for posfun calculation
    int idx;
    double lower; double upper;
    dvariable fPenLower;
    dvariable fPenUpper;
    if (debug<0) cout<<"list(";//start of list
    for (int i=pDevs.indexmin();i<=pDevs.indexmax();i++){
        if (pDevs(i).get_phase_start()){
            idx = devs(i).indexmax();//index of final dev
            lower = pDevs(i).get_minb();
            upper = pDevs(i).get_maxb();
            fPenLower.initialize();
            fPenUpper.initialize();
            posfun2(devs(i,idx)-lower,scale,fPenLower);
            posfun2(upper-devs(i,idx),scale,fPenUpper);
            objFun += penWgt*(fPenLower+fPenUpper);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"pen="<<fPenLower+fPenUpper<<cc<<"objfun="<<penWgt*(fPenLower+fPenUpper)<<cc<<
                                                          "val="<<devs(i,idx)<<cc<<"minb="<<lower<<cc<<"maxb="<<upper<<cc<<"fPenL="<<fPenLower<<cc<<"fPenU="<<fPenUpper<<"),"<<endl;
        }
    }
    if (debug<0) cout<<tb<<tb<<"NULL)";//end of penalties list
    
}

dvar_vector model_parameters::calc1stDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector cp; cp = d;
    dvar_vector r(mn,mx-1);
    r = cp(mn+1,mx).shift(mn)-cp(mn,mx-1);
    RETURN_ARRAYS_DECREMENT();
    return r;
}

dvar_vector model_parameters::calc2ndDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector r = calc1stDiffs(calc1stDiffs(d));
    RETURN_ARRAYS_DECREMENT();
    return r;
}

void model_parameters::calcNLLs_Recruitment(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Starting calcNLLs_Recruitment"<<endl;
    double nllWgtRecDevs = 1.0;//TODO: read in from input file (as vector?))
    nllRecDevs.initialize();
    if (debug<0) cout<<"list("<<endl;
    if (debug<0) cout<<tb<<"recDevs=list("<<endl;
    for (int pc=1;pc<npcRec;pc++){
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<cc<<endl;
        }
    }//pc
    {
        int pc = npcRec;
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<endl;
        }
    }//pc
    if (debug<0) cout<<tb<<")";//recDevs
    if (debug<0) cout<<")";
    if (debug>=dbgObjFun) cout<<"Finished calcNLLs_Recruitment"<<endl;
}

void model_parameters::calcObjFun(int debug, ostream& cout)
{
    if ((debug>=dbgObjFun)||(debug<0)) cout<<"Starting calcObjFun"<<endl;
    //objective function penalties
    calcPenalties(debug,cout);
    //prior likelihoods
    calcAllPriors(debug,cout);
    //recruitment component
    calcNLLs_Recruitment(debug,cout);
    
    //data components
    calcNLLs_Fisheries(debug,cout);
    calcNLLs_Surveys(debug,cout);
    
    if (debug<0) cout<<"total objFun = "<<objFun<<endl;
    if ((debug>=dbgObjFun)||(debug<0)) cout<<"Finished calcObjFun"<<endl<<endl;
    
}

void model_parameters::calcNorm2NLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNorm2NLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y]);
        }
    }
    nll += 0.5*norm2(zscr);
    double wgt = 1.0;//TODO: implement likelihood weights
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='norm2',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNorm2NLL()"<<endl;
    
}

void model_parameters::calcNormalNLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
   if (debug>=dbgAll) cout<<"Starting calcNormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (obs[i]-mod[y])/stdv[i];
                //  nll += log(stdv[i]);
            }
        }
    }
    nll += 0.5*norm2(zscr);
    double wgt = 1.0;//TODO: implement likelihood weights
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcNormalNLL()"<<endl;
    
}

void model_parameters::calcLognormalNLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcLognormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (log(obs[i]+smlVal)-log(mod[y]+smlVal))/stdv[i];
                //nll += log(stdv[i]);
            }
        }
    }
    nll += 0.5*norm2(zscr);
    double wgt = 1.0;//TODO: implement likelihood weights
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='lognormal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcLognormalNLL()"<<endl;
    
}

void model_parameters::calcMultinomialNLL(dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcMultinomialNLL()"<<endl;
    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
    double wgt = 1.0;//TODO: incorporate weights
    objFun += wgt*nll;
    if (debug<0){
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = 0.0;
        if (ss>0) effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(nll.type='multinomial',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
        cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcMultinomialNLL()"<<endl;
 
}

void model_parameters::calcNLL(int llType, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    switch (llType){
        case tcsam::LL_NONE:
            break;
        case tcsam::LL_LOGNORMAL:
            calcLognormalNLL(mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORMAL:
            calcNormalNLL(mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORM2:
            calcNorm2NLL(mod,obs,stdv,yrs,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(1)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }    
}

void model_parameters::calcNLL(int llType, dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
{
    switch (llType){
        case tcsam::LL_NONE:
            break;
        case tcsam::LL_MULTINOMIAL:
            calcMultinomialNLL(mod,obs,ss,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(2)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }
   
}

void model_parameters::calcNLLs_AggregateCatch(AggregateCatchData* ptrAB, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_AggregateCatch()"<<endl;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvar_vector tAB_y(mny,mxy);
    int isBio = ptrAB->type==AggregateCatchData::KW_BIOMASS_DATA;
    if (debug>=dbgAll) cout<<"isBio="<<isBio<<tb<<"type="<<ptrAB->type<<endl;
    if (debug<0) cout<<"list(fit.type='"<<tcsam::getFitType(ptrAB->optFit)<<"',fits=list("<<endl;
    if (ptrAB->optFit==tcsam::FIT_BY_TOT){
        tAB_y.initialize();
        if (isBio){
            for (int x=1;x<=nSXs;x++) {
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            }//x
        } else {
            for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y));//sum over x,m,s,z
        }
        if (debug>=dbgAll) cout<<"FIT_BY_TOT: "<<tAB_y<<endl;
        if (debug<0) {
            cout<<"list(";
            cout<<"x="<<qt<<tcsam::getSexType(ALL_SXs)     <<qt<<cc;
            cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
            cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
            cout<<"nll=";
        }
        calcNLL(ptrAB->llType, tAB_y, ptrAB->C_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout);                
        if (debug<0) cout<<")";
        if (debug<0) cout<<")";
    } else if (ptrAB->optFit==tcsam::FIT_BY_X){
        for (int x=1;x<=nSXs;x++){
            tAB_y.initialize();
            if (isBio){
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            } else {
                for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x));//sum over m,s,z
            }
            if (debug>=dbgAll) cout<<"FIT_BY_X("<<x<<"): "<<tAB_y<<endl;
            if (debug<0) {
                cout<<"list(";
                cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
                cout<<"nll=";
            }
            calcNLL(ptrAB->llType, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout); 
            if (debug<0) cout<<"),"<<endl;
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XM){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                tAB_y.initialize();
                if (isBio){
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                } else {
                    for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m));//sum over s,z
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)        <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(m)   <<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)<<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, tAB_y, ptrAB->C_xmsy(x,m,ALL_SCs), ptrAB->sd_xmsy(x,m,ALL_SCs), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XS){
        for (int x=1;x<=nSXs;x++){
            for (int s=1;s<=nSCs;s++){
                tAB_y.initialize();
                if (isBio){
                    for (int m=1;m<=nMSs;m++) {
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//m
                } else {
                    for (int m=1;m<=nMSs;m++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over m,z
                        }//y
                    }//m
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(s)         <<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,s), ptrAB->sd_xmsy(x,ALL_MSs,s), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//s
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XMS){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    tAB_y.initialize();
                    if (isBio){
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                    } else {
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over z
                    }
                    if (debug<0) {
                        cout<<"list(";
                        cout<<"x="<<qt<<tcsam::getSexType(x)     <<qt<<cc;
                        cout<<"m="<<qt<<tcsam::getMaturityType(m)<<qt<<cc;
                        cout<<"s="<<qt<<tcsam::getShellType(s)   <<qt<<cc;
                        cout<<"nll=";
                    }
                    calcNLL(ptrAB->llType, tAB_y, ptrAB->C_xmsy(x,m,s), ptrAB->sd_xmsy(x,m,s), ptrAB->yrs, debug, cout); 
                    if (debug<0) cout<<"),"<<endl;
                }//s
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else {
        std::cout<<"Calling calcNLLs_AggregateCatch with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrAB->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<")";
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_AggregateCatch()"<<endl;
    }
}

void model_parameters::calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    if (debug>=dbgNLLs) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    if (ptrZFD->optFit==tcsam::FIT_NONE) return;
    ivector yrs = ptrZFD->yrs;
    int y;
    double ss;
    dvariable nT;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvector     oP_z;//observed size comp.
    dvar_vector mP_z;//model size comp.
    if (ptrZFD->optFit==tcsam::FIT_BY_XE){
        oP_z.allocate(1,nSXs*nZBs);
        mP_z.allocate(1,nSXs*nZBs);
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
        oP_z.allocate(1,nMSs*nZBs);
        mP_z.allocate(1,nMSs*nZBs);
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_XME){
        oP_z.allocate(1,nSXs*nMSs*nZBs);
        mP_z.allocate(1,nSXs*nMSs*nZBs);
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
    } else 
    {
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
    }
    if (debug<0) cout<<"list("<<endl;
    for (int iy=1;iy<=yrs.size();iy++) {
        y = yrs[iy];
        if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
        if ((mny<=y)&&(y<=mxy)) {
            if (ptrZFD->optFit==tcsam::FIT_BY_TOT){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=ALL_SXs;x++){
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                    }
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++)mP_z += mA_yxmsz(y,x,m,s);
                        }
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    if (debug<0) {
                        cout<<"'"<<y<<"'=list(";
                        cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                        cout<<"y="<<y<<cc;
                        cout<<"x='"<<tcsam::getSexType(ALL_SXs)<<"'"<<cc;
                        cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                        cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                        cout<<"fit=";
                    }
                    calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                    if (debug<0) cout<<")"<<cc<<endl;
                }
                //FIT_BY_TOT
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
                //FIT_BY_X
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XM){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XM
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XS){
                for (int x=1;x<=nSXs;x++) {
                    for (int s=1;s<=nSCs;s++){
                        ss = 0;
                        nT.initialize();
                        for (int m=1;m<=nMSs;m++) nT += sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int m=1;m<=nMSs;m++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XS
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XMS){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) {
                            ss = 0;
                            nT = sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                            if (value(nT)>0){
                                oP_z.initialize();//observed size comp.
                                mP_z.initialize();//model size comp.                            
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                if (sum(oP_z)>0) oP_z /= sum(oP_z);
                                if (debug>=dbgNLLs){
                                    cout<<"ss = "<<ss<<endl;
                                    cout<<"oP_Z = "<<oP_z<<endl;
                                }
                                mP_z += mA_yxmsz(y,x,m,s);
                                mP_z /= nT;//normalize model size comp
                                if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//nT>0
                        }//s
                    }//m
                }//x
                //FIT_BY_XMS
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XE){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                        }
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        dvar_vector mPt = mP_z(mnz,mxz);
                        dvector oPt = oP_z(mnz,mxz);
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                }//nT>0
                //FIT_BY_XE
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//s
                            for (int s=1;s<=nSCs;s++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                        }//m
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//nT>0
                }//x
                //FIT_BY_X_ME
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//m
                            for (int m=1;m<=nMSs;m++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//m
                        }//s
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//s
                    }//nT>0
                }//x
                //FIT_BY_X_SE
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_XME){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            if (m<=nMSs) {for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);}
                        }//m
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                }//nT>0
                //FIT_BY_XME
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++) {
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                dvar_vector mPt = mP_z(mnz,mxz);
                                dvector oPt = oP_z(mnz,mxz);
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XM_SE
            } else 
            {
                std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
                std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
                std::cout<<"Aborting..."<<endl;
                exit(-1);
            }
        } //if ((mny<=y)&&(y<=mxy))
    } //loop over iy
    if (debug<0) cout<<"NULL)";
    if (debug>=dbgNLLs) cout<<"Finished calcNLLs_CatchNatZ()"<<endl;
}

void model_parameters::calcNLLs_Fisheries(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Fisheries()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug>=dbgAll) cout<<"calculating NLLs for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<ptrMC->lblsFsh[f]<<"=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppFsh[f-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasN && ptrObs->ptrRCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---retained catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrN,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrB,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasN && ptrObs->ptrTCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---total catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrN,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrB,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasN && ptrObs->ptrDCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---discard catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrN,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrB,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Fisheries()"<<endl;
}

void model_parameters::calcNLLs_Surveys(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Surveys()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug>=dbgAll) cout<<"calculating NLLs for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<ptrMC->lblsSrv[v]<<"=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppSrv[v-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasN && ptrObs->ptrICD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---index catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrN,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasB && ptrObs->ptrICD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---index catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrB,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---index catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrICD->ptrZFD,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Surveys()"<<endl;
}

void model_parameters::calcAllPriors(int debug, ostream& cout)
{
    if (debug>=dbgPriors) cout<<"Starting calcAllPriors()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    //recruitment parameters
    if (debug<0) cout<<tb<<"recruitment=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,  pLnR,  debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRCV,pLnRCV,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLgtRX,pLgtRX,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRa, pLnRa, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRb, pLnRb, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //natural mortality parameters
    if (debug<0) cout<<tb<<"'natural mortality'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnM,   pLnM,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMT, pLnDMT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMX, pLnDMX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMM, pLnDMM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMXM,pLnDMXM,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //growth parameters
    if (debug<0) cout<<tb<<"growth=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrA,   pLnGrA,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrB,   pLnGrB,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrBeta,pLnGrBeta,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //maturity parameters
    if (debug<0) cout<<tb<<"maturity=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrMat->pLgtPrMat,pLgtPrMat,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //selectivity parameters
    if (debug<0) cout<<tb<<"'selectivity functions'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS1,pS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS2,pS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS3,pS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS4,pS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS5,pS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS6,pS6,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS1,devsS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS2,devsS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS3,devsS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS4,devsS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS5,devsS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS6,devsS6,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //fishing mortality parameters
    if (debug<0) cout<<tb<<"fisheries=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC,    pLnC,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCT,  pLnDCT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCX,  pLnDCX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCM,  pLnDCM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCXM, pLnDCXM,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //survey catchability parameters
    if (debug<0) cout<<tb<<"surveys=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnQ,    pLnQ,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQT,  pLnDQT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQX,  pLnDQX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQM,  pLnDQM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQXM, pLnDQXM,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<endl;
    
    if (debug<0) cout<<")"<<endl;
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;
}

void model_parameters::ReportToR_Data(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Data(...)"<<endl;
    ptrMDS->writeToR(os,"data",0);
    if (debug) cout<<"Finished ReportToR_Data(...)"<<endl;
}

void model_parameters::ReportToR_Params(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Params(...)"<<endl;
    ptrMPI->writeToR(os);
    if (debug) cout<<"Finished ReportToR_Params(...)"<<endl;
}

void model_parameters::ReportToR_ModelProcesses(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelProcesses(...)"<<endl;
    
    os<<"mp=list("<<endl;
        os<<"M_cxm         ="; wts::writeToR(os,value(M_cxm),   adstring("pc=1:"+str(npcNM )),xDms,mDms);   os<<cc<<endl;
        os<<"prMolt2Mat_cz ="; wts::writeToR(os,value(prMat_cz),adstring("pc=1:"+str(npcMat)),zbDms);       os<<cc<<endl;
        os<<"sel_cz        ="; wts::writeToR(os,value(sel_cz),  adstring("pc=1:"+str(npcSel)),zbDms);       os<<cc<<endl;
        
        os<<"M_yxmsz        =";  wts::writeToR(os,wts::value(M_yxmsz),   yDms,xDms,mDms,sDms,zbDms);  os<<cc<<endl;
        os<<"prMolt2Mat_yxz =";  wts::writeToR(os,     value(prMat_yxz), yDms,xDms,zbDms);            os<<cc<<endl;
        os<<"T_list=list("<<endl;
            os<<"mnZAM_cz   ="; wts::writeToR(os,value(mnGrZ_cz),adstring("pc=1:"+str(npcGr )),zbDms);       os<<cc<<endl;
            os<<"T_czz      ="; wts::writeToR(os,value(prGr_czz),adstring("pc=1:"+str(npcGr )),zbDms,zpDms); os<<cc<<endl;
            os<<"mnZAM_yxsz =";  wts::writeToR(os,wts::value(mnGrZ_yxsz),yDms,xDms,sDms,zbDms);              os<<cc<<endl;
            os<<"T_yxszz    =";  wts::writeToR(os,wts::value(prGr_yxszz),yDms,xDms,sDms,zbDms,zpDms);        os<<endl;
        os<<")"<<cc<<endl;
        os<<"R_list=list("<<endl;
            os<<"R_y  ="; wts::writeToR(os,value(R_y), yDms);                                os<<cc<<endl;
            os<<"R_yx ="; wts::writeToR(os,value(R_yx),yDms,xDms);                           os<<cc<<endl;
            os<<"R_yz ="; wts::writeToR(os,value(R_yz),yDms,zbDms);                          os<<cc<<endl;
            os<<"Rx_c ="; wts::writeToR(os,value(Rx_c),adstring("pc=1:"+str(npcRec)));       os<<cc<<endl;
            os<<"R_cz ="; wts::writeToR(os,value(R_cz),adstring("pc=1:"+str(npcRec)),zbDms); os<<endl;
        os<<")"<<cc<<endl;
        d6_array tmF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
        for (int f=1;f<=nFsh;f++){
            for (int y=mnYr;y<=mxYr;y++){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            tmF_fyxmsz(f,y,x,m,s) = value(rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s));
                        }
                    }
                }
            }
        }
        os<<"Eff_list=list("<<endl;
        for (int f=1;f<=nFsh;f++){
            int fd = mapM2DFsh(f);//index of corresponding fishery data object
            os<<"`"<<ptrMDS->ppFsh[fd-1]->name<<"`=list("<<endl;
            if (ptrMDS->ppFsh[fd-1]->ptrEff) {
                adstring yDmsp = "y="+str(eff_fy(f).indexmin())+":"+str(eff_fy(f).indexmax());
                os<<"optFcAvg="<<optsFcAvg(f)<<cc;
                os<<"avgEff="<<avgEff(f)<<cc<<endl;
                os<<"eff_y ="; wts::writeToR(os,eff_fy(f),yDmsp); os<<cc<<endl;
                os<<"cpF_xmsy ="; wts::writeToR(os,wts::value(cpF_fxmsy(f)), xDms,mDms,sDms,yDmsp); os<<cc<<endl;
                os<<"avgFc_xms ="; wts::writeToR(os,    value(avgFc_fxms(f)),xDms,mDms,sDms); os<<endl;
            }
            os<<")"<<cc<<endl;
        }
        os<<"NULL)"<<cc<<endl;
        os<<"F_list=list("<<endl;
            //handling mortality
            os<<"hm_fy     ="; wts::writeToR(os,     value(hmF_fy),    fDms,yDms);                     os<<cc<<endl;
            os<<"cpF_fyxms ="; wts::writeToR(os,wts::value(cpF_fyxms ),fDms,yDms,xDms,mDms,sDms);      os<<cc<<endl;
            os<<"sel_fyxmsz="; wts::writeToR(os,wts::value(sel_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"ret_fyxmsz="; wts::writeToR(os,wts::value(ret_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_yxmsz ="; wts::writeToR(os,wts::value(tmF_yxmsz),      yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"cpF_fyxmsz="; wts::writeToR(os,wts::value(cpF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmF_fyxmsz="; wts::writeToR(os,wts::value(rmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmF_fyxmsz="; wts::writeToR(os,wts::value(dmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_fyxmsz="; wts::writeToR(os,            tmF_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
            os<<"sel_vyxmsz="; wts::writeToR(os,wts::value(s_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"Q_vyxms   ="; wts::writeToR(os,wts::value(q_vyxms), vDms,ypDms,xDms,mDms,sDms);       os<<cc<<endl;
            os<<"Q_vyxmsz  ="; wts::writeToR(os,wts::value(q_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelProcesses(...)"<<endl;
}

void model_parameters::ReportToR_ModelResults(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelResults(...)"<<endl;
    
    d3_array ones(1,nSXs,1,nMSs,1,nZBs);//weighting for simple sums
    for (int x=1;x<=nSXs;x++) ones(x) = 1.0;
    
    //population numbers, biomass
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    d4_array n_yxms   = tcsam::calcYXMSfromYXMSZ(vn_yxmsz,ones);
    d4_array b_yxms   = tcsam::calcYXMSfromYXMSZ(vn_yxmsz,ptrMDS->ptrBio->wAtZ_xmz);
        
    //numbers, biomass captured (NOT mortality)
    d6_array vcpN_fyxmsz = wts::value(cpN_fyxmsz);
    d5_array cpN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vcpN_fyxmsz,ones);
    d5_array cpB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vcpN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discards (NOT mortality)
    d6_array vdsN_fyxmsz = wts::value(dsN_fyxmsz);
    d5_array dsN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdsN_fyxmsz,ones);
    d5_array dsB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdsN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass retained (mortality)
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    d5_array rmN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vrmN_fyxmsz,ones);
    d5_array rmB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discard mortality
    d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
    d5_array dmN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdmN_fyxmsz,ones);
    d5_array dmB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //survey numbers, biomass
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d5_array n_vyxms   = tcsam::calcIYXMSfromIYXMSZ(vn_vyxmsz,ones);
    d5_array b_vyxms   = tcsam::calcIYXMSfromIYXMSZ(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    os<<"mr=list("<<endl;
        os<<"iN_xmsz ="; wts::writeToR(os,vn_yxmsz(mnYr),xDms,mDms,sDms,zbDms); os<<cc<<endl;
        os<<"P_list=list("<<endl;
            os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms);                       os<<cc<<endl;
            os<<"B_yxms   ="; wts::writeToR(os,       b_yxms,ypDms,xDms,mDms,sDms);             os<<cc<<endl;
            os<<"N_yxmsz  ="; wts::writeToR(os,     vn_yxmsz,ypDms,xDms,mDms,sDms,zbDms);        os<<cc<<endl;
            os<<"nmN_yxmsz="; wts::writeToR(os,wts::value(nmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmN_yxmsz="; wts::writeToR(os,wts::value(tmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;    
        os<<"F_list=list("<<endl;
            os<<"cpB_fyxms ="; wts::writeToR(os,cpB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
            os<<"dsB_fyxms ="; wts::writeToR(os,dsB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
            os<<"rmB_fyxms ="; wts::writeToR(os,rmB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
            os<<"dmB_fyxms ="; wts::writeToR(os,dmB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
            os<<"cpN_fyxmsz="; wts::writeToR(os,vcpN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dsN_fyxmsz="; wts::writeToR(os,vdsN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmN_fyxmsz="; wts::writeToR(os,vrmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmN_fyxmsz="; wts::writeToR(os,vdmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
           os<<"MB_vyx  ="; wts::writeToR(os,value(mb_vyx),vDms,ypDms,xDms);            os<<cc<<endl;
           os<<"B_vyxms ="; wts::writeToR(os,      b_vyxms,vDms,ypDms,xDms,mDms,sDms);  os<<cc<<endl;
           os<<"N_vyxmsz="; wts::writeToR(os,    vn_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<endl;
       os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelResults(...)"<<endl;
    
}

void model_parameters::ReportToR_ModelFits(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    //recalc objective function components and and write results to os
    os<<"model.fits=list("<<endl;
        os<<tb<<"penalties="; calcPenalties(-1,os);      os<<cc<<endl;
        os<<tb<<"priors=";    calcAllPriors(-1,os);      os<<cc<<endl;
        os<<tb<<"components=list("<<endl;
            os<<tb<<tb<<"recruitment="; calcNLLs_Recruitment(-1,os); os<<endl;
        os<<tb<<")"<<cc<<endl;
        os<<tb<<"fisheries="; calcNLLs_Fisheries(-1,os);  os<<cc<<endl; 
        os<<tb<<"surveys=";   calcNLLs_Surveys(-1,os);    os<<endl;  
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelFits(...)"<<endl;
}

void model_parameters::updateMPI(int debug, ostream& cout)
{
    if (debug) cout<<"Starting updateMPI(...)"<<endl;
    //recruitment parameters
    ptrMPI->ptrRec->pLnR->setFinalVals(pLnR);
    ptrMPI->ptrRec->pLnRCV->setFinalVals(pLnRCV);
    ptrMPI->ptrRec->pLgtRX->setFinalVals(pLgtRX);
    ptrMPI->ptrRec->pLnRa->setFinalVals(pLnRa);
    ptrMPI->ptrRec->pLnRb->setFinalVals(pLnRb);
    //cout<<"setting final vals for pDevsLnR"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) (*ptrMPI->ptrRec->pDevsLnR)[p]->setFinalVals(pDevsLnR(p));
     
    //natural mortality parameters
    ptrMPI->ptrNM->pLnM->setFinalVals(pLnM);
    ptrMPI->ptrNM->pLnDMT->setFinalVals(pLnDMT);
    ptrMPI->ptrNM->pLnDMX->setFinalVals(pLnDMX);
    ptrMPI->ptrNM->pLnDMM->setFinalVals(pLnDMM);
    ptrMPI->ptrNM->pLnDMXM->setFinalVals(pLnDMXM);
    
    //growth parameters
    ptrMPI->ptrGr->pLnGrA->setFinalVals(pLnGrA);
    ptrMPI->ptrGr->pLnGrB->setFinalVals(pLnGrB);
    ptrMPI->ptrGr->pLnGrBeta->setFinalVals(pLnGrBeta);
    
    //maturity parameters
    //cout<<"setting final vals for pLgtPrMat"<<endl;
    for (int p=1;p<=npLgtPrMat;p++) (*ptrMPI->ptrMat->pLgtPrMat)[p]->setFinalVals(pLgtPrMat(p));
    
    //selectivity parameters
    ptrMPI->ptrSel->pS1->setFinalVals(pS1);
    ptrMPI->ptrSel->pS2->setFinalVals(pS2);
    ptrMPI->ptrSel->pS3->setFinalVals(pS3);
    ptrMPI->ptrSel->pS4->setFinalVals(pS4);
    ptrMPI->ptrSel->pS5->setFinalVals(pS5);
    ptrMPI->ptrSel->pS6->setFinalVals(pS6);
    //cout<<"setting final vals for pDevsS1"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS1->getSize();p++) (*ptrMPI->ptrSel->pDevsS1)[p]->setFinalVals(pDevsS1(p));
    //cout<<"setting final vals for pDevsS2"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS2->getSize();p++) (*ptrMPI->ptrSel->pDevsS2)[p]->setFinalVals(pDevsS2(p));
    //cout<<"setting final vals for pDevsS3"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS3->getSize();p++) (*ptrMPI->ptrSel->pDevsS3)[p]->setFinalVals(pDevsS3(p));
    //cout<<"setting final vals for pDevsS4"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS4->getSize();p++) (*ptrMPI->ptrSel->pDevsS4)[p]->setFinalVals(pDevsS4(p));
    //cout<<"setting final vals for pDevsS5"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS5->getSize();p++) (*ptrMPI->ptrSel->pDevsS5)[p]->setFinalVals(pDevsS5(p));
    //cout<<"setting final vals for pDevsS6"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS6->getSize();p++) (*ptrMPI->ptrSel->pDevsS6)[p]->setFinalVals(pDevsS6(p));
     
    //fully-selected fishing capture rate parameters
    ptrMPI->ptrFsh->pHM->setFinalVals(pHM);
    ptrMPI->ptrFsh->pLnC->setFinalVals(pLnC);
    ptrMPI->ptrFsh->pLnDCT->setFinalVals(pLnDCT);
    ptrMPI->ptrFsh->pLnDCX->setFinalVals(pLnDCX);
    ptrMPI->ptrFsh->pLnDCM->setFinalVals(pLnDCM);
    ptrMPI->ptrFsh->pLnDCXM->setFinalVals(pLnDCXM);
    //cout<<"setting final vals for pDevsLnC"<<endl;
    for (int p=1;p<=ptrMPI->ptrFsh->pDevsLnC->getSize();p++) (*ptrMPI->ptrFsh->pDevsLnC)[p]->setFinalVals(pDevsLnC(p));
    
    //survey catchability parameters
    ptrMPI->ptrSrv->pLnQ->setFinalVals(pLnQ);
    ptrMPI->ptrSrv->pLnDQT->setFinalVals(pLnDQT);
    ptrMPI->ptrSrv->pLnDQX->setFinalVals(pLnDQX);
    ptrMPI->ptrSrv->pLnDQM->setFinalVals(pLnDQM);
    ptrMPI->ptrSrv->pLnDQXM->setFinalVals(pLnDQXM);
    
    if (debug) cout<<"Finished updateMPI(...)"<<endl;
}

void model_parameters::ReportToR(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR(...)"<<endl;
    updateMPI(debug,cout);
        
    os<<"res=list("<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        
        //model data
        ptrMDS->writeToR(os,"data",0); os<<","<<endl;
        
        //parameter values
        ReportToR_Params(os,debug,cout); os<<","<<endl;
        
        //model processes
        ReportToR_ModelProcesses(os,debug,cout); os<<","<<endl;
        
        //model results
        ReportToR_ModelResults(os,debug,cout); os<<","<<endl;
        //model fit quantities
        ReportToR_ModelFits(os,debug,cout); os<<","<<endl;
        
        //simulated model data
        createSimData(debug, cout, 0, ptrSimMDS);//deterministic
        ptrSimMDS->writeToR(os,"sim.data",0); os<<","<<endl;
        
        //do OFL calculations
        calcOFL(mxYr,0,cout);//updates ptrOFL
        ptrOFL->writeToR(os,"oflResults",0);
    os<<")"<<endl;
    if (debug) cout<<"Finished ReportToR(...)"<<endl;
}

void model_parameters::writeParameters(ofstream& os,int toR, int willBeActive)
{
    os<<"index, phase, idx.mn, idx.mx, min, max, value, name, type"<<endl;
    //recruitment parameters
    wts::writeParameter(os,pLnR,toR,willBeActive);      
    wts::writeParameter(os,pLnRCV,toR,willBeActive);      
    wts::writeParameter(os,pLgtRX,toR,willBeActive);      
    wts::writeParameter(os,pLnRa,toR,willBeActive);      
    wts::writeParameter(os,pLnRb,toR,willBeActive);      
    wts::writeParameter(os,pDevsLnR,toR,willBeActive);      
    
    //natural mortality parameters
    wts::writeParameter(os,pLnM,toR,willBeActive);      
    wts::writeParameter(os,pLnDMT,toR,willBeActive);      
    wts::writeParameter(os,pLnDMX,toR,willBeActive);      
    wts::writeParameter(os,pLnDMM,toR,willBeActive);      
    wts::writeParameter(os,pLnDMXM,toR,willBeActive);      
    
    //growth parameters
    wts::writeParameter(os,pLnGrA,toR,willBeActive);      
    wts::writeParameter(os,pLnGrB,toR,willBeActive);      
    wts::writeParameter(os,pLnGrBeta,toR,willBeActive);      
    
    //maturity parameters
    wts::writeParameter(os,pLgtPrMat,toR,willBeActive);      
    
    //selectivity parameters
    wts::writeParameter(os,pS1,toR,willBeActive);      
    wts::writeParameter(os,pS2,toR,willBeActive);      
    wts::writeParameter(os,pS3,toR,willBeActive);      
    wts::writeParameter(os,pS4,toR,willBeActive);      
    wts::writeParameter(os,pS5,toR,willBeActive);      
    wts::writeParameter(os,pS6,toR,willBeActive);      
    wts::writeParameter(os,pDevsS1,toR,willBeActive);      
    wts::writeParameter(os,pDevsS2,toR,willBeActive);      
    wts::writeParameter(os,pDevsS3,toR,willBeActive);      
    wts::writeParameter(os,pDevsS4,toR,willBeActive);      
    wts::writeParameter(os,pDevsS5,toR,willBeActive);      
    wts::writeParameter(os,pDevsS6,toR,willBeActive);      
    
    //fishery parameters
    wts::writeParameter(os,pHM,toR,willBeActive);      
    wts::writeParameter(os,pLnC,toR,willBeActive);      
    wts::writeParameter(os,pLnDCT,toR,willBeActive);      
    wts::writeParameter(os,pLnDCX,toR,willBeActive);      
    wts::writeParameter(os,pLnDCM,toR,willBeActive);      
    wts::writeParameter(os,pLnDCXM,toR,willBeActive);      
    wts::writeParameter(os,pDevsLnC,toR,willBeActive);      
    
    //survey parameters
    wts::writeParameter(os,pLnQ,toR,willBeActive);      
    wts::writeParameter(os,pLnDQT,toR,willBeActive);      
    wts::writeParameter(os,pLnDQX,toR,willBeActive);      
    wts::writeParameter(os,pLnDQM,toR,willBeActive);      
    wts::writeParameter(os,pLnDQXM,toR,willBeActive);      
    
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
        
    //write active parameters to rpt::echo
    rpt::echo<<"Finished phase "<<current_phase()<<endl;
    if (last_phase()) {
        //write report as R file
        ReportToR(report,1,rpt::echo);
        //write parameter values to csv
        ofstream os("TCSAM2015.final_params.active.csv", ios::trunc);
        writeParameters(os,0,1);
        os.close();
    }
    
}

void model_parameters::between_phases_calculations(void)
{
    rpt::echo<<endl<<endl<<"#---------------------"<<endl;
    rpt::echo<<"Starting phase "<<current_phase()<<" of "<<initial_params::max_number_phases<<endl;
}

void model_parameters::final_calcs()
{
    {cout<<"writing model sim data to file"<<endl;
        ofstream echo1; echo1.open("ModelSimData.dat", ios::trunc);
        writeSimData(echo1,0,cout,ptrSimMDS);
    }
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::app);
        mcmc<<"NULL)"<<endl;
        mcmc.close();
    }
    
    long hour,minute,second;
    double elapsed_time;
    
    time(&finish); 
    elapsed_time = difftime(finish,start);
    
    hour = long(elapsed_time)/3600;
    minute = long(elapsed_time)%3600/60;
    second = (long(elapsed_time)%3600)%60;
    cout << endl << endl << "Starting time: " << ctime(&start);
    cout << "Finishing time: " << ctime(&finish);
    cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
    
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000,5000,5000,5000,5000,5000,10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{0.5,0.1,.01,.001,1e-3,1e-4}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 1000000000; //must be smaller than 2,147,483,647
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(40000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(4000);
  time(&start);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
