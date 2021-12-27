#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "PSO_for_cur.h"

int main(int argc,char* argv[]) {

    char *id = (char*)"068";
    int pool = 5;
    bool fixsampling8=true;
    bool fixcurvature8=true;
    bool fixsampling20=false;
    bool fixcurvature20=false;
    bool fixTO=false;
    bool fixc=true;
    bool usepast=false;
    bool tofit=false;
    bool iftest=false;
 
    for(int i=1;i<argc;i++){
      if(!strcmp(argv[i],"-id")){
        id=argv[++i];
      }else if(!strcmp(argv[i],"-pool")){
        pool=stoi(argv[++i]);
      }else if(!strcmp(argv[i],"-fixsampling8")){
        fixsampling8 = true;
      }else if(!strcmp(argv[i],"-fixsampling8")){
        fixcurvature8 = true;
      }else if(!strcmp(argv[i],"-fixsampling20")){
        fixsampling20 = true;
      }else if(!strcmp(argv[i],"-fixsampling20")){
        fixcurvature20 = true;
      }else if(!strcmp(argv[i],"-fixTO")){
        fixTO = true;
      }else if(!strcmp(argv[i],"-fixc")){
        fixc = true;
      }else if(!strcmp(argv[i],"-past")){
        usepast = true;
      }else if(!strcmp(argv[i],"-fit")){
        tofit = true;
      }else if(!strcmp(argv[i],"-test")){
        iftest = true;
      }else{
      std::cout<<"Can't recognize "<<argv[i]<<std::endl;
      exit(1);
      }
    }


    raven cuteone(id,pool,&fixsampling8,&fixcurvature8,&fixsampling20,&fixcurvature20,&fixTO,&fixc,&usepast);

    if (iftest){
      while (strcmp(cuteone.feeling,"craptacular"))
      {
        cuteone.flytest();
        cuteone.move_r68();
      }
    }else{
      if (!tofit){
        if (usepast){
          cuteone.recon();
          cuteone.superfit();
        }else{
          while (strcmp(cuteone.feeling,"craptacular"))
          {
          // cuteone.flytest();
          cuteone.recon();
          // cuteone.sigmafit();

          // cuteone.skymap();
          // cuteone.psffit();
          cuteone.superfit();
          // cuteone.starfit();

          cuteone.move_psf();
          // cuteone.move_r68();
          }
        }
      }else{
        // cuteone.superfit();
        cuteone.starfit();
        // cuteone.search();
      }
    }
    
    cout << "I'm tired!" <<endl<<endl<<"REALLY!!"<<endl;
}