#include <iostream>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include <TFile.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>
#include <TTree.h>
#include <TH2D.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>

//using namespace ROOT::Math;
using namespace std;

//typedef ROOT::Math::SVector>float,4>    SVector4;

class raven 
{
    public:
        raven(char* _id,int _pool,bool* _fixsampling8,bool* _fixcurvature8,bool* _fixsampling20,bool* _fixcurvature20,bool* _fixTO,bool* _fixc,bool* _usepast)  //const float uprange[4], const float lowrange[4]
        {
            if(init(_id, _pool, _fixsampling8, _fixcurvature8, _fixsampling20, _fixcurvature20, _fixTO, _fixc, _usepast))
            {
                cerr<<"Failed to initiate the particle!"<<endl;
                exit(1);
            }
        }

        // TLorentzVector x;
        // TLorentzVector v;
        TLorentzVector x8;
        TLorentzVector v8;
        TLorentzVector x20;
        TLorentzVector v20;
        float vc=0;
        double z=0;
        int init(char* _id,int _pool,bool* _fixsampling8,bool* _fixcurvature8,bool* _fixsampling20,bool* _fixcurvature20,bool* _fixTO,bool* _fixc,bool* _usepast);

        float cur8[4] = {-1.5, -0.3, 0.01, 0.001};
        // float cur20[4] = {-0.550683, -0.450342, 0.0638406, 0.000681859};  // r_2.5
        float cur20[4] = {-1.5, -0.3, 0.01, 0.00085};
        // float cur20[4] = { -1.87393, 0.309068, 0.0638406, 0.000681859}; // sam_15
        // float cur20[4] = { -1.87393, 0.309068, 0.10288, 0.000369}; // psf_new
        // float cur20[4] = { -2.6841, -0.688197, 0.10288, 0.000369}; // str
        // float cur20[4] = {-5.511463, -0.488383, 0.030715, 0.000415}; // best

        double uprange8[4] = {1, 1, 0.16, 0.0013};
        double lowrange8[4] = {-15, -4, 0, 0};

        // double uprange8[4] = {-3, 1, 0.16, 0.0013}; // poo1_1
        // double lowrange8[4] = {-8, -1.5, 0, 0};

        // double uprange8[4] = {1, 1, 0.16, 0.0013};  // pool1_2
        // double lowrange8[4] = {-10, -1.5, 0, 0};

        // double uprange8[4] = {6, 5, 0.16, 0.0013};  // pool1_4
        // double lowrange8[4] = {-11, -5, 0, 0};

        // double uprange8[4] = {2, 2, 0.16, 0.0013};
        // double lowrange8[4] = {-4, -2, 0, 0};

        // double uprange8[4] = {12, 3, 0.16, 0.0013};
        // double lowrange8[4] = {-12, -3, 0, 0};   

        // double uprange8[4] = {5.12, 5.12, 5.12, 5.12};        //for flytest
        // double lowrange8[4] = {-5.12, -5.12, -5.12, -5.12};       //

        // double uprange20[4] = {5.12, 5.12, 5.12, 5.12};        //for flytest
        // double lowrange20[4] = {-5.12, -5.12, -5.12, -5.12};       //

        // double uprange20[4] = {1, 1, 0.2, 0.0013};
        // double lowrange20[4] = {-2.0, -1.0, 0, 0.0004};

        // double uprange20[4] = {1.0, 0.8, 0.2, 0.0015};          //str
        // double lowrange20[4] = {-3.0, -0.8, 0, 0};

        // double uprange20[4] = {1.0, 0.8, 0.16, 0.0013};         //
        // double lowrange20[4] = {-4.0, -1, 0, 0};

        // double uprange20[4] = {-3, 1, 0.16, 0.0013};
        // double lowrange20[4] = {-8, -1.5, -0.05, -0.0006};

        // double uprange20[4] = {-1, 0, 0.16, 0.0013};
        // double lowrange20[4] = {-15, -4, 0, 0};

        double uprange20[4] = {1, 1, 0.16, 0.0013};
        double lowrange20[4] = {-15, -4, 0, 0};

        float cur8c = 0;
        float cur8cu = 1;
        float cur8cl = -7;
        // float cur8cu = 5.12;
        // float cur8cl = -5.12;
        float bestcur8cm;
        float bestcur8cr;

        Double_t sigma=0.274936;

        TVector3 TO;
        TVector3 TOm;
        TVector3 TOr;
        TVector3 TOv;

        // float TOz[3] = {11.81, 14.71, 12.71};
        // float TOu[3] = {13.31, 16.21, 14.21};
        // float TOl[3] = {10.31, 13.21, 11.21};

        float TOz[3] = {12, 15, 13};
        float TOu[3] = {13.5, 16.5, 14.5};
        float TOl[3] = {10.5, 13.5, 11.5};
        // float TOu[3] = {5.12, 5.12, 5.12};
        // float TOl[3] = {-5.12, -5.12, -5.12};

        double uprange[4];
        double lowrange[4];

        double cur[4];

        int pool;
        char *id;
        bool fixsampling8;
        bool fixcurvature8;
        bool fixsampling20;
        bool fixcurvature20;
        bool usepast;
        bool fixTO;
        bool fixc;
        bool move_line=false;

        char *feeling = (char*)"very good";
        int heart=35;
        int myheart=heart;

        int recon();
        int skymap();
        int psffit();
        int sigmafit();
        int flytest();
        int move_psf();
        int move_r68();
        int search();
        int superfit();
        int starfit();
        int fixline();
        string Trim(string& str);
        float line(double x);
        TLorentzVector bestmemory8;
        TLorentzVector bestmemory20;
        TLorentzVector bestraven8;
        TLorentzVector bestraven20;
        TLorentzVector array_n;
        TLorentzVector array3_n;
        TLorentzVector array2_n;
        TLorentzVector array4_n;
        float omiga = 0.8;           //Inertia factor
        float c1 = 1;                //accelerate const
        float c2 = 1;
        float vscale = 0.3;
        float water = 1.1;
        float fire = 1.2;
        TRandom *rd = new TRandom(0);
        float bestfood = 100;
        float bestmemoryfood = 100;
        float bestsig = 0;
        float bestmsig = 0;
        float food_chain[19];
        float food_chain2[19];

        float rc_cutL = 0.5;
        float rc_cutH = 1.5;
        float h_cutL = 3000;
        float h_cutH = 5000;
        float sigma_cut = 0.1;
};

int raven::init(char* _id,int _pool,bool* _fixsampling8,bool* _fixcurvature8,bool* _fixsampling20,bool* _fixcurvature20,bool* _fixTO,bool* _fixc,bool* _usepast)  //const float uprange[4], const float lowrange[4]
{
    id = _id;
    pool = _pool;
    fixcurvature8 = *_fixcurvature8;
    fixcurvature20 = *_fixcurvature20;
    usepast = *_usepast;
    fixsampling8 = *_fixsampling8;
    fixsampling20 = *_fixsampling20;
    fixTO = *_fixTO;
    fixc = *_fixc;
    array_n.SetXYZT(0.,0.,-1.,1/125);
    array3_n.SetXYZT(1.,1.,1.,1/125);

    cout<<"I'm a little raven! my id is "<<id<<" i'm on the pool "<< pool << "; i will fixcur8? " <<fixcurvature8<<"; i will fixsamp8? "<<fixsampling8<<endl;
    cout<<"i will fixcur20? "<<fixcurvature20<<"; i will fixsamp20? "<<fixsampling20<<endl<<endl;

    double vmaxc = cur8cu-cur8cl;



    double vmax80 = uprange8[0]-lowrange8[0];
    double vmax81 = uprange8[1]-lowrange8[1];
    double vmax82 = uprange8[2]-lowrange8[2];
    double vmax83 = uprange8[3]-lowrange8[3];

    double vmax200 = uprange20[0]-lowrange20[0];
    double vmax201 = uprange20[1]-lowrange20[1];
    double vmax202 = uprange20[2]-lowrange20[2];
    double vmax203 = uprange20[3]-lowrange20[3];

    double TOmax[3] = {TOu[0]-TOl[0], TOu[1]-TOl[1], TOu[2]-TOl[2]};

    if (!fixsampling8 && !fixcurvature8){
        x8.SetXYZT(rd->Uniform(uprange8[0], lowrange8[0]), rd->Uniform(uprange8[1], lowrange8[1]), rd->Uniform(uprange8[2], lowrange8[2]), rd->Uniform(uprange8[3], lowrange8[3]));

        v8.SetXYZT(vscale*rd->Uniform(-vmax80, vmax80), vscale*rd->Uniform(-vmax81, vmax81), vscale*rd->Uniform(-vmax82, vmax82), vscale*rd->Uniform(-vmax83, vmax83));

    }else if(fixsampling8 && !fixcurvature8){
        x8.SetXYZT(cur8[0], cur8[1], rd->Uniform(uprange8[2], lowrange8[2]), rd->Uniform(uprange8[3], lowrange8[3]));

        v8.SetXYZT(0, 0, vscale*rd->Uniform(-vmax82, vmax82), vscale*rd->Uniform(-vmax83, vmax83));
    }else if(!fixsampling8 && fixcurvature8){
        x8.SetXYZT(rd->Uniform(uprange8[0], lowrange8[0]), rd->Uniform(uprange8[1], lowrange8[1]), cur8[2], cur8[3]);

        v8.SetXYZT(vscale*rd->Uniform(-vmax80, vmax80), vscale*rd->Uniform(-vmax81, vmax81), 0, 0);
    }else{
        x8.SetXYZT(cur8[0], cur8[1], cur8[2], cur8[3]);

        v8.SetXYZT(0, 0, 0, 0);
    }

    if (!fixsampling20 && !fixcurvature20){
        x20.SetXYZT(rd->Uniform(uprange20[0], lowrange20[0]), rd->Uniform(uprange20[1], lowrange20[1]), rd->Uniform(uprange20[2], lowrange20[2]), rd->Uniform(uprange20[3], lowrange20[3]));

        v20.SetXYZT(vscale*rd->Uniform(-vmax200, vmax200), vscale*rd->Uniform(-vmax201, vmax201), vscale*rd->Uniform(-vmax202, vmax202), vscale*rd->Uniform(-vmax203, vmax203));

    }else if(fixsampling20 && !fixcurvature20){
        x20.SetXYZT(cur20[0], cur20[1], rd->Uniform(uprange20[2], lowrange20[2]), rd->Uniform(uprange20[3], lowrange20[3]));

        v20.SetXYZT(0, 0, vscale*rd->Uniform(-vmax202, vmax202), vscale*rd->Uniform(-vmax203, vmax203));
    }else if(!fixsampling20 && fixcurvature20){
        x20.SetXYZT(rd->Uniform(uprange20[0], lowrange20[0]), rd->Uniform(uprange20[1], lowrange20[1]), cur20[2], cur20[3]);

        v20.SetXYZT(vscale*rd->Uniform(-vmax200, vmax200), vscale*rd->Uniform(-vmax201, vmax201), 0, 0);
    }else{
        x20.SetXYZT(cur20[0], cur20[1], cur20[2], cur20[3]);

        v20.SetXYZT(0, 0, 0, 0);
    }


    if (!fixTO){
        TO.SetXYZ(rd->Uniform(TOl[0],TOu[0]),rd->Uniform(TOl[1],TOu[1]),rd->Uniform(TOl[2],TOu[2]));
        TOv.SetXYZ(rd->Uniform(-TOmax[0],TOmax[0]), rd->Uniform(-TOmax[1],TOmax[1]), rd->Uniform(-TOmax[2],TOmax[2]));
    }else{
        TO.SetXYZ(TOz[0],TOz[1],TOz[2]);
        TOv.SetXYZ(0, 0, 0);
    }
    
    if (pool == 1){
        vc = 0;
    }else if(pool == 5){
        if (!fixc){
            vc = rd->Uniform(-vmaxc,vmaxc);
            cur8c = rd->Uniform(cur8cl, cur8cu);
        }else{
            vc = 0;
        }
    }else{
        vc = 0;
    }

    bestmemory20 = x20;
    bestraven20 = x20;
    bestmemory8 = x8;
    bestraven8 = x8;
    TOm = TO;
    TOr = TO;
    bestcur8cm = cur8c;
    bestcur8cr = cur8c;

    cout<<"The coordinate of raven "<< id << " is (" << x8.X() << "," << x8.Y() << "," << x8.Z() << "," << x8.T() << "," <<x20.X() << "," << x20.Y() << "," << x20.Z() << "," << x20.T()<<","<< cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ")" <<endl;
    cout<<"The velocity of raven "<< id << " is   (" << v8.X() << "," << v8.Y() << "," << v8.Z() << "," << v8.T() << "," <<v20.X() << "," << v20.Y() << "," << v20.Z() << "," << v20.T()<<","<< vc << "," << TOv[0] << "," << TOv[1] << "," << TOv[2] << ") "<<endl<<endl;
    return 0;
}

int raven::move_psf()
{
    if (!usepast){
    string line2="000, 0, 0, 0, 0, 100, 100";
    string lastline2="000, 0, 0, 0, 0, 100, 100";
    string number;
    ifstream fp2(Form("./magic/raven%s.csv", id),ios::in); //定义声明一个ifstream对象，指定文件路径
    //getline(fp,line); //跳过列名，第一行不做处理
    while (getline(fp2,line2)){ //循环读取每行数据
        //vector<float> food_chain{2,0,0,0,0,2};
        //将一行数据按'，'分割
        string::size_type idx = line2.find(",");
        string::size_type idy = line2.find(".");
        string::size_type idz = line2.find("#");
        if ( idx == string::npos || idy == string::npos || idz != string::npos) {
            continue;
        }
        lastline2 = line2;
    }

    istringstream readstr2(lastline2); //string数据流化
    for(int j = 0;j < 19;j++){ //可根据数据的实际情况取循环获取
        getline(readstr2,number,','); //循环读取数据
        if (number != ""){
            //cout<<number<<endl;
            food_chain[j]=stof(Trim(number)); //字符串传int .c_str()
        }
    }

    //compare with sky.csv
    ifstream fp("./sky.csv",ios::in); //定义声明一个ifstream对象，指定文件路径
    string line = "000, 0, 0, 0, 0, 100, 100";
    string lastline = "000, 0, 0, 0, 0, 100, 100";
    //getline(fp,line); //跳过列名，第一行不做处理
    while (getline(fp,line)){ //循环读取每行数据
        string::size_type idx = line.find(",");
        string::size_type idy = line.find(".");
        string::size_type idz = line.find("#");
        if ( idx == string::npos || idy == string::npos || idz != string::npos) {
            continue;
        }
        lastline = line;
    }

    istringstream readstr(lastline); //string数据流化
    for(int j = 0;j < 19;j++){ //可根据数据的实际情况取循环获取
    getline(readstr,number,','); //循环读取数据
        if (number != ""){
        //cout<<number<<endl;
        food_chain2[j]=stof(Trim(number)); //字符串传int .c_str()
        }
    }

    if( food_chain2[13]<bestfood) {      // || food_chain2[18]>bestsig
        cout << "That piece of food is the biggest! The lucky raven is at ("<< food_chain2[1] << "," << food_chain2[2]<< "," <<  food_chain2[3]<< "," <<  food_chain2[4] << "," <<food_chain2[5]<< "," << food_chain2[6] << "," << food_chain2[7] << "," << food_chain2[8] << "," << food_chain2[9] << "," << food_chain2[10] << "," << food_chain2[11] << "," << food_chain2[12] << ")" <<endl;
        cout << "The best food is : " << food_chain2[13] <<endl<<endl;

        bestfood = food_chain2[13];
        bestsig = food_chain2[18];
        bestcur8cr = food_chain2[9];
        TOr.SetXYZ(food_chain2[10], food_chain2[11], food_chain2[12]);
        bestraven8.SetXYZT(food_chain2[1], food_chain2[2], food_chain2[3], food_chain2[4]);
        bestraven20.SetXYZT(food_chain2[5], food_chain2[6], food_chain2[7], food_chain2[8]);
    }


    if (food_chain[13] < bestmemoryfood  && food_chain[13] >sigma_cut && food_chain[14]>rc_cutL && food_chain[14]<rc_cutH && food_chain[15] > h_cutL && food_chain[15]<h_cutH)   //|| food_chain[18] > bestmsig)
    {
        cout << "That piece of food is the biggest in my memory! It's at ("<< food_chain[1] << "," << food_chain[2]<< "," <<  food_chain[3]<< "," <<  food_chain[4] << "," <<food_chain[5]<< "," << food_chain2[6] << "," << food_chain2[7] << "," << food_chain2[8] << "," << food_chain2[9] << "," << food_chain2[10] << "," << food_chain2[11] << "," << food_chain2[12] << ")" <<endl;
        cout << "The best memory is : " << food_chain[13] <<endl<<endl;


        bestmemoryfood = food_chain[13];
        bestmsig = food_chain[18];
        bestmemory8.SetXYZT(food_chain[1], food_chain[2], food_chain[3], food_chain[4]);
        bestmemory20.SetXYZT(food_chain[5], food_chain[6], food_chain[7], food_chain[8]);
        bestcur8cm = food_chain[9];
        TOm.SetXYZ(food_chain[10], food_chain[11], food_chain[12]);
        vscale = vscale/water;
        if (vscale <= 0.1) vscale = 0.1;


        if (bestmemoryfood < food_chain2[13]){    // || bestmsig > food_chain2[18]
            ofstream outFile;
            outFile.open("./sky.csv", ios::app);
            outFile << id << ", " << food_chain[1] << ", "  << food_chain[2] << ", "  << food_chain[3] << ", "  << food_chain[4] << ", "  <<food_chain[5]<< ", " << food_chain[6] << ", "<< food_chain[7]<<", "<<food_chain[8]<<", "<<food_chain[9]<<", "<<food_chain[10]<<", "<<food_chain[11]<<", "<< food_chain[12] << ", " << food_chain[13] << ", " << food_chain[14] << ", " << food_chain[15] << ", " << food_chain[16] << ", " << food_chain[17] << ", " << food_chain[18] <<endl;
            outFile.close();
        }

        myheart += 1;
        if (myheart<0.75*heart&& myheart>0.5*heart){
            feeling = (char*)"nice";
        }else if (myheart<0.5*heart && myheart>0.25*heart){
            feeling = (char*)"OK";
        }else if (myheart>0&&myheart<0.25*heart){
            feeling = (char*)"bad";
        }else if (myheart == 0){
            feeling = (char*)"craptacular";
        }else if (myheart>0.75*heart){
            feeling = (char*)"very good";
        }
    }
    else{
        cout << "No bigger food found in my memory!"<<endl<<endl;
        vscale *= fire;
        if (vscale >= 2) vscale = 2;
        myheart -= 1;
        if (myheart<0.75*heart&& myheart>0.5*heart){
            feeling = (char*)"nice";
        }else if (myheart<0.5*heart && myheart>0.25*heart){
            feeling = (char*)"OK";
        }else if (myheart>0&&myheart<0.25*heart){
            feeling = (char*)"bad";
        }else if (myheart == 0){
            feeling = (char*)"craptacular";
        }else if (myheart>0.75*heart){
            feeling = (char*)"very good";
        }
        cout << Form("I'm feeling %s",feeling) << endl<<endl;
    }

    cout<<"The coordinate of raven "<< id << " is (" << x8.X() << "," << x8.Y() << "," << x8.Z() << "," << x8.T() << "," <<x20.X() << "," << x20.Y() << "," << x20.Z() << "," << x20.T()<<","<< cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ")" <<endl;
    cout<<"The velocity of raven "<< id << " is   (" << v8.X() << "," << v8.Y() << "," << v8.Z() << "," << v8.T() << "," <<v20.X() << "," << v20.Y() << "," << v20.Z() << "," << v20.T()<<","<< vc << "," << TOv[0] << "," << TOv[1] << "," << TOv[2] << ")"<<endl<<endl;

    if (move_line){
        array2_n = vscale*(c1*rd->Uniform(0,1)*(bestmemory8 - x8) + c2*rd->Uniform(0,1)*(bestraven8 - x8));
        array4_n = array2_n;
        array2_n = (array2_n*array_n)*array3_n;
        array2_n(0) = array4_n(0);
        array2_n(1) = array4_n(1);
        v8 =  omiga*v8 + array2_n;
        v20 =  omiga*v20 + array2_n;
    }else{
        v8 =  omiga*v8 + vscale*(c1*rd->Uniform(0,1)*(bestmemory8 - x8) + c2*rd->Uniform(0,1)*(bestraven8 - x8));
        v20 =  omiga*v20 + vscale*(c1*rd->Uniform(0,1)*(bestmemory20 - x20) + c2*rd->Uniform(0,1)*(bestraven20 - x20));
        vc = omiga*vc + vscale*(c1*rd->Uniform(0,1)*(bestcur8cm - cur8c) + c2*rd->Uniform(0,1)*(bestcur8cr - cur8c));
        TOv = omiga*TOv + vscale*(c1*rd->Uniform(0,1)*(TOm - TO) + c2*rd->Uniform(0,1)*(TOr - TO));
    }

    x8 =  x8 + v8;
    x20 = x20 + v20;
    cur8c = cur8c + vc;
    TO = TO + TOv;
    for (int i=0; i<4;i++){
        while(x8[i]>uprange8[i] || x8[i]<lowrange8[i]){
            if(x8[i]>uprange8[i]){
                cout<<"boom"<<endl;
                x8[i]=uprange8[i]-(x8[i]-uprange8[i]);
                v8[i]*=-1;
            }else if(x8[i]<lowrange8[i]){
                cout<<"boom"<<endl;
                x8[i]=lowrange8[i]+(lowrange8[i]-x8[i]);
                v8[i]*=-1;
            }
        }
    }

    for (int i=0; i<4;i++){
        while(x20[i]>uprange20[i] || x20[i]<lowrange20[i]){
            if(x20[i]>uprange20[i]){
                cout<<"boom"<<endl;
                x20[i]=uprange20[i]-(x20[i]-uprange20[i]);
                v20[i]*=-1;
            }else if(x20[i]<lowrange20[i]){
                cout<<"boom"<<endl;
                x20[i]=lowrange20[i]+(lowrange20[i]-x20[i]);
                v20[i]*=-1;
            }
        }
    }

    for (int i=0; i<3;i++){
        while(TO[i]>TOu[i] || TO[i]<TOl[i]){
            if(TO[i]>TOu[i]){
                cout<<"boom"<<endl;
                TO[i]=TOu[i]-(TO[i]-TOu[i]);
                TOv[i]*=-1;
            }else if(TO[i]<TOl[i]){
                cout<<"boom"<<endl;
                TO[i]=TOl[i]+(TOl[i]-TO[i]);
                TOv[i]*=-1;
            }
        }
    }

    while(cur8c>cur8cu || cur8c<cur8cl){
        if(cur8c>cur8cu){
            cout<<"boomc"<<endl;
            cur8c=cur8cu-(cur8c-cur8cu);
            vc*=-1;
        }else if(cur8c<cur8cl){
            cout<<"boomc"<<endl;
            cur8c=cur8cl+(cur8cl-cur8c);
            vc*=-1;
        }
    }

    cout<<"i fly;"<<endl;
    cout<<"The coordinate of raven "<< id << " is (" << x8.X() << "," << x8.Y() << "," << x8.Z() << "," << x8.T() << "," <<x20.X() << "," << x20.Y() << "," << x20.Z() << "," << x20.T()<<","<< cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ")" <<endl;
    cout<<"The velocity of raven "<< id << " is   (" << v8.X() << "," << v8.Y() << "," << v8.Z() << "," << v8.T() << "," <<v20.X() << "," << v20.Y() << "," << v20.Z() << "," << v20.T()<<","<< vc << "," << TOv[0] << "," << TOv[1] << "," << TOv[2] << ")"<<endl<<endl;
    }
    return 0;
}

int raven::move_r68()
{
    string line2="000, 0, 0, 0, 0, 100";
    string lastline2="000, 0, 0, 0, 0, 100";
    string number;
    ifstream fp2(Form("./magic/raven%s.csv", id),ios::in); //定义声明一个ifstream对象，指定文件路径
    //getline(fp,line); //跳过列名，第一行不做处理
    while (getline(fp2,line2)){ //循环读取每行数据
        //vector<float> food_chain{2,0,0,0,0,2};
        //将一行数据按'，'分割
        string::size_type idx = line2.find(",");
        string::size_type idy = line2.find(".");
        string::size_type idz = line2.find("#");
        if ( idx == string::npos || idy == string::npos || idz != string::npos) {
            continue;
        }
        lastline2 = line2;
    }

    istringstream readstr2(lastline2); //string数据流化
    for(int j = 0;j < 14;j++){ //可根据数据的实际情况取循环获取
        getline(readstr2,number,','); //循环读取数据
        if (number != ""){
            //cout<<number<<endl;
            food_chain[j]=stof(Trim(number)); //字符串传int .c_str()
        }
    }

    //compare with sky.csv
    ifstream fp("./sky.csv",ios::in); //定义声明一个ifstream对象，指定文件路径
    string line = "000, 0, 0, 0, 0, 1000";
    string lastline = "000, 0, 0, 0, 0, 1000";
    //getline(fp,line); //跳过列名，第一行不做处理
    while (getline(fp,line)){ //循环读取每行数据
        string::size_type idx = line.find(",");
        string::size_type idy = line.find(".");
        string::size_type idz = line.find("#");
        if ( idx == string::npos || idy == string::npos || idz != string::npos) {
            continue;
        }
        lastline = line;
    }

    istringstream readstr(lastline); //string数据流化
    for(int j = 0;j < 14;j++){ //可根据数据的实际情况取循环获取
    getline(readstr,number,','); //循环读取数据
        if (number != ""){
        //cout<<number<<endl;
        food_chain2[j]=stof(Trim(number)); //字符串传int .c_str()
        }
    }

    if( food_chain2[13]<bestfood ) {
        cout << "That piece of food is the biggest! The lucky raven is at ("<< food_chain2[1] << "," << food_chain2[2]<< "," <<  food_chain2[3]<< "," <<  food_chain2[4] << "," <<food_chain2[5]<< "," << food_chain2[6] << "," << food_chain2[7] << "," << food_chain2[8] << "," << food_chain2[9] << "," << food_chain2[10] << "," << food_chain2[11] << "," << food_chain2[12] << ")" <<endl;
        cout << "The best food is : " << food_chain2[13] <<endl<<endl;

        bestfood = food_chain2[13];
        bestcur8cr = food_chain2[9];
        TOr.SetXYZ(food_chain2[10], food_chain2[11], food_chain2[12]);
        bestraven8.SetXYZT(food_chain2[1], food_chain2[2], food_chain2[3], food_chain2[4]);
        bestraven20.SetXYZT(food_chain2[5], food_chain2[6], food_chain2[7], food_chain2[8]);
    }


    if (food_chain[13] < bestmemoryfood)
    {
        cout << "That piece of food is the biggest in my memory! It's at ("<< food_chain[1] << "," << food_chain[2]<< "," <<  food_chain[3]<< "," <<  food_chain[4] << "," <<food_chain[5]<< "," << food_chain2[6] << "," << food_chain2[7] << "," << food_chain2[8] << "," << food_chain2[9] << "," << food_chain2[10] << "," << food_chain2[11] << "," << food_chain2[12] << ")" <<endl;
        cout << "The best memory is : " << food_chain[13] <<endl<<endl;


        bestmemoryfood = food_chain[13];
        bestmemory8.SetXYZT(food_chain[1], food_chain[2], food_chain[3], food_chain[4]);
        bestmemory20.SetXYZT(food_chain[5], food_chain[6], food_chain[7], food_chain[8]);
        bestcur8cm = food_chain[9];
        TOm.SetXYZ(food_chain[10], food_chain[11], food_chain[12]);
        vscale = vscale/water;
        if (vscale <= 0.1) vscale = 0.1;

        if (bestmemoryfood < food_chain2[13] ){
            ofstream outFile;
            outFile.open("./sky.csv", ios::app);
            outFile << id << ", " << food_chain[1] << ", "  << food_chain[2] << ", "  << food_chain[3] << ", "  << food_chain[4] << ", "  << food_chain[5] << ", " << food_chain[6] << "," << food_chain[7] << "," << food_chain[8] << "," << food_chain[9] << "," << food_chain[10] << "," << food_chain[11] << "," << food_chain[12] << "," << food_chain[13] << endl;
            outFile.close();
        }

        myheart += 1;
        if (myheart<0.75*heart&& myheart>0.5*heart){
            feeling = (char*)"nice";
        }else if (myheart<0.5*heart && myheart>0.25*heart){
            feeling = (char*)"OK";
        }else if (myheart>0&&myheart<0.25*heart){
            feeling = (char*)"bad";
        }else if (myheart == 0){
            feeling = (char*)"craptacular";
        }else if (myheart>0.75*heart){
            feeling = (char*)"very good";
        }
    }
    else{
        cout << "No bigger food found in my memory!"<<endl<<endl;
        vscale *= fire;
        if (vscale >= 2) vscale = 2;
        myheart -= 1;
        if (myheart<0.75*heart&& myheart>0.5*heart){
            feeling = (char*)"nice";
        }else if (myheart<0.5*heart && myheart>0.25*heart){
            feeling = (char*)"OK";
        }else if (myheart>0&&myheart<0.25*heart){
            feeling = (char*)"bad";
        }else if (myheart == 0){
            feeling = (char*)"craptacular";
        }else if (myheart>0.75*heart){
            feeling = (char*)"very good";
        }
        cout << Form("I'm feeling %s",feeling) << endl<<endl;
    }

    cout<<"The coordinate of raven "<< id << " is (" << x8.X() << "," << x8.Y() << "," << x8.Z() << "," << x8.T() << "," <<x20.X() << "," << x20.Y() << "," << x20.Z() << "," << x20.T()<<","<< cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ")" <<endl;
    cout<<"The velocity of raven "<< id << " is   (" << v8.X() << "," << v8.Y() << "," << v8.Z() << "," << v8.T() << "," <<v20.X() << "," << v20.Y() << "," << v20.Z() << "," << v20.T()<<","<< vc << "," << TOv[0] << "," << TOv[1] << "," << TOv[2] << ") "<<endl<<endl;    
    
    if (move_line){
        array2_n = vscale*(c1*rd->Uniform(0,1)*(bestmemory8 - x8) + c2*rd->Uniform(0,1)*(bestraven8 - x8));
        array4_n = array2_n;
        array2_n = (array2_n*array_n)*array3_n;
        array2_n(0) = array4_n(0);
        array2_n(1) = array4_n(1);
        v8 =  omiga*v8 + array2_n;
        v20 =  omiga*v20 + array2_n;
    }else{
        v8 =  omiga*v8 + vscale*(c1*rd->Uniform(0,1)*(bestmemory8 - x8) + c2*rd->Uniform(0,1)*(bestraven8 - x8));
        v20 =  omiga*v20 + vscale*(c1*rd->Uniform(0,1)*(bestmemory20 - x20) + c2*rd->Uniform(0,1)*(bestraven20 - x20));
        vc = omiga*vc + vscale*(c1*rd->Uniform(0,1)*(bestcur8cm - cur8c) + c2*rd->Uniform(0,1)*(bestcur8cr - cur8c));
        TOv = omiga*TOv + vscale*(c1*rd->Uniform(0,1)*(TOm - TO) + c2*rd->Uniform(0,1)*(TOr - TO));
    }
    x8 =  x8 + v8;
    x20 = x20 + v20;
    cur8c = cur8c + vc;
    TO = TO + TOv;
    for (int i=0; i<4;i++){
        while(x8[i]>uprange8[i] || x8[i]<lowrange8[i]){
            if(x8[i]>uprange8[i]){
                cout<<"boom"<<endl;
                x8[i]=uprange8[i]-(x8[i]-uprange8[i]);
                v8[i]*=-1;
            }else if(x8[i]<lowrange8[i]){
                cout<<"boom"<<endl;
                x8[i]=lowrange8[i]+(lowrange8[i]-x8[i]);
                v8[i]*=-1;
            }
        }
    }

    for (int i=0; i<4;i++){
        while(x20[i]>uprange20[i] || x20[i]<lowrange20[i]){
            if(x20[i]>uprange20[i]){
                cout<<"boom"<<endl;
                x20[i]=uprange20[i]-(x20[i]-uprange20[i]);
                v20[i]*=-1;
            }else if(x20[i]<lowrange20[i]){
                cout<<"boom"<<endl;
                x20[i]=lowrange20[i]+(lowrange20[i]-x20[i]);
                v20[i]*=-1;
            }
        }
    }

    for (int i=0; i<3;i++){
        while(TO[i]>TOu[i] || TO[i]<TOl[i]){
            if(TO[i]>TOu[i]){
                cout<<"boom"<<endl;
                TO[i]=TOu[i]-(TO[i]-TOu[i]);
                TOv[i]*=-1;
            }else if(TO[i]<TOl[i]){
                cout<<"boom"<<endl;
                TO[i]=TOl[i]+(TOl[i]-TO[i]);
                TOv[i]*=-1;
            }
        }
    }

    while(cur8c>cur8cu || cur8c<cur8cl){
        if(cur8c>cur8cu){
            cout<<"boomc"<<endl;
            cur8c=cur8cu-(cur8c-cur8cu);
            vc*=-1;
        }else if(cur8c<cur8cl){
            cout<<"boomc"<<endl;
            cur8c=cur8cl+(cur8cl-cur8c);
            vc*=-1;
        }
    }

    cout<<"i fly;"<<endl;
    cout<<"The coordinate of raven "<< id << " is (" << x8.X() << "," << x8.Y() << "," << x8.Z() << "," << x8.T() << "," <<x20.X() << "," << x20.Y() << "," << x20.Z() << "," << x20.T()<<","<< cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ")" <<endl;
    cout<<"The velocity of raven "<< id << " is   (" << v8.X() << "," << v8.Y() << "," << v8.Z() << "," << v8.T() << "," <<v20.X() << "," << v20.Y() << "," << v20.Z() << "," << v20.T()<<","<< vc << "," << TOv[0] << "," << TOv[1] << "," << TOv[2] << ")"<<endl<<endl;
    return 0;
}

int raven::recon()
{   
    const char *rprogramm = "/afs/ihep.ac.cn/users/c/caowy/home/recon_for_psf";
    const char *eosprefix="root://eos01.ihep.ac.cn/";
    const char *outdir=Form("/eos/user/c/caowy/wcda/recon_data/cur/%d", pool);
    // const char *infile="/eos/user/z/zwang/wcda/data/pool123_rec/v20_source_filter_loop_cure_cut/tot_0518_0820_r1_theta55_reduced.root";
    // const char *infile="/eos/user/z/zwang/wcda/data/pool123_rec/v20_source_filter_pincCutCao/tot_0518_1110_crab_r7_pincCutCao_reduced.root";
    const char *infile="/eos/user/z/zwang/wcda/data/pool123_rec/v20_source_filter_pincCutCao/tot_0701_1110_crab_r7_pincCutCao_reduced.root";
    const char *fname="tot_0518_0820_r1_theta55_reduced";
    char *rec;

    if ( !usepast ){
        if (pool == 1){
            rec = Form("%s/recon -i %s/%s -o %s/%s/%s_pool%d_rec%s.root -s %d -c %s/conf/configure_1021_20201110_gaobo_st_hy_charge_pool123_2021_h.txt -8inch %f %f %f %f %f %f -20inch %f %f %f %f %f %f -TO %f %f %f", rprogramm, eosprefix, infile, eosprefix, outdir, fname, pool, id, pool, rprogramm, 0., x8[0], x8[1], x8[2], x8[3], 0., 0., cur20[0], cur20[1], cur20[2], cur20[3], 0., TO[0], TO[1], TO[2]);
        }else if(pool ==5){
            rec = Form("%s/recon -i %s/%s -o %s/%s/%s_pool%d_rec%s.root -s %d -c %s/conf/configure_1021_20201110_gaobo_st_hy_charge_pool123_2021_h.txt -8inch %f %f %f %f %f %f -20inch %f %f %f %f %f %f -TO %f %f %f", rprogramm, eosprefix, infile, eosprefix, outdir, fname, pool, id, pool, rprogramm, cur8c, x20[0], x20[1], x20[2], x20[3], 0., 0., x20[0], x20[1], x20[2], x20[3], 0., TO[0], TO[1], TO[2]); //cur8[0], cur8[1], cur8[2], cur8[3],
        }else{
            rec = Form("%s/recon -i %s/%s -o %s/%s/%s_pool%d_rec%s.root -s %d -c %s/conf/configure_1021_20201110_gaobo_st_hy_charge_pool123_2021_h.txt -8inch %f %f %f %f %f %f -20inch %f %f %f %f %f %f -TO %f %f %f", rprogramm, eosprefix, infile, eosprefix, outdir, fname, pool, id, pool, rprogramm, 0., cur8[0], cur8[1], cur8[2], cur8[3], 0., 0., x20[0], x20[1], x20[2], x20[3], 0., TO[0], TO[1], TO[2]); //cur8[0], cur8[1], cur8[2], cur8[3],
        }
    }else{
        rec = Form("%s/recon -i %s/%s -o %s/%s/%s_pool%d_rec_past.root -s %d -c %s/conf/configure_1021_20201110_gaobo_st_hy_charge_pool123_2021_h.txt -8inch %f %f %f %f %f %f -20inch %f %f %f %f %f %f -TO %f %f %f", rprogramm, eosprefix, infile, eosprefix, outdir, fname, pool, pool, rprogramm, 0., -1.5, -0.3, 0.01, 0.001, 0., 0., -1.5, -0.3, 0.01, 0.00085, 0., 16., 17., 17.);
    }

    FILE *rr = popen(rec, "r"); // build pipe
	if (!rr)
		return 1;

    char tmp[1024];
	while (fgets(tmp, sizeof(tmp), rr) != NULL)
		std::cout << tmp << std::endl; // can join each line as string
	pclose(rr);

	return 0;
}

int raven::skymap()
{
    const char *dir="/home/lhaaso/caowy/sky_map";
    const char *datadir=Form("/eos/user/c/caowy/wcda/recon_data/cur/%d",pool);
    const char *outdirs="/home/lhaaso/caowy/forcurvature_V10/magic";
    const char *name=Form("Crab%s_%f", id, x20[0]);

    system(Form("file=`eos ls %s/*%s.root`; echo /eos/user/c/caowy/wcda/recon_data/cur/%d/${file} > %s/runshell/%s.txt", datadir, id, pool, dir, name));

    char *sky = Form("%s/src/sky_for_wz_cur %s/runshell/%s.txt %s/%s.root > ./magic/raven%s.log", dir, dir, name, outdirs, name, id);

    FILE *ss = popen(sky, "r"); // build pipe
	if (!ss)
		return 1;

    char tmp[1024];
	while (fgets(tmp, sizeof(tmp), ss) != NULL)
		std::cout << tmp << std::endl; // can join each line as string
	pclose(ss);

    return 0;
}

int raven::psffit()
{
    //const char *dirf="/eos/user/c/caowy/caowy/for_wz/cur";

    char *fitting = Form("./psfsrc/fly -id %s -cur %f %f %f %f", id, x8.X(), x8.Y(), x8.Z(), x8.T());

    FILE *ff = popen(fitting, "r"); // build pipe
	if (!ff)
		return 1;

    char tmp[1024];
	while (fgets(tmp, sizeof(tmp), ff) != NULL)
		std::cout << tmp << std::endl; // can join each line as string
	pclose(ff);


    return 0;
}

int raven::flytest()
{   
    ofstream outFile;
    outFile.open(Form("./magic/raven%s.csv", id), ios::app);

    z = pow(x8[2],2) - cos(2*M_PI*x8[2]) + pow(x8[3],2) - cos(2*M_PI*x8[3]) + pow(x8[0],2) - cos(2*M_PI*x8[0]) + pow(x8[1],2) - cos(2*M_PI*x8[1]) + pow(cur8c,2) - cos(2*M_PI*cur8c) + pow(TO[0],2) - cos(2*M_PI*TO[0]) + pow(TO[1],2) - cos(2*M_PI*TO[1]) + pow(TO[2],2) - cos(2*M_PI*TO[2]) + pow(x20[0],2) - cos(2*M_PI*x20[0]) + pow(x20[1],2) - cos(2*M_PI*x20[1]) + pow(x20[2],2) - cos(2*M_PI*x20[2]) + pow(x20[3],2) - cos(2*M_PI*x20[3]);
    
    outFile << id << ", " << x8[0] << ", "  << x8[1] << ", "  << x8[2] << ", "  << x8[3] << ", " << x20[0] << "," << x20[1] << "," << x20[2] << "," << x20[3] << "," << cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << ", " << z <<endl;
    outFile.close();

    return 0;
}


int raven::search()
{
    ifstream fp("./sky.csv",ios::in); //定义声明一个ifstream对象，指定文件路径
    string line;
    //getline(fp,line); //跳过列名，第一行不做处理
    while (getline(fp,line)){ //循环读取每行数据
        //vector<float> food_chain{2,0,0,0,0,2};
        string number;
        istringstream readstr(line); //string数据流化
        //将一行数据按'，'分割
        string::size_type idx = line.find(",");
        string::size_type idy = line.find(".");
        string::size_type idz = line.find("#");
        if ( idx == string::npos || idy == string::npos || idz != string::npos) {
            continue;
        }

        for(int j = 0;j < 8;j++){ //可根据数据的实际情况取循环获取
            getline(readstr,number,','); //循环读取数据
            //cout<<number<<endl;
            food_chain[j]=stof(Trim(number)); //字符串传int .c_str()
        }
        if (food_chain[5] < bestfood)
        {
            ofstream outFile;
            outFile.open("final_result.csv", ios::out);
            outFile << "That piece of food is the biggest! The lucky raven is at ("<< food_chain[1] << "," << food_chain[2]<< "," <<  food_chain[3]<< "," <<  food_chain[4] << ")" <<endl;
            outFile << "The best food is : " << food_chain[6] <<endl<<endl;
            outFile.close();
        }
    }
    return 0;
}

string raven::Trim(string& str)
{
	//str.find_first_not_of(" \t\r\n"),在字符串str中从索引0开始，返回首次不匹配"\t\r\n"的位置
	str.erase(0,str.find_first_not_of(" \t\r\n^@"));
	str.erase(str.find_last_not_of(" \t\r\n^@") + 1);
	return str;
}

int raven::sigmafit()
{
    string prefix = "root://eos01.ihep.ac.cn";
    TFile *fin = TFile::Open(Form("%s//eos/user/c/caowy/wcda/recon_data/cur/%d/tot_0518_0820_r1_theta55_reduced_pool%d_rec%s.root",prefix.data(), pool, pool, id));
    // TFile *fin = TFile::Open(Form("%s//eos/user/c/caowy/wcda/recon_data/cur/%d/tot_0518_0820_r1_theta55_reduced_pool%d_rec%s.root",prefix.data(), pool, pool, id));

    TTree *trec=0;
    TObject *tob=(TObject*)fin->Get("rec");
    trec=(TTree*)tob;

    double dec_src = 22.02;
    double ra_src  = 83.63;   
    Float_t ra; Float_t dec;

    int num = trec->GetEntries();
    //Double_t x_max=trec->GetMaximum("ra")-ra_src;
	//Double_t y_max=trec->GetMaximum("dec")-dec_src;
	//Double_t x_min=trec->GetMinimum("ra")-ra_src;
	//Double_t y_min=trec->GetMinimum("dec")-dec_src;
    trec->SetBranchAddress("ra",&ra);
    trec->SetBranchAddress("dec",&dec);
    TH1D *map=new TH1D("cur","for psf fiting",100000,0,1);

    float sigma;

    int max = 0;

    for (int ievent=0; ievent<num; ievent++){
        trec->GetEntry(ievent);
        map->Fill(abs(ra-ra_src));
        map->Fill(abs(dec-dec_src));
        if (abs(ra-ra_src)<1){
            max += 1;
        }
        if (abs(dec-dec_src)<1){
            max += 1;
        }
    }

    double all=0;

    cout<<"[";
    for (int i=0; i<100000; i++)
    {
        if (i%100 == 0){
        cout<<"#";
        }
        all += map->GetBinContent(i);
        sigma = i*0.00001;
        if (all/num/2 > 0.1){
            break;
        }
    }
    cout<<"]"<<endl<<endl;
    cout<<all<<" : "<<num<<"   sigma = "<< sigma <<endl<<endl;

    ofstream outFile;
    outFile.open(Form("./magic/raven%s.csv",id), ios::app);

    outFile<< id << "," << x20[0] << ","  << x20[1] << ","  << x20[2] << ","  << x20[3] << ","  << sigma <<endl;

    outFile.close();
    return 0;
}

float raven::line(double x)
{
    double y = 0.0012 - 0.08*x;
    return y;
}

int raven::fixline()
{   
    double sx = rd->Uniform(uprange[0], lowrange[0]);
    double cx = rd->Uniform(uprange[2], lowrange[2]);
    x20.SetXYZT(sx, line(sx), cx, line(cx));
    return 0;
}

int raven::superfit()
{
    int numofP = 8;
    int numofO = 50;
    double range = 2;

    double dec_src = 22.02;
    double ra_src  = 83.63;
    // double dec_src = 21.72;
    // double ra_src  = 83.53; 
    double r_tr = 1;
    Double_t on = 0;
    Double_t bg = 0;
    Float_t ra,dec;
    UInt_t nhit;
    Double_t mjd;
    Float_t zenc;
    Float_t xc;
    Float_t yc;
    UInt_t dnsec;
    string prefix = "root://eos01.ihep.ac.cn";

    TH2D *map=new TH2D("map","sky map",500,-range,range,500,-range,range);
    TH2D *map2=new TH2D("sky","sky map",500,-5,5,500,-5,5);
    TH1D *ring = new TH1D("Parameters","for psf fitting",numofO,0,range); //h2->GetXaxis()->GetXmax()-0.1);
    TF1 *func = new TF1("func","[7]*([2]/(2*TMath::Pi()*[1]*[1])*TMath::Gaus(x,[0],[1])/([2]+[3])+[3]/(2*TMath::Pi()*[5]*[5])*TMath::Gaus(x,[4],[5])/([2]+[3]))+[6]",0,range); //-[4]*x


    TFile *fin=0;
    TTree *trec=0;

    // if (!usepast){
    //     fin = TFile::Open(Form("%s//eos/user/c/caowy/wcda/recon_data/cur/%d/tot_0518_0820_r1_theta55_reduced_pool%d_rec%s.root",prefix.data(), pool, pool, id));
    // }else{
    //     fin = TFile::Open(Form("%s//eos/user/c/caowy/wcda/recon_data/cur/%d/tot_0518_0820_r1_theta55_reduced_pool%d_rec_past.root",prefix.data(), pool, pool));
    // }
    // fin = TFile::Open("/eos/user/c/caowy/wcda/recon_data/cur/5/tot_0518_1110_crab_r7_pincCutCao_reduced_pool5_rec_past.root");
    // fin = TFile::Open("/eos/user/c/caowy/wcda/recon_data/cur/5/tot_0518_1110_crab_r7_pincCutCao_reduced_pool5_rec_best.root");
    fin = TFile::Open("/eos/user/c/caowy/wcda/recon_data/cur/5/tot_0518_1110_crab_r7_pincCutCao_reduced_pool5_rec_bestc5.root");
    // fin = TFile::Open("/eos/user/z/zwang/wcda/data/pool123_rec/v20_source_filter_pincCutCao/rec/tot_0518_1110_crab_r7_pincCutCao_reduced_pool5_rec.root");
    // fin = TFile::Open("/eos/user/z/zwang/wcda/data/pool123_rec/v20_source_filter_pincCutCao/rec/tot_0518_1110_crab_r7_pincCutCao_reduced_pool5_rec_v21.root");

    TObject *tob=(TObject*)fin->Get("rec");
    trec=(TTree*)tob;
    trec->SetBranchAddress("ra",&ra);
    trec->SetBranchAddress("dec",&dec);
    trec->SetBranchAddress("nhit",&nhit);
    trec->SetBranchAddress("mjd",&mjd);
    trec->SetBranchAddress("theta",&zenc);
    trec->SetBranchAddress("xc",&xc);
    trec->SetBranchAddress("yc",&yc);
    trec->SetBranchAddress("dnsec",&dnsec);

    int Nrec = trec->GetEntries();
    for (int ievent=0; ievent<Nrec; ievent++){
        trec->GetEntry(ievent);
        if (dnsec<1000) continue;
        // if (nhit>200) continue;
        // if (xc>130 || xc<20 || yc>0 || yc<-110) continue;
        double xx = (ra - ra_src)*cos(dec*M_PI/180);
        double yy = dec - dec_src;
        Double_t r = sqrt(pow(xx,2)+pow(yy,2));

        map2->Fill(xx,yy);
        if (r>range) continue;
        // int n = ceil(r/(1.5/numofO));
        // double rr = n*(1.5/numofO);
        // double rl = (n-1)*(1.5/numofO);
        ring->Fill(r,1/(2*M_PI*r));
        // ring->Fill(r,1/(M_PI*pow(rr,2)-M_PI*pow(rl,2)));
        map->Fill(xx,yy);
    }

    TCanvas *c = new TCanvas("c1","PSF",500,600);
    c->Divide(1,2);
    c->cd(1);
    gPad -> SetGrid(1,1);

    func->SetParameters(0,0.3,ring->GetBinContent(1),ring->GetBinContent(1),0,0.1,ring->GetBinContent(1)/4);
    func->FixParameter(0,0);
    func->SetParName(0,"Mean1");
    func->FixParameter(4,0);
    func->SetParName(4,"Mean2");
    func->SetParLimits(1,0.1,range/2);
    func->SetParName(1,"Sigma1");
    func->SetParLimits(5,0.05,range/2);
    func->SetParName(5,"Sigma2");
    func->SetParLimits(2,0,5000);
    func->SetParName(2,"A");
    func->SetParLimits(3,0,5000);
    func->SetParName(3,"B");
    func->SetParLimits(6,0,2000);
    func->SetParName(6,"Const");
    func->SetParLimits(7,0,10000);
    func->SetParName(7,"h_gaus");
     //func->SetParLimits(4,0,10000);
    ring->Fit(func);

    gStyle -> SetOptFit(1);

    Double_t chi2 = func->GetChisquare();
    Double_t reduced_chi2 = chi2/(numofO-numofP);
    Double_t sigma1 = func->GetParameter(1);
    Double_t sigma2 = func->GetParameter(5);
    Double_t cstt = func->GetParameter(6);
    //Double_t gra = func->GetParameter(4);
    Double_t A = func->GetParameter(2);
    Double_t B = func->GetParameter(3);
    Double_t h = func->GetParameter(7);
    sigma = (A+B)/(A/sigma1+B/sigma2);
    r_tr = 1.51*sigma;//sigma; //1.5*
    std::cout << "chi2 = " << chi2 << std::endl;
    std::cout << "reduced_chi2 = " << reduced_chi2 << std::endl;
    std::cout << "sigma = "<< sigma << std::endl;   
    //TF1 *cst = new TF1("const",Form("%f",cstt),0,h2->GetXaxis()->GetXmax()-0.1);
    TF1 *line = new TF1("line","[0]-[1]*x",0,range); //-[1]*x
    TF1 *Dgaus = new TF1("Dgaus","[6]*([2]/(2*TMath::Pi()*[1]*[1])*TMath::Gaus(x,[0],[1])/([2]+[3])+[3]/(2*TMath::Pi()*[5]*[5])*TMath::Gaus(x,[4],[5])/([2]+[3]))",0,range);
    line->SetParameters(cstt,0); //, gra);
    Dgaus->SetParameters(0,sigma1,A,B,0,sigma2,h);

    Double_t hgaus = Dgaus->Eval(0);

    for (int ievent=0; ievent<Nrec; ievent++){
        trec->GetEntry(ievent);
        if (dnsec<1000) continue;
        // if (nhit<200) continue;
        // if (xc>130 || xc<20 || yc>0 || yc<-110) continue;
        double xx = (ra - ra_src)*cos(dec*M_PI/180);
        double yy = dec - dec_src;
        Double_t r = sqrt(pow(xx,2)+pow(yy,2));

        if (r<=r_tr) {
            on += 1;
        }
    }

    bg = cstt*M_PI*pow(r_tr,2)*numofO/range;

    Double_t alpha = M_PI*pow(r_tr,2)/pow(range,2);
    Double_t off = bg/alpha;
    Double_t sig = sqrt(2)*sqrt(on*log((1+alpha)/alpha*on/(on+off)) + off*log((1+alpha)*off/(on+off)));

    std::cout << "R_tr = "<< r_tr << "  On = " << on << "  Bg = " << bg << "  Sig = " << sig << endl;

    Dgaus->SetLineColor(43);
    line->SetLineColor(42);
    func->SetLineColor(45);
    ring->Draw("E1");

    ring->SetTitle(Form("Psf fiting of raven%s",id)); //id
    // ring->SetTitle(Form("Psf fiting of raven%s","129"));
    ring->GetXaxis()->SetTitle("#psi^{2}[deg^{2}]");
    ring->GetYaxis()->SetTitle("dN/d#Omega[ev/sr]");
    ring->GetXaxis() -> CenterTitle();
    ring->GetYaxis() -> CenterTitle();

    line->Draw("same");
    Dgaus->Draw("same");
    TPaveText *pt = new TPaveText(0.3,0.6,0.6,0.9,"brNDC");
    pt->AddText(Form("Sigma = %f", sigma));
    pt->AddText("sam8_1  :    sam8_2   :   cur8_1    :  cur8_2");
    pt->AddText(Form("(%f, %f, %f, %f)",x8(0),x8(1),x8(2),x8(3)));
    pt->AddText("sam20_1  :    sam20_2   :   cur20_1    :  cur20_2");
    pt->AddText(Form("(%f, %f, %f, %f)",x20(0),x20(1),x20(2),x20(3)));
    pt->AddText("c   :     TO2    :    TO3     :   TO4");
    pt->AddText(Form("(%f, %f, %f, %f)",cur8c,TO[0],TO[1],TO[2]));
    pt->AddText(Form("r_tr   :     on    :    bg     :   sig"));
    pt->AddText(Form("(%f, %f, %f, %f)",r_tr,on,bg,sig));
    // pt->AddText(Form("(%f, %f, %f, %f)",-1.64157, 0.188066, 0.0843894, 0.000497484));
    pt->Draw();

    c->cd(2);
    map->Draw("colz");

    TFile *fout1 = TFile::Open(Form("./magic/raven%s_%f.root", id, sigma),"recreate");
    fout1->cd();
    ring->Write();
    Dgaus->Write();
    map->Write();
    map2->Write();   
    fout1->Close();

    // system(Form("mkdir ./magic/raven%s_%f",id,sigma));
    // char *rstarfit = Form("/home/lhaaso/caowy/2dfit_use/starfit ./magic/raven%s_%f.root -usehist -histname sky -outdir ./magic/raven%s_%f/ >./magic/raven%s_%f/raven%f.log -tillend",id, sigma,id,sigma,id,sigma,sigma);
    // // char *rstarfit = Form("/lhaasofs/user/gmxiang/Rc_crab/spectrum/20190901_20200229/2dfit_1/starfit ./magic/raven%s_skymap.root -usehist -histname sky -outdir ./magic/%s/ >./magic/%s/%s.log -tillend",id,id,id,id);

    // FILE *ff = popen(rstarfit, "r"); // build pipe
	// if (!ff)
	// 	return 1;

    // char tmp[1024];
	// while (fgets(tmp, sizeof(tmp), ff) != NULL)
	// 	std::cout << tmp << std::endl; // can join each line as string
	// pclose(ff);

    c -> SaveAs(Form("./magic/raven%s_%f.eps",id, sigma));

    ofstream outFile;
    outFile.open(Form("./magic/raven%s.csv",id), ios::app);

    outFile<< id << "," << x8[0] << ","  << x8[1] << ","  << x8[2] << ","  << x8[3] << "," << x20[0] << "," << x20[1] << "," << x20[2] << "," << x20[3] << "," << cur8c << "," << TO[0] << "," << TO[1] << "," << TO[2] << "," << sigma << "," << reduced_chi2 <<","<< hgaus+cstt << "," <<hgaus << "," << cstt << "," << sig <<endl;

    outFile.close();
    return 0;
}

int raven::starfit()
{    
    system(Form("mkdir ./magic/raven%s_%f",id,sigma));
    char *rstarfit = Form("/home/lhaaso/caowy/2dfit_use/starfit ./magic/raven%s_%f.root -usehist -histname sky -outdir ./magic/raven%s_%f/ >./magic/raven%s_%f.log -tillend",id, sigma, id,sigma,id,sigma);
    // char *rstarfit = Form("/lhaasofs/user/gmxiang/Rc_crab/spectrum/20190901_20200229/2dfit_1/starfit ./magic/raven%s_skymap.root -usehist -histname sky -outdir ./magic/%s/ >./magic/%s/%s.log -tillend",id,id,id,id);

    FILE *ff = popen(rstarfit, "r"); // build pipe
	if (!ff)
		return 1;

    char tmp[1024];
	while (fgets(tmp, sizeof(tmp), ff) != NULL)
		std::cout << tmp << std::endl; // can join each line as string
	pclose(ff);


    return 0;
}
