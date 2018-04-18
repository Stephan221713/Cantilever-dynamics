//This is a program to calculate the cantilever motion in dynamic contact mode (e.g. for ESM experiments). It was written by Stephan Bradler using a bit of previous code by Stefan Kachel.
//Version 4, 2018-04-18, Stephan Bradler. For details see Readme.pdf

//Header
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
using namespace std;

//Define functions
double ExtractParameter(string variable);
void calculate_A_and_phi();
int sgn(double sgnarg);
double q(double x);
void load_input();
void list_parameters();
void converge_DHO();
void calculate_DHO();
void calculate_resonance_curve();
void AdjustContactDamping();
void AdjustAirDamping();
void load_input_experimental();
void list_parameters_experimental();
void list_adjusted_parameters();
void AdjustMasses();
void Predict_ffree2();
void AdjustContactStiffness();
void calculate_sensitivities();

//Define global parameters

//input
double k, L, H, zeta, psi, w, VACxDV, fStart, fStop, fStep, Df, mT, mL, eta, k1, k2, d1, d2, f0, g1, g2, epsilon_r, TE, u0;
double Precision, AdjustDf, AdjustQAir, AdjustQContact, IntegrationSteps, LaserPosition, AdjustParametersFromExperiment, ResonanceCurve, ContactStiffnessVariation, LeverShape;

//input experimental
double ffree1, Qfree1, ffree2, Qfree2, fcr, Qcr, csr, Min_k1, Max_k1;

//calculate lever motion
double pi, eps0, qa, EI, theta, f, omega, A, A_max, f_max, phi, NIh, argu1, argu2, q0, q1, q2, qL;
complex<double> compi, alpha, k1c, k2c, c1L, c2L, c3L, c4L, kcxL, z0sL, z1sL, z2sL, z3sL, pa2, pa4, pa0, pb2, pb4, pb0;
complex<double> B1, B2, B4, z1L, z0L, z2L, z3L, c_an, Av_an, Ae_an, Aq_an, N_an, T, X, U, Fvert, Flon, cp, cm, sp, sm;
int NIn,NIi;

//DHO calculation
double f1, f2, A1, A2, phi1, phi2, Omega_DHO, Phi_G, X1, X2, f0_DHO, Q_DHO, A_drive, phi_drive, A_drive_mem, phi_drive_mem, f0_DHO_mem, Q_DHO_mem;
int convDHOcounter;

//Lever shape
double x, qx, xIh;
complex<double> c1x, c2x, c3x, c4x, z0sx, z1sx, z0x, z1x;
int xIi, xIn;

//Adjust damping parameters
double g1inp, g2inp, g1a, g2a, g1b, g2b, g2g1ratio, Qair, Qa, Qb, Slope, etainp, etaa, etab, Qcont;
int DLI;

//Adjust masses and contact stiffness
double EMR, m_eff, mTm_eff, mTa, mTb, ffree2a, ffree2b, SlopemT, mLa, mLb, ffree1a, ffree1b, SlopemL, k1a, k1b, fcra, fcrb, Slopek1;
int mTLI, mLLI, k1LI;

//calculate sensitivites
double k1mem, k2mem, d1mem, d2mem, f0mem, VACxDVmem, u0mem, StaticSens, f0_DHOmem, g1mem, g2mem, Q_DHOmem;

ofstream results_file;
bool silent=true;
int ParmLoopCounter, PLCAmp;




//-----------------------Main----------------------
int main(){
    
    //Create outputfile and header
    string dateiname = "Calculated cantilever dynamics.txt";
    results_file.open(dateiname.c_str(), ios::out);
    if(!results_file){
        cerr << "File could not be created." << endl;
        getchar();
    }
    results_file << "Results from the numerical calculation of the cantilever motion" << endl;
    
    load_input();
    list_parameters();
    
    pi = 3.141592653589793;
    compi = complex<double> (0,1.0);
    eps0= 8.854187817e-12;

    if (AdjustParametersFromExperiment>0){
        load_input_experimental();
        list_parameters_experimental();
        AdjustMasses();//Assuming only local electrostatic force
        AdjustContactStiffness();//Assuming only vertical sample displacement
        AdjustQAir=0;
        fStart=fcr*0.8;
        fStop=fcr*1.2;
        list_adjusted_parameters();}

    if (ResonanceCurve>0){
        silent=false;}
    if (d1==0&&d2==0&&f0==0&&VACxDV==0&&u0==0){d1=1;}
    calculate_resonance_curve();
    
    results_file << endl << "DHO Parameters and sensitivities" << endl << endl;
    if (ContactStiffnessVariation>0){
        results_file << "k1\tk2\tg1\tg2\t";}
    results_file << "A_drive\tphi_drive\tf0_DHO\tQ_DHO";
    results_file << "\tSensd1(pm/pm)\tSensd2(pm/pm)\tSensf0(pm/nN)\tSensq(pm/V^2)\tSensu0(pm/pm)";
    results_file << endl;
    
    //Begin parameter loop
    k1mem=k1;k2mem=k2;PLCAmp=0;g1mem=g1;g2mem=g2;
    if (ContactStiffnessVariation>0){
        PLCAmp=int(ContactStiffnessVariation+0.00000001);
        f0_DHO*=(1.0-ContactStiffnessVariation/100.0);
        }
    for (ParmLoopCounter=-PLCAmp; ParmLoopCounter<=PLCAmp; ParmLoopCounter+=1){
        k1=k1mem*pow(1.1,ParmLoopCounter);k2=k2mem*pow(1.1,ParmLoopCounter);g1=g1mem*pow(1.1,ParmLoopCounter);g2=g2mem*pow(1.1,ParmLoopCounter);
    
    f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
    if (AdjustQContact>0){
        if (AdjustQAir>0){cout << "Cannot adjust both air and contact damping, only adjusting contact damping" << endl; AdjustQAir=0;}
        AdjustContactDamping();}
    if (AdjustQAir>0){AdjustAirDamping();}
    
    converge_DHO();
    if (ContactStiffnessVariation>0){
        results_file << k1 << '\t' << k2 << '\t' << g1 << '\t' << g2 << '\t';}
    results_file << A_drive << '\t' << phi_drive << '\t' << f0_DHO << '\t' << Q_DHO;
    if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;g1=g1mem*pow(1.1,ParmLoopCounter);g2=g2mem*pow(1.1,ParmLoopCounter);}
    calculate_sensitivities();
    results_file << endl;
    
    }//End parameter loop
    k1=k1mem;k2=k2mem;
    if (AdjustQContact>0){
        if (AdjustQAir>0){cout << "Cannot adjust both air and contact damping, only adjusting contact damping" << endl; AdjustQAir=0;}
        AdjustContactDamping();}
    if (AdjustQAir>0){AdjustAirDamping();}
    converge_DHO();
    
    if (LeverShape>0){
        results_file << endl << "Lever shape at the resonance frequency" << endl << endl << "x\tAz0\tphiz0\tAz1\tphiz1\tShape" << endl;
        for (LaserPosition=0; LaserPosition<=1.0001; LaserPosition+=0.01){
            omega=f0_DHO*2.0*pi;
            calculate_A_and_phi();
            results_file << x << '\t' << abs(z0x) << '\t' << -arg(z0x) << '\t' << A << '\t' << phi << '\t' << -imag(z0x) << endl;
        }
    }
    
    results_file.close();
    cout << "Calculation finished" << endl;
    return EXIT_SUCCESS;                //end program
}//end main

//-----------------------Functions----------------------
void load_input(){
     k = ExtractParameter("k");
     L = ExtractParameter("L");
     H = ExtractParameter("H");
     zeta = ExtractParameter("zeta");
     psi = ExtractParameter("psi");
     w = ExtractParameter("w");
     VACxDV = ExtractParameter("VACxDV");
     fStart = ExtractParameter("fStart");
     fStop = ExtractParameter("fStop");
     fStep = ExtractParameter("fStep");
     Df = ExtractParameter("Df");
     mT = ExtractParameter("mT");
     mL = ExtractParameter("mL");
     eta = ExtractParameter("eta");
     k1 = ExtractParameter("k1");
     k2 = ExtractParameter("k2");
     d1 = ExtractParameter("d1");
     d2 = ExtractParameter("d2");
     f0 = ExtractParameter("f0");
     u0 = ExtractParameter("u0");
     g1 = ExtractParameter("g1");
     g2 = ExtractParameter("g2");
     epsilon_r = ExtractParameter("epsilon_r");
     TE = ExtractParameter("TE");
     Precision = ExtractParameter("Precision");
     AdjustDf = ExtractParameter("AdjustDf");
     AdjustQAir = ExtractParameter("AdjustQAir");
     AdjustQContact = ExtractParameter("AdjustQContact");
     IntegrationSteps = ExtractParameter("IntegrationSteps");
     LaserPosition = ExtractParameter("LaserPosition");
     AdjustParametersFromExperiment = ExtractParameter("AdjustParametersFromExperiment");
     ResonanceCurve = ExtractParameter("ResonanceCurve");
     ContactStiffnessVariation = ExtractParameter("ContactStiffnessVariation");
     LeverShape = ExtractParameter("LeverShape");
}//end load_input

double ExtractParameter(string variable){//Read parameter values from input file. Written by Stefan Kachel.
     string number;
     bool found=false;
     fstream file("Input.txt",ios::in);

     while(file >> number){
     if(!number.compare(variable)){
          file >> number;
          found=true;
          stringstream stream;
          stream << number;
          double number; 
          stream >> number;
          return number;
          }
     }
     if(!found) cout << "Parameter not found!\n";
     file.close();
}//end ExtractParameter

void list_parameters(){
     results_file << endl << "Input" << endl << endl << "Parameter\tvalue in SI units" << endl;
     results_file << "L\t" << L <<endl;
     results_file << "H\t" << H <<endl;
     results_file << "w\t" << w <<endl;
     results_file << "zeta\t" << zeta <<endl;
     results_file << "psi\t" << psi <<endl;
     results_file << "k\t" << k <<endl;
     results_file << "mT\t" << mT <<endl;
     results_file << "mL\t" << mL <<endl;
     results_file << "eta\t" << eta <<endl;
     results_file << "k1\t" << k1 <<endl;
     results_file << "k2\t" << k2 <<endl;
     results_file << "g1\t" << g1 <<endl;
     results_file << "g2\t" << g2 <<endl;
     results_file << "d1\t" << d1 <<endl;
     results_file << "d2\t" << d2 <<endl;
     results_file << "f0\t" << f0 <<endl;
     results_file << "VACxDV\t" << VACxDV <<endl;
     results_file << "u0\t" << u0 <<endl;
     results_file << "epsilon_r\t" << epsilon_r <<endl;
     results_file << "LaserPosition\t" << LaserPosition <<endl;
     results_file << "Df\t" << Df <<endl;
     results_file << "TE\t" << TE <<endl << endl;
     results_file << "fStart\t" << fStart <<endl;
     results_file << "fStop\t" << fStop <<endl;
     results_file << "fStep\t" << fStep <<endl;
     results_file << "AdjustDf\t" << AdjustDf <<endl;
     results_file << "AdjustQAir\t" << AdjustQAir <<endl;
     results_file << "AdjustQContact\t" << AdjustQContact <<endl;
     results_file << "Precision\t" << Precision <<endl;
     results_file << "IntegrationSteps\t" << IntegrationSteps <<endl << endl;
     results_file << "AdjustParametersFromExperiment\t" << AdjustParametersFromExperiment <<endl;
     results_file << "ResonanceCurve\t" << ResonanceCurve <<endl;
     results_file << "ContactStiffnessVariation\t" << ContactStiffnessVariation <<endl;
     results_file << "LeverShape\t" << LeverShape <<endl;
}//end list_parameters

void calculate_resonance_curve(){//Loops through frequencies between fStart and fStop. Finds the resonance peak.
     A_max=0;
     if(!silent){results_file << endl << "Amplitude and phase as a function of frequency" << endl << endl << "f\tA\tphi" << endl;}
     for(f=fStart; f<=fStop; f+=fStep){
          omega=f*2.0*pi;
          calculate_A_and_phi();
          if(!silent){results_file << f << '\t' << A << '\t' << phi << endl;}
          if(A>A_max){
               A_max=A;
               f_max=f;
          }
     }
     f1=f_max-Df/2.0;
     f2=f1+Df;
     calculate_DHO();
}//end calculate_resonance_curve

void calculate_A_and_phi(){//This calculates the lever slope at the laser position for a given set of input parameters. This is where most of the math happens.

     //Step 1: Some precalculations
          theta = zeta+psi;
          qa = -eps0*epsilon_r*w*VACxDV*cos(zeta);
          NIn = int(IntegrationSteps);
          EI = k*L*L*L/3.0;
          alpha=complex<double>(1.0,-eta/omega);
          alpha*=3.0*omega*omega/L/L/L/L/k*mL;
          alpha=pow(alpha,0.25);
          k1c=complex<double> (k1,omega*g1);
          k2c=complex<double> (k2,omega*g2);
       
     //Step 2: Calculate solution for z0(0)=z1(0)=z2(0)=z3(0)=0

          //using numerical integration via Simpson method to calculate coefficients cxL
          if(VACxDV!=0){//only necessary if non-local electrostatic forces are present
          NIh=L/NIn;
          q0=q(0.0);
          qL=q(L);
          c1L=exp(-alpha*L)*qL-q0;
          c2L=compi*(exp(-compi*alpha*L)*qL-q0);
          c3L=-(exp(alpha*L)*qL-q0);
          c4L=-compi*(exp(compi*alpha*L)*qL-q0);
          for(NIi=0; NIi<NIn; NIi+=1){
               argu1=NIh*NIi+NIh/2.0;
               argu2=NIh*NIi;
               q1=q(argu1);
               q2=q(argu2);
               c1L+=4.0*exp(-alpha*argu1)*q1+2.0*exp(-alpha*argu2)*q2;
               c2L+=compi*(4.0*exp(-compi*alpha*argu1)*q1+2.0*exp(-compi*alpha*argu2)*q2);
               c3L+=-(4.0*exp(alpha*argu1)*q1+2.0*exp(alpha*argu2)*q2);
               c4L+=-compi*(4.0*exp(compi*alpha*argu1)*q1+2.0*exp(compi*alpha*argu2)*q2);
          }
          kcxL=NIh/24.0/alpha/alpha/alpha/EI;//common factor for cxL
          c1L*=kcxL; c2L*=kcxL; c3L*=kcxL; c4L*=kcxL;

          //Calculate values for x=L of the first inhomogenous solution
          z0sL=c1L*exp(alpha*L)+c2L*exp(compi*alpha*L)+c3L*exp(-alpha*L)+c4L*exp(-compi*alpha*L);
          z1sL=alpha*(c1L*exp(alpha*L)+c2L*compi*exp(compi*alpha*L)-c3L*exp(-alpha*L)-c4L*compi*exp(-compi*alpha*L));
          z2sL=alpha*alpha*(c1L*exp(alpha*L)-c2L*exp(compi*alpha*L)+c3L*exp(-alpha*L)-c4L*exp(-compi*alpha*L));
          z3sL=alpha*alpha*alpha*(c1L*exp(alpha*L)-c2L*compi*exp(compi*alpha*L)-c3L*exp(-alpha*L)+c4L*compi*exp(-compi*alpha*L));
          } else {z0sL=0; z1sL=0; z2sL=0; z3sL=0;}//for VACxDV=0, the specific solution z0sL and its derivatives are 0

     //Step 3: Calculate B2 and B4 from pa2*B2+pa4*B4+pa0=0 and pb2*B2+Pb4*B4+Pb0=0 using the boundary conditions

          T=H*H/EI*(k1c*sin(theta)*sin(theta)+k2c*cos(theta)*cos(theta))-H*H/16.0*mT/mL*L*pow(alpha,4.0);
          X=H/EI*sin(theta)*cos(theta)*(k2c-k1c);
          U=1.0/EI*(k1c*cos(theta)*cos(theta)+k2c*sin(theta)*sin(theta))-mT/mL*L*pow(alpha,4.0);
          Fvert=k1c*d1+f0;
          Flon=k2c*d2;
          B1=u0/2.0;
       
          //A short notation for mixed trigonometric functions
          cp=cos(alpha*L)+cosh(alpha*L);
          cm=cos(alpha*L)-cosh(alpha*L);
          sp=sin(alpha*L)+sinh(alpha*L);
          sm=sin(alpha*L)-sinh(alpha*L);
       
          //These are the actual boundary conditions from which B2 and B4 are calculated
          pa0=-z2sL+H/EI*cos(theta)*Flon-T*z1sL-H/EI*sin(theta)*Fvert-X*z0sL+B1*alpha*alpha*cm+B1*T*alpha*sm-B1*X*cp;
          pa2=alpha*alpha*cp+T*alpha*sp-X*cm;
          pa4=alpha*alpha*sp-T*alpha*cm-X*sm;
          pb0=-z3sL+U*z0sL-1.0/EI*cos(theta)*Fvert+X*z1sL-1.0/EI*sin(theta)*Flon-B1*alpha*alpha*alpha*sp+B1*U*cp-B1*X*alpha*sm;
          pb2=-alpha*alpha*alpha*sm+U*cm-X*alpha*sp;
          pb4=alpha*alpha*alpha*cp+U*sm+X*alpha*cm;
          B2=(pb0/pb4-pa0/pa4)/(pa2/pa4-pb2/pb4);
          B4=(pb0/pb2-pa0/pa2)/(pa4/pa2-pb4/pb2);

     //Step 4: Calculate z1L and from that A and phi

          //Numerical solution
          z0L=z0sL+B1*cp+B2*cm+B4*sm;
          z1L=z1sL+alpha*(-B1*sm-B2*sp+B4*cm);
          z2L=z2sL+alpha*alpha*(-B1*cm-B2*cp-B4*sp);
          z3L=z3sL+alpha*alpha*alpha*(B1*sp+B2*sm-B4*cp);
          A=abs(z1L);
          phi=-arg(z1L);
          
     //Step 5 (optional): Calculate the lever displacement and its derivative at the laser spot position
     //Note: This is a flexible way of calculating the lever displacement and slope at the laser spot position.
     //If you want to calculate the lever shape (displacement and slope for many x values), this is not computationally efficient.
     //It would be better to integrate to the first x value and reuse this integration as a starting point for the next x value.
     //The calculation is still reasonably fast, so I did not include an optimized code for the lever shape calculation.
     //Also, this only includes the absolute displacement and its first derivative. For higher derivatives just adapt the last lines of step 2 as well as step 4.
               
          if(LaserPosition!=1.0){
          xIn=int(NIn*LaserPosition);
          x=L*LaserPosition;
          if(VACxDV!=0){
          xIh=x/xIn;
          q0=q(0.0);
          qx=q(x);
          c1x=exp(-alpha*x)*qx-q0;
          c2x=compi*(exp(-compi*alpha*x)*qx-q0);
          c3x=-(exp(alpha*x)*qx-q0);
          c4x=-compi*(exp(compi*alpha*x)*qx-q0);
          for(xIi=0; xIi<xIn; xIi+=1){
               argu1=xIh*xIi+xIh/2.0;
               argu2=xIh*xIi;
               q1=q(argu1);
               q2=q(argu2);
               c1x+=4.0*exp(-alpha*argu1)*q1+2.0*exp(-alpha*argu2)*q2;
               c2x+=compi*(4.0*exp(-compi*alpha*argu1)*q1+2.0*exp(-compi*alpha*argu2)*q2);
               c3x+=-(4.0*exp(alpha*argu1)*q1+2.0*exp(alpha*argu2)*q2);
               c4x+=-compi*(4.0*exp(compi*alpha*argu1)*q1+2.0*exp(compi*alpha*argu2)*q2);
          }
          c1x*=kcxL; c2x*=kcxL; c3x*=kcxL; c4x*=kcxL;
          z0sx=c1x*exp(alpha*x)+c2x*exp(compi*alpha*x)+c3x*exp(-alpha*x)+c4x*exp(-compi*alpha*x);
          z1sx=alpha*(c1x*exp(alpha*x)+compi*c2x*exp(compi*alpha*x)-c3x*exp(-alpha*x)-compi*c4x*exp(-compi*alpha*x));
          }else{z0sx=0;z1sx=0;}//We can skip the integration for VACxDV=0.
          z0x=z0sx+B1*(cos(alpha*x)+cosh(alpha*x))+B2*(cos(alpha*x)-cosh(alpha*x))+B4*(sin(alpha*x)-sinh(alpha*x));
          z1x=z1sx-B1*alpha*(sin(alpha*x)-sinh(alpha*x))-B2*alpha*(sin(alpha*x)+sinh(alpha*x))+B4*alpha*(cos(alpha*x)-cosh(alpha*x));
          if(LaserPosition==0){z0x=complex<double> (u0,0);z1x=complex<double> (0,0);}
          A=abs(z1x);
          phi=-arg(z1x);
          }

}//end calculate_A_and_phi

double q(double x){//Calculates q(x)
     return qa*pow(H*cos(zeta)+(L-x)*sin(zeta),-2.0);
}//end q

void calculate_DHO(){//Performs a theoretical DART analysis at the frequencies f1 and f2. Formulas from Gannepalli et al. 2011; doi:10.1088/0957-4484/24/15/159501
     omega=f1*2.0*pi;
     calculate_A_and_phi();
     A1=A;
     phi1=phi;
     omega=f2*2.0*pi;
     calculate_A_and_phi();
     A2=A;
     phi2=phi;
     Omega_DHO = f1*A1/(f2*A2);
     Phi_G = tan(phi2 - phi1);
     X1 = (-1+sgn(Phi_G)*sqrt(1+Phi_G*Phi_G)/Omega_DHO)/Phi_G;
     X2 = (1-sgn(Phi_G)*Omega_DHO*sqrt(1+Phi_G*Phi_G))/Phi_G;
     f0_DHO = sqrt(f1*f2*(f2*X1-f1*X2)/(f1*X1-f2*X2));
     Q_DHO = sqrt(f1*f2*(f2*X1-f1*X2)*(f1*X1-f2*X2))/(f2*f2-f1*f1);
     A_drive = A1*sqrt(pow((f0_DHO*f0_DHO-f1*f1),2)+pow((f0_DHO*f1/Q_DHO),2))/(f0_DHO*f0_DHO);
     phi_drive = phi1-atan2((f0_DHO*f1),(Q_DHO*(f0_DHO*f0_DHO-f1*f1)));
}//end calculate_DHO

int sgn(double sgnarg){
     if(sgnarg>0){return 1;}
     else if (sgnarg<0){return (-1);}
     else{return 0;};
}//end sgn

void converge_DHO(){//Optimizes test frequencies f1 and f2 for optimal DART analysis.
     for(convDHOcounter=0; convDHOcounter<100; convDHOcounter+=1){ 
          A_drive_mem=A_drive; phi_drive_mem=phi_drive; f0_DHO_mem=f0_DHO; Q_DHO_mem=Q_DHO;
          if (AdjustDf!=0){Df=f0_DHO/sqrt(2.0)/Q_DHO*AdjustDf;}
          f1=f0_DHO-Df/2.0; f2=f0_DHO+Df/2.0;
          calculate_DHO();
          if (abs((A_drive-A_drive_mem)/A_drive_mem)<Precision &&
              abs((phi_drive-phi_drive_mem)/phi_drive_mem)<Precision &&
              abs((f0_DHO-f0_DHO_mem)/f0_DHO_mem)<Precision &&
              abs((Q_DHO-Q_DHO_mem)/Q_DHO_mem)<Precision)
              {convDHOcounter+=200;}
          else if (convDHOcounter==99){cout << "DHO parameters did not converge!" << endl;}
          if (TE!=0){f1=f0_DHO-Df/2.0+TE; f2=f0_DHO+Df/2.0+TE;calculate_DHO();}
     }
}//end converge_DHO

void AdjustContactDamping(){//Automatically adjusts g1 and g2 to match predefined Q factor AdjustQContact.
     if (eta==0){//This facilitates the optimization procedure.
            if (g1==0 && g2==0){//Initial guesses need to be calculated.
                g1=k1/(2.0*pi*f0_DHO*AdjustQContact);
                g2=g1*k2/k1;
                silent=true;
                calculate_resonance_curve();//With no damping input, no peak was found yet, so we need to find a peak here.
                silent=false;
                }
            converge_DHO();
            for (DLI=0; DLI<100; DLI+=1){//Optimizing g1 and g2;
                g1*=Q_DHO/AdjustQContact;
                g2*=Q_DHO/AdjustQContact;
                converge_DHO();
                if (abs((Q_DHO-AdjustQContact)/AdjustQContact)<Precision){DLI+=200;cout<<"Contact damping parameters converged"<<endl;}
                if (DLI==99){cout << "Contact damping parameters did not converge" << endl;}
                }
        }
        else {
            g1inp=g1; g2inp=g2; g1=0; g2=0;
            converge_DHO();
            if (Q_DHO<AdjustQContact){cout << "eta is too high for AdjustQContact, setting eta to 0" << endl; eta=0;
            //Contact damping only increases total damping, so air damping may not be stronger than desired total damping. In this case, eta is set to 0 instead and the calculation is the same as above.
                                      g1=g1inp; g2=g2inp;
                                      if (g1==0 && g2==0){
                                          g1=k1/(2.0*pi*f0_DHO*AdjustQContact);
                                          g2=g1*k2/k1;
                                          silent=true;
                                          calculate_resonance_curve();
                                          silent=false;
                                          }
                                      converge_DHO();
                                      for (DLI=0; DLI<100; DLI+=1){
                                          g1*=Q_DHO/AdjustQContact;
                                          g2*=Q_DHO/AdjustQContact;
                                          converge_DHO();
                                          if (abs((Q_DHO-AdjustQContact)/AdjustQContact)<Precision){DLI+=200;cout<<"Contact damping parameters converged"<<endl;}
                                          if (DLI==99){cout << "Contact damping parameters did not converge" << endl;}
                                          }
            }
            else{
                 Qair=Q_DHO;
                 g1=g1inp; g2=g2inp;
                 if (g1==0 && g2==0){//Generating initial guess, if not given in input.
                     g1=k1/(2.0*pi*f0_DHO*AdjustQContact);
                     g2=g1*k2/k1;
                     }
                 g1a=g1; g2a=g2;
                 g2g1ratio=g2/g1;
                 converge_DHO();
                 Qa=Q_DHO;
                 g1=g1a*(1.0/AdjustQContact-1.0/Qair)/(1.0/Qa-1.0/Qair);//Generating second guess using 1/Qtotal=1/Qair+1/Qcontact
                 g2=g1*g2g1ratio;
                 converge_DHO();
                 Qb=Q_DHO;
                 g1b=g1; g2b=g2;
                 Slope=(Qb-Qa)/(g1b-g1a);
                 for (DLI=0; DLI<100; DLI+=1){//Optimizing contact damping
                     g1=g1a+(AdjustQContact-Qa)/Slope;
                     g2=g1*g2g1ratio;
                     g1a=g1; g2a=g2;
                     converge_DHO();
                     Qa=Q_DHO;
                     g1=g1a+(AdjustQContact-Qa)/Slope;
                     g2=g1*g2g1ratio;
                     g1b=g1; g2b=g2;
                     converge_DHO();
                     Qb=Q_DHO;
                     Slope=(Qb-Qa)/(g1b-g1a);
                     if (abs((Q_DHO-AdjustQContact)/AdjustQContact)<Precision){DLI+=200;cout<<"Contact damping parameters converged"<<endl;}
                     if (DLI==99){cout << "Contact damping parameters did not converge" << endl;}
                 }
            }
        }
}//end AdjustContactDamping

void AdjustAirDamping(){//Adjusting eta such that experimental Q value is reproduced. Very similar to AdjustContactDamping.
     if (g1==0 && g2==0){//Simple case: no contact damping.
            if (eta==0){//Initial guess not given, needs to be calculated.
                eta=2.0*pi*f0_DHO/AdjustQAir;
                silent=true;
                calculate_resonance_curve();//With no damping input, no peak was found yet, so we need to find a peak here.
                silent=false;
                }
            converge_DHO();
            for (DLI=0; DLI<100; DLI+=1){//Optimizing eta.
                eta*=Q_DHO/AdjustQAir;
                converge_DHO();
                if (abs((Q_DHO-AdjustQAir)/AdjustQAir)<Precision){DLI+=200;cout<<"Eta converged"<<endl;}
                if (DLI==99){cout << "Eta did not converge" << endl;}
                }
        }
        else {
            etainp=eta; eta=0;
            converge_DHO();
            if (Q_DHO<AdjustQAir){cout << "g1 and g2 are too high for AdjustQAir, setting g1 and g2 to 0" << endl; g1=0; g2=0;
            //Air damping only increases total damping, so contact damping may not be stronger than desired total damping. In this case, g1 and g2 are set to 0 instead and the calculation is the same as above.
                                      eta=etainp;
                                      if (eta==0){
                                          eta=2.0*pi*f0_DHO/AdjustQAir;
                                          silent=true;
                                          calculate_resonance_curve();
                                          silent=false;
                                          }
                                      converge_DHO();
                                      for (DLI=0; DLI<100; DLI+=1){
                                          eta*=Q_DHO/AdjustQAir;
                                          converge_DHO();
                                          cout << Q_DHO << endl;
                                          if (abs((Q_DHO-AdjustQAir)/AdjustQAir)<Precision){DLI+=200;cout<<"Eta converged"<<endl;}
                                          if (DLI==99){cout << "Eta did not converge" << endl;}
                                          }
            }
            else{
                 Qcont=Q_DHO;
                 eta=etainp;
                 if (eta==0){//If no initial eta is given, we need to calculated a first guess here.
                     eta=2.0*pi*f0_DHO/AdjustQAir;
                     silent=true;
                     calculate_resonance_curve();
                     silent=false;
                     }
                 etaa=eta;
                 converge_DHO();
                 Qa=Q_DHO;
                 eta=etaa*(1.0/AdjustQAir-1.0/Qcont)/(1.0/Qa-1.0/Qcont);//Generating second guess using 1/Qtotal=1/Qair+1/Qcontact
                 converge_DHO();
                 Qb=Q_DHO;
                 etab=eta;
                 Slope=(Qb-Qa)/(etab-etaa);
                 for (DLI=0; DLI<100; DLI+=1){//Optimzing eta.
                     eta=etaa+(AdjustQAir-Qa)/Slope;
                     etaa=eta;
                     converge_DHO();
                     Qa=Q_DHO;
                     eta=etaa+(AdjustQAir-Qa)/Slope;
                     etab=eta;
                     converge_DHO();
                     Qb=Q_DHO;
                     Slope=(Qb-Qa)/(etab-etaa);
                     if (abs((Q_DHO-AdjustQAir)/AdjustQAir)<Precision){DLI+=200;cout<<"Eta converged"<<endl;}
                     if (DLI==99){cout << "Eta did not converge" << endl;}
                 }
            }
        }
}//end AdjustAirDamping

void load_input_experimental(){
     ffree1 = ExtractParameter("ffree1");
     Qfree1 = ExtractParameter("Qfree1");
     ffree2 = ExtractParameter("ffree2");
     Qfree2 = ExtractParameter("Qfree2");
     fcr = ExtractParameter("fcr");
     Qcr = ExtractParameter("Qcr");
     csr = ExtractParameter("csr");
     Min_k1 = ExtractParameter("Min_k1");
     Max_k1 = ExtractParameter("Max_k1");
}//end load_input_experimental

void list_parameters_experimental(){
     results_file << endl << "Input_Experimental" << endl << endl << "Parameter\tvalue in SI units" << endl;
     results_file << "ffree1\t" << ffree1 <<endl;
     results_file << "Qfree1\t" << Qfree1 <<endl;
     results_file << "ffree2\t" << ffree2 <<endl;
     results_file << "Qfree2\t" << Qfree2 <<endl;
     results_file << "fcr\t" << fcr <<endl;
     results_file << "Qcr\t" << Qcr <<endl;
     results_file << "csr\t" << csr <<endl;
     results_file << "Min_k1\t" << Min_k1 <<endl;
     results_file << "Max_k1\t" << Max_k1 <<endl;
}//end list_parameters_experimental

void list_adjusted_parameters(){
     results_file << endl << "Adjusted parameters" << endl << endl << "Parameter\tvalue in SI units" << endl;
     results_file << "mT\t" << mT <<endl;
     results_file << "mL\t" << mL <<endl;
     results_file << "eta\t" << eta <<endl;
     results_file << "k1\t" << k1 <<endl;
     results_file << "k2\t" << k2 <<endl;
     results_file << "g1\t" << g1 <<endl;
     results_file << "g2\t" << g2 <<endl;
}//end list_parameters

void AdjustMasses(){//Adjusts masses and eta to match first two free resonance frequencies. Eta may be different for both cases. The eta for the first free resonance frequency is used for subsequent calculation.
     EMR=0.24267204503035;//setting a few parameters
     m_eff=k/(2.0*pi*ffree1)/(2.0*pi*ffree1);
     omega=2.0*pi*ffree2;
     mTa=0; mTb=0; ffree2a=0; ffree2b=0;
     k1=0; k2=0; g1=0; g2=0;
     d1mem=d1; d2mem=d2; f0mem=f0; VACxDVmem=VACxDV; u0mem=u0;
     d1=0; d2=0; f0=-1.0e-10; VACxDV=0; u0=0;
     eta=2.0*pi*ffree2/Qfree2;
     for(mTm_eff=0; mTm_eff<=0.5; mTm_eff+=0.01){//looping through combinations of mT and mL that should give the correct first resonance frequency.
                                                           //and testing if the second resonance frequency also matches. This gives initial guesses for the tip mass.
          mT=mTm_eff*m_eff;
          mL=m_eff*(1.0-mTm_eff)/EMR;
          calculate_A_and_phi();
          if (A>ffree2a){
               ffree2b=ffree2a;
               mTb=mTa;
               ffree2a=A;
               mTa=mTm_eff;
          }else if (A>ffree2b){
               ffree2b=A;
               mTb=mTm_eff;
          }else{}
     }
     
     //Exact calculation of the second resonance frequency from the tip mass (see Predict_ffree2()) is used to optimize tip mass.
     mTa*=m_eff; mTb*=m_eff;
     mT=mTa;
     Predict_ffree2();
     ffree2a=f0_DHO;
     mT=mTb;
     Predict_ffree2();
     ffree2b=f0_DHO;
     SlopemT=(ffree2b-ffree2a)/(mTb-mTa);
     for(mTLI=0; mTLI<100; mTLI+=1){
          mT=mTa+(ffree2-ffree2a)/SlopemT;
          mTa=mT;
          Predict_ffree2();
          ffree2a=f0_DHO;
          mT=mTa+(ffree2-ffree2a)/SlopemT;
          mTb=mT;
          cout << mT << endl;
          Predict_ffree2();
          ffree2b=f0_DHO;
          SlopemT=(ffree2b-ffree2a)/(mTb-mTa);
          if (abs((f0_DHO-ffree2)/ffree2)<Precision){mTLI+=200;cout<<"Tip mass converged"<<endl;}
          if (mTLI==99){cout << "Tip mass did not converge" << endl;}
     }
     f0_DHO=ffree1;
     Q_DHO=Qfree1;
     AdjustQAir=Qfree1;
     AdjustAirDamping();//Using the eta at the first resonance frequency for subsequent calculations.
     d1=d1mem; d2=d2mem; f0=f0mem; VACxDV=VACxDVmem; u0=u0mem;
}//end AdjustMasses

void Predict_ffree2(){//This uses a given tip mass, then adjusts the lever mass so that the first free resonance frequency is matched.
                             //Afterwards the resulting second free resonance frequency is calculated. At every step, eta is optimized to match the experimental Q factor.
     mLa=(m_eff-mT)/EMR;//First guess.
     mLb=mLa*(1.0-1.0/Qfree1);//Second guess.
     f0_DHO=ffree1;//Setting parameters for first free resonance frequency.
     Q_DHO=Qfree1;
     AdjustQAir=Qfree1;
     mL=mLa;
     AdjustAirDamping();
     ffree1a=f0_DHO;
     mL=mLb;
     AdjustAirDamping();
     ffree1b=f0_DHO;
     SlopemL=(ffree1b-ffree1a)/(mLb-mLa);
     for(mLLI=0; mLLI<100; mLLI+=1){//Optimzing mL.
          mL=mLa+(ffree1-ffree1a)/SlopemL;
          mLa=mL;
          AdjustAirDamping();
          ffree1a=f0_DHO;
          mL=mLa+(ffree1-ffree1a)/SlopemL;
          mLb=mL;
          AdjustAirDamping();
          ffree1b=f0_DHO;
          SlopemL=(ffree1b-ffree1a)/(mLb-mLa);
          if (abs((f0_DHO-ffree1)/ffree1)<Precision){mLLI+=200;cout<<"Lever mass converged"<<endl;}
          if (mLLI==99){cout << "Lever mass did not converge" << endl;}
     }
     f0_DHO=ffree2;//Calculating the second free resonance frequency
     Q_DHO=Qfree2;
     AdjustQAir=Qfree2;
     AdjustAirDamping();
}//end Predict_ffree2

void AdjustContactStiffness(){//Adjusts contact parameters to match experimental contact resonance properties.
     d1mem=d1; d2mem=d2; f0mem=f0; VACxDVmem=VACxDV; u0mem=u0;
     d1=1.0e-12; d2=0; f0=0; VACxDV=0; u0=0;
     omega=2.0*pi*fcr;
     k1a=0; k1b=0; fcra=0; fcrb=0;
     for(k1=Min_k1; k1<=Max_k1; k1*=1.2){//Generating initial guesses.
          k2=k1*csr;
          calculate_A_and_phi();
          if (A>fcra){
               fcrb=fcra;
               k1b=k1a;
               fcra=A;
               k1a=k1;
          }else if (A>fcrb){
               fcrb=A;
               k1b=k1;
          }else{}
     }
     f0_DHO=fcr;
     Q_DHO=Qcr;
     AdjustQContact=Qcr;
     g1=k1a/(2.0*pi*f0_DHO*AdjustQContact);
     g2=g1*csr;
     k1=k1a;
     k2=k1*csr;
     AdjustContactDamping();
     fcra=f0_DHO;
     k1=k1b;
     k2=k1*csr;
     AdjustContactDamping();
     fcrb=f0_DHO;
     Slopek1=(fcrb-fcra)/(k1b-k1a);
     for(k1LI=0; k1LI<100; k1LI+=1){//Optimizing k1 and k2. At every step, g1 and g2 are optimized.
          k1=k1a+(fcr-fcra)/Slopek1;
          k2=k1*csr;
          k1a=k1;
          AdjustContactDamping();
          fcra=f0_DHO;
          k1=k1a+(fcr-fcra)/Slopek1;
          k2=k1*csr;
          k1b=k1;
          cout << k1 << endl;
          AdjustContactDamping();
          fcrb=f0_DHO;
          Slopek1=(fcrb-fcra)/(k1b-k1a);
          if (abs((f0_DHO-fcr)/fcr)<Precision){k1LI+=200;cout<<"Contact stiffness converged"<<endl;}
          if (k1LI==99){cout << "Contact stiffness did not converge" << endl;}
     }
     d1=d1mem; d2=d2mem; f0=f0mem; VACxDV=VACxDVmem; u0=u0mem;
}//end AdjustContactStiffness

void calculate_sensitivities(){
     StaticSens=1.5e-12/L/cos(zeta)*LaserPosition*(2.0-LaserPosition);//It is assumed that the laser spot positioning is not changed between calibration and measurement.
     d1mem=d1; d2mem=d2; f0mem=f0; VACxDVmem=VACxDV; u0mem=u0;
     d1=1.0e-12; d2=0; f0=0; VACxDV=0; u0=0;
     f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
     converge_DHO();
     if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;}
     results_file << '\t' << A_drive/StaticSens;
     d1=0; d2=1.0e-12;
     f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
     converge_DHO();
     if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;}
     results_file << '\t' << A_drive/StaticSens;
     d2=0; f0=-1.0e-9;
     f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
     converge_DHO();
     if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;}
     results_file << '\t' << A_drive/StaticSens;
     f0=0; VACxDV=1.0;
     f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
     converge_DHO();
     if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;}
     results_file << '\t' << A_drive/StaticSens;
     VACxDV=0; u0=1.0e-12;
     f0_DHOmem=f0_DHO;Q_DHOmem=Q_DHO;
     converge_DHO();
     if(f0_DHO>0){}else{f0_DHO=f0_DHOmem;Q_DHO=Q_DHOmem;}
     results_file << '\t' << A_drive/StaticSens;
     d1=d1mem; d2=d2mem; f0=f0mem; VACxDV=VACxDVmem; u0=u0mem;
}//end calculate_sensitivities
