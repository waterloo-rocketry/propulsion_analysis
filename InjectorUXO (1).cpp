/*
 * rocket_sim.cxx
 *
 * Copyright 2016 Nicholas Christopher <nicholas@Mackenstein>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */


#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

double data_grab(string property,double value, string TorP,double (&a)[30][100], string (&header)[30]); //propertry you want, value of pressure or temperature, Constant Temp or Pressure intervals (use "P" or "T")
//callable values listed in the function





void data_gather(ifstream &infile,double (&a2)[30][100], string (&header2)[30]);


int main()
{
	/*
	 This code is used to determine when critical mass flow will occur in our injector, and to properly size the injector holes. It uses the Homogeneous Equilibrium Model.
	*/
	double g=9.81;//gravity
	double constT[30][100],constP[30][100];
        string headT[30],headP[30];

ifstream Pfile("N20_100_1000PSI.txt");
ifstream Tfile("N20_Neg30_35T.txt");

data_gather(Pfile, constP, headP);
cout <<endl;
data_gather(Tfile, constT, headT);
cout <<endl;
//****Parameters*********
	double T1=0;//degrees C
    for (int i=0;i<3;i++)
    {
        if(i==1)
            T1=15;
        else if(i==2)
            T1=20;
	double P1=data_grab("Pressure (psia)",T1,"T",constT, headT);//psi,absolute pre-injector pressure
	double P2=0;//psi, absolute post-injector pressure
	double P2F=500;//max downstream pressure
	double P2int=50;//step intervals


	double Cd=1;//
	double Cvc=1;//0.64;//vena contracta effect 0.64 for small L/D, 1 for large L/D
	double k=0.4;//local loss coeff
	double Ks = 0.0003;//roughness value (mm);
	double Aratio=0;//A1/A2

	double L=25.4;//orifice length (mm)
	double Di=2.18;//orifice diameter (mm)
	double N=30;//number of orifices

	//double L=3.175;//current injector
	//double Di=2.9464;//current injector

	double ul=0.0774;//dynamic viscosity at 0 deg mNs/m^2 (for determining estimate of reynolds number
	double ug=15.2;//dynamic visc at o deg uNs/m^2
	ul/=1000;
	ug/=1000000;

	Di/=1000;//(m)
	L/=1000;//(m)
	Ks/=1000;//(m)

	double A= N*M_PI*Di*Di/4;//m^3

	double f = pow(1/(2*log10(3.7/(Ks/Di))),2); //friction factor, assume complete turbulence. effectively 0 for small L/D


	Cd=Cvc/sqrt(1+N*(k+(f*L/Di))+Aratio*Aratio); //discharge coeff for large L/D (ie new injector)
	cout<<"P1: "<<P1<<endl;
	cout<<"D: "<<Di<<endl;
	cout<<"A: "<<A<<endl;
	cout<<"Cd: "<<Cd<<endl;
	cout<<"f: "<<f<<endl;

//****Unit conversion****
//	T1+=273;//Kelvin
//	P1*=6894.76;//Pascals
//	P2*=6894.76;//Pascals
double x2;
double h2;
double rho2;
double mHEM;
double v2;
double u2;
double Re;
	std::cout <<"P2"<< " "<< "x2" << " "<< "Re" <<" " <<"h1 rho"<< " "<< "mHEM@"<< P1<<"PSIA mSPI@"<< P1<< "PSIA mDYER@"<<P1<<"PSIA"<<endl;



	double s1l =data_grab("Entropy (l, J/g*K)",T1,"T",constT, headT);//entropy at 1
	double h1l =data_grab("Enthalpy (l, kJ/kg)",T1,"T",constT, headT);//enthalpy at 1
	h1l*=1000;
	double rhol1 =data_grab("Density (l, kg/m3)",T1,"T",constT, headT);//density at 1
//	cout<<s1l<<endl;
	while(P2<(P2F+1)){
	double s2l =data_grab("Entropy (l, J/g*K)",P2,"P",constP, headP);//liquid entropy at 2
	double s2g =data_grab("Entropy (v, J/g*K)",P2,"P",constP, headP);//gas entropy at 2
	double h2l =data_grab("Enthalpy (l, kJ/kg)",P2,"P",constP, headP);//liquid enthalpy at 2
	double h2g =data_grab("Enthalpy (v, kJ/kg)",P2,"P",constP, headP);//gas enthalpy at 2
	h2l*=1000;
	h2g*=1000;
	double rho2l =data_grab("Density (l, kg/m3)",P2,"P",constP, headP);//liquid density at 2
	double rho2g =data_grab("Density (v, kg/m3)",P2,"P",constP, headP);//gas density at 2
	x2 = (s1l-s2l)/(s2g-s2l);
	h2 = h2g*x2+(1-x2)*h2l;
	rho2 = 1/(x2/rho2g+(1-x2)/rho2l);
	mHEM = Cd*A*rho2*sqrt(2*(h1l-h2));
	v2 = mHEM/(rho2*A);
	u2 = ug*x2+(1-x2)*ul;
	Re = v2*Di*rho2/u2;

	double P1Pa=P1*6894.76;//Pascals
	double P2Pa=P2*6894.76;//Pascals
	double mSPI = Cd*A*sqrt(2*rhol1*(P1Pa-P2Pa));
    double kappa = 1;//sqrt((460-P2)/(520-P2)); //dyer et al factor
	double mDYER = ((mSPI*kappa)/(1.0+kappa) + (mHEM)/(1.0+kappa));
	std::cout <<P2<< " "<< x2 << " "<< Re <<" "<<h1l<<" "<<rho2<< " "<< mHEM<<" "<<mSPI<< " "<< mDYER<< endl;
	P2+=P2int;
}
    }
//	system("PAUSE");
	//cin>>Mt;//so the code doesn't end (couldnt get system("pause") to work
	return 0;
}




double data_grab(string property,double value, string TorP,double (&a)[30][100], string (&header)[30])
{

int ref_col;

if(TorP=="P")
{
ref_col =1;
}

else
{
ref_col =0;
}

int find1 =0;
int find2 =0;
while(property != header[find1])
{
find1++;
if (find1>30)
	return 0;
}

while(a[ref_col][find2]<value)
{
find2++;
}
//cout<< "value is at: " << a[ref_col][find2]<<endl;
return a[find1][find2];
}


/* List of all callable Values:

Temperature (C),Pressure (psia),Density (l, kg/m3),Volume (l, m3/kg)
Internal Energy (l, kJ/kg),Enthalpy (l, kJ/kg),Entropy (l, J/g*K)

Density (v, kg/m3),Volume (v, m3/kg),Internal Energy (v, kJ/kg),
Enthalpy (v, kJ/kg),Entropy (v, J/g*K)

Cv (l, J/g*K),Cp (l, J/g*K),Sound Spd. (l, m/s),Joule-Thomson (l, F/psia),Viscosity (l, uPa*s),Therm. Cond. (l, W/m*K),Surf. Tension (l, N/m),Cv (v, J/g*K),Cp (v, J/g*K)	Sound Spd. (v, m/s),Joule-Thomson (v, F/psia),Viscosity (v, uPa*s),Therm. Cond. (v, W/m*K)
*/


void data_gather(ifstream &infile, double (&a2)[30][100], string (&header2)[30])
{
string line;
getline(infile, line);
string delimiter = "\t";
size_t pos = 0;
int e=0;
while ((pos = line.find(delimiter)) != std::string::npos) {
    header2[e] = line.substr(0, pos);
    line.erase(0, pos + delimiter.length());
	e++;
}


int col = 0;
int row = 0;
string transfer;
while (getline(infile, line))
{
while (e>=col && (pos = line.find(delimiter)) != std::string::npos) {
    transfer = line.substr(0, pos);
    istringstream isss(transfer);
    isss >> a2[col][row];

    line.erase(0, pos + delimiter.length());
	col++;
}
col=0;

//cout<<a[col][row];
//cout<<endl;
    row++;
}
return;
}

