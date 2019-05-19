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

using namespace std;

int main()
{
	/*
	 * This is a simple program to approximate nozzle parameters based on equations in the book "Rocket Propulsion Elements" 
	 * 
	*/
	double g=9.81;//gravity
	double OF=5.5;//mass of oxidizer / mass of fuel grain
	double k=1.231;//cp/cv
    double T=2910;//degrees K
	double P=500;//psi,absolute chamber pressure
    double Pe=14.7;//psi,absolute exit pressure
    double R=.342;//KJ/(kg*k)
    double rho = 3.447;//kg/m^3
	double Dt=1.60;//nozzle throat diameter (inch)
	double At=0.001297*Dt*Dt*0.25*M_PI;//nozzle throat area (m^2)
    double Cf =1.4680;
    R*=1000;
	P*=6894.76;//Pascals
    Pe*=6894.76;//Pascals

	double numerator=pow(2/(k+1),(k+1)/(k-1));
	double denom=k*R*T;
    double a = 1017;
    double V = 1.0/rho;
    double specV2 = V*pow(P/Pe,1/k);
    double mdotn=At*P*k*sqrt(numerator/denom);//At*a/V;
    double mdoti=mdotn*OF/(1+OF);
    double v2 = sqrt(2*k*R*T*(1-pow(Pe/P,(k-1)/k))/(k-1));
    cout<<"optimal v2: "<<v2<<" m/s";
    cout<<"optimal A2: "<< mdotn*specV2/(v2*At)<<" m^2";
    cout<<"throat: "<<At<<" m^2";

	//cout<<"mass flow rate out of rocket: "<<mdotn<<" kg/s\n";
	//cout<<"mass flow rate out of injector: "<<mdoti<<" kg/s\n";

	//cout<<"\ntotal distance (feet) "<<yt<<;
	
//	system("PAUSE");
	//cin>>Mt;//so the code doesn't end (couldnt get system("pause") to work
	return 0;
}

