// FCComp Test File
//Compile with:
//	g++ -o FCComptry FCComptry.cpp ../FCComps/FCComp.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../FCComp.h"
#include "../genlib.h"
//#include "../FCComps/MassStream.h"
//#include "../FCComps/isoname.h"
#include <map>
#include <iostream>
#include <string>
#include <set>

using namespace std;

void print( string s)
{
	cout << s << "\n";
}

void printStringSet(set<string> ss)
{
	for (set<string>::iterator i = ss.begin(); i != ss.end(); i++)
	{
		print("\t" + *i);
	}
}

int main()
{
	print("Testing Empty Fuel Cycle Component:");
	FCComp emptyFCC;
	print("\tEmpty Name: " + emptyFCC.name);
	print("\tEmpty Figure Directory: " + emptyFCC.figdir);
	print("\tEmpty Pass Number: " + to_str(emptyFCC.PassNum) );
	print("");

	print("Testing Isotope Initialized Fuel Cycle Component:");
	set<int> isoset ((int []) {922350}, (int []) {922350}+1);
	FCComp isoFCC(isoset, "FCC_Iso");
	print("\tIso Name: " + isoFCC.name);
	print("\tIso Figure Directory: " + isoFCC.figdir);
	print("\tIso Pass Number: " + to_str(isoFCC.PassNum) );
	print("");

	print("Testing Isotope & Paramter Initialized Fuel Cycle Component:");
	string parr [] = {"mass"};
	set<string> parset (parr, parr+length_array(parr));
	FCComp ipFCC(isoset, parset, "FCC_IsoPar");
	print("\tIsotope & Paramter Name: " + ipFCC.name);
	print("\tIsotope & Paramter Figure Directory: " + ipFCC.figdir);
	print("\tIsotope & Paramter Pass Number: " + to_str(ipFCC.PassNum) );
	print("");

	print("Testing Figure Adding:");
	ipFCC.addToFigsList(ipFCC.figdir + "922350.png", "iso");
	ipFCC.addToFigsList(ipFCC.figdir + "mass.png", "param");
	ipFCC.addToFigsList(ipFCC.figdir + "other.png", "other");
	print("\t\t\tDone!");
	print("");

	print("Testing setParams:");
	ipFCC.setParams();
	print("\tTesting Input Params:");
	for (ParamDictIter p = ipFCC.ParamsIn.begin(); p !=  ipFCC.ParamsIn.end(); p++)
	{
		print("\t\t" + p->first + "\t" + to_str(p->second));
	}
	print("\tTesting Output Params:");
	for (ParamDictIter p = ipFCC.ParamsOut.begin(); p !=  ipFCC.ParamsOut.end(); p++)
	{
		print("\t\t" + p->first + "\t" + to_str(p->second));
	}
	print("");

	print("Testing File Pass Writing:");
	//First Pass
	ipFCC.writeIsoPass();
	ipFCC.writeParamPass();

	//Second Pass	
	ipFCC.IsosIn.mass = 1.0;
	ipFCC.IsosIn.comp[922350]  = 0.05;
	ipFCC.IsosOut.mass = 1.0;
	ipFCC.IsosOut.comp[922350] = 0.008;
	ipFCC.writeIsoPass();
	ipFCC.ParamsIn["mass"]  = 10.0;
	ipFCC.ParamsOut["mass"] = 5.0;
	ipFCC.writeParamPass();

	//Third Pass, same data as second
	ipFCC.writeout();

	//Also make sure writeout() bypasses the Params when params2track is empty
	isoFCC.writeout();
	print("\t\t\tDone!");

	return 0;
}
