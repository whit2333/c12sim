/// \file utils.h
/// gemc utility functions \n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef utils_H
#define utils_H 1

// Qt headers
//#include <QtSql>
//#include <QtWidgets>

// C++ headers
#include <iostream>
#include <string>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"

// gemc headers
#include "string_utilities.h"
#include "options.h"



// merging two <string, string>
void mergeMaps(map<string, string>& lhs, const map<string, string>& rhs);


// calculate rotation matrix from input string (MYSQL version)
G4RotationMatrix calc_rotation(string, string);

// calculate shift from input string (MYSQL version)
G4ThreeVector    calc_position(string);


// returns a G4Colour from a string
G4Colour gcol(string);


// utility class used as interface between the various factories
// and the gemc classes
// All leading and trailing spaces are removed
class gtable
{
	public:
		gtable(vector<string> d)
		{
			for(unsigned int i=0; i<d.size(); i++)
				data.push_back(TrimSpaces(d[i]));
		}
		gtable(){data.clear();}
		~gtable(){;}
		
		vector<string> data;
		
		void add_data(QVariant input)
		{
			data.push_back(TrimSpaces(qv_tostring(input)));
		}
		
		void add_data(string input)
		{
			data.push_back(TrimSpaces(input));
		}
		
	
	// Overloaded "<<" for gtable class. Dumps infos on screen.
	friend ostream &operator<<(ostream &stream, gtable);                          

};



// returns a map of files (name is key) and their types (value)
map<string, string> getFilesInDirectory(string);



// returns a string with the current time
string timeStamp();

// returns value + best units as a string
string bestValueUnits(double, string);

// convert a vector of int into a vector of double
vector<double> convertVintVdouble(vector<int> input);

#endif





