// Qt headers
//#include <QApplication>
//#include <QSplashScreen>
//#include <QtWidgets>


// gemc headers
#include "utils.h"
#include "string_utilities.h"
//#include "splash.h"

// G4 headers
#include "G4UnitsTable.hh"


// C++ headers
#include "dirent.h"



// merging two <string, string>
// The rhs overwrites what's in the lhs
void mergeMaps(map<string, string>& lhs, const map<string, string>& rhs)
{
	for(map<string, string>::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	{
		lhs[it->first] = it->second;
	}
}


// calculate rotation matrix from input string (MYSQL version)
G4RotationMatrix calc_rotation(string r, string dname)
{
	G4RotationMatrix rot(G4ThreeVector(1, 0, 0),
						 G4ThreeVector(0, 1, 0),
						 G4ThreeVector(0, 0, 1));
	
	stringstream vars(r);
	string var;
	vars >> var;
	
	//	cout << dname
	
	if(var != "ordered:")
	{
		             rot.rotateX(get_number(var,1));
		vars >> var; rot.rotateY(get_number(var,1));
		vars >> var; rot.rotateZ(get_number(var,1));
	}
	else
	{
		string order;
		vars >> order;
		if(order == "xzy")
		{
			vars >> var; rot.rotateX(get_number(var,1));
			vars >> var; rot.rotateZ(get_number(var,1));
			vars >> var; rot.rotateY(get_number(var,1));
		}
		else if(order == "yxz")
		{
			vars >> var; rot.rotateY(get_number(var,1));
			vars >> var; rot.rotateX(get_number(var,1));
			vars >> var; rot.rotateZ(get_number(var,1));
		}
		else if(order == "yzx")
		{
			vars >> var; rot.rotateY(get_number(var,1));
			vars >> var; rot.rotateZ(get_number(var,1));
			vars >> var; rot.rotateX(get_number(var,1));
		}
		else if(order == "zxy")
		{
			vars >> var; rot.rotateZ(get_number(var,1));
			vars >> var; rot.rotateX(get_number(var,1));
			vars >> var; rot.rotateY(get_number(var,1));
		}
		else if(order == "zyx")
		{
			vars >> var; rot.rotateZ(get_number(var,1));
			vars >> var; rot.rotateY(get_number(var,1));
			vars >> var; rot.rotateX(get_number(var,1));
		}
		else
		{
			cout << "     >> ERROR: Ordered rotation <" << order << "> for " << dname << " is wrong, it's none of the following:"
			<< " xzy, yxz, yzx, zxy or zyx. Exiting." << endl;
			exit(0);
		}
	}
	
	return rot;
	
}


// calculate shift from input string (MYSQL version)
G4ThreeVector calc_position(string v)
{
	
	G4ThreeVector pos(0, 0, 0);
	stringstream vars(v);
	string var;
	
	vars >> var; pos.setX(get_number(var,1));
	vars >> var; pos.setY(get_number(var,1));
	vars >> var; pos.setZ(get_number(var,1));
	
	return pos;
}


// returns a G4Colour from a string
G4Colour gcol(string cvar)
{
	G4Colour thisCol;
	
	// if color is 6 digits then it's only rrggbb. Setting transparency to zero
	if(cvar.size() == 6)
		thisCol = G4Colour(strtol(cvar.substr(0, 2).c_str(), NULL, 16)/255.0,
						   strtol(cvar.substr(2, 2).c_str(), NULL, 16)/255.0,
						   strtol(cvar.substr(4, 2).c_str(), NULL, 16)/255.0,
						   1);
	
	// Transparency 0 to 5 where 5=max transparency  (default is 0 if nothing is specified)
	else if(cvar.size() == 7)
		thisCol = G4Colour(strtol(cvar.substr(0, 2).c_str(), NULL, 16)/255.0,
						   strtol(cvar.substr(2, 2).c_str(), NULL, 16)/255.0,
						   strtol(cvar.substr(4, 2).c_str(), NULL, 16)/255.0,
						   1.0 - stringToDouble(cvar.substr(6, 1))/5.0);
	
	return thisCol;
}




map<string, string> getFilesInDirectory(string directory)
{
	map<string, string> filesMap;
	
	DIR *dir;
	struct dirent *ent;
	dir = opendir(directory.c_str());
	if (dir != NULL)
	{
		int len;
		while ((ent = readdir (dir)) != NULL)
		{
			len = strlen(ent->d_name);
			
			// checking various extensions
			if(strcmp(".dat", &(ent->d_name[len - 4])) == 0)
				filesMap[directory + "/" + ent->d_name] = "ASCII" ;
			
			if(strcmp(".txt", &(ent->d_name[len - 4])) == 0)
				filesMap[directory + "/" + ent->d_name] = "ASCII" ;
			
		}
		closedir (dir);
	}
	else
	{
		cout << "    Error: directory " << directory << " could not be opened." << endl;
	}
	
	return filesMap;
}
#include "ctime"
string timeStamp()
{
	
	time_t now = time(NULL);
	struct tm * ptm = localtime(&now);
	char buffer[32];
	// Format: Mo, 15.06.2009 20:20:00
	strftime (buffer, 32, "%a, %m.%d.%Y %H:%M:%S", ptm);
	
	return string(buffer);
}



// returns value + best units as a string
string bestValueUnits(double value, string unitType)
{
	stringstream bestVU;
	bestVU << G4BestUnit(value, unitType);
	string var, uni;
	bestVU >> var >> uni;
	
	return var + "*" + uni;
	
}



ostream &operator<<(ostream &stream, gtable gt)
{
	cout  << endl;
	
	for(unsigned i=0; i<gt.data.size(); i++)
		cout << "   data # " << i << ": " << gt.data[i] << endl;
	
	return stream;
}







vector<double> convertVintVdouble(vector<int> input)
{
	vector<double> out;
	
	for(unsigned i=0; i<input.size(); i++)
		out.push_back(input[i]);
	
	return out;
}








