#ifndef CMCDD_RESULT_H
#define CMCDD_RESULT_H
//-----------------------------------------------------------------------------
// Copyright (c) 2020 CMCDD Research Group, Lester Sigauke, Kevin Lobb
//                    Chemistry Dept, Rhodes University
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//-----------------------------------------------------------------------------

#include <openbabel/mol.h>
#include <mathlib.h>

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
#define ORIGH -2   //we want these substituted
#define IGNOH 0    //we want these H's to remain
#define DELEH -1   //we want these deleted because they are overlapping the bond!
#define NOTAH 1


//local manipulation tools
void Translate(OpenBabel::OBMol * thesubst, Vector3 from, Vector3 to);
void Rotate(OpenBabel::OBMol * thesubst, Vector3 mainvector_norm, Vector3 substvector_norm);

void BondRotate(OpenBabel::OBMol * thesubst, Vector3 axis,float angle);
double AtomDistance(OpenBabel::OBAtom * atoma,OpenBabel::OBAtom * atomb);
int badcontacts(OpenBabel::OBMol * mol1,OpenBabel::OBMol * mol2);


class counter
{
  public:
     counter()
     {
        number=0;
     };
     ~counter(){};
     int requestnumber()
     {
        number++;
        return number;
     };
     int currentnumber(){return number;}

  private:
     int number;
};

//define the class and the functions associated with the class
class multiplestructuremaker
{
  public:
    multiplestructuremaker(){};
    ~multiplestructuremaker(){};

   void openfile(std::string filename);
   void closefile();
   void writetofile(std::string stringtowrite);
   private:
   std::ofstream myfile;

};



class errorfilemaker
{ 
   public:
    errorfilemaker(){};
    ~errorfilemaker(){};
    
    void openerrorfile(std::string filename);
    void closeerrorfile();
    void writetoerrorfile(std::string stringtowrite);
    
    private:
     std::ofstream myefile;
};


class logfilemaker
{
   public:
    logfilemaker(){};
    ~logfilemaker(){};

    void openlogfile(std::string filename);
    void closelogfile();
    void writetologfile(std::string stringtowrite);
    void writetologfile(int inttowrite);
    void writetologfile(float floattowrite);

   private:
    std::ofstream myfile;
};


class result
{
//this class contains a molecule and an array telling about the Hydrogens
public:
   result() //this is the constructor, when created this called automatically
   { 
      theoutputformat=std::string("xyz") ;
      theoutputdirectory=std::string("./");   
      optimize=false;
      nocurrentlysubst=0;
      force=true;
      chonly=true;
      confsearch=0;
      HBA=0;HBD=0;LOGP=0.0;MW=0.0;TPSA=0.0;      
   };
   ~result(){themolecule.Clear();};//this is called when object is destroyed

   std::string theoutputdirectory;
   std::string theoutputformat;
   std::string formula;

   OpenBabel::OBMol themolecule;
   int accounting[1000];
   int substupto;//first h to be substituted
   int maxh;//last h to be substituted
   int nocurrentlysubst;//current number of substitutions that have taken place


   //Details for calculation
   int maximumsubs;//maximum allowed number of substitutions
   int maximumoutput;//maximum allowed number of outputs
   
   int confsearch;//level of conformational search to be performed on final structure
   bool optimize;//MMFF optimization?   
   bool physico;//calculate physicochemical props?

   bool force;//force every marked hydrogen to be substituted at the end - allow no unsubstituted product
   bool chonly;//only allow H's bonded to C to be substitutable (ignored for h list)

   //// Descriptors
   int HBD; //no of Hydrogen Bond Donors
   int HBA; //no of Hydrogen Bond Acceptors    
   float LOGP; // octanol/water partition coefficient
   float MW;
   float TPSA;

   //float approxmass;

   // make an instance of the variable counter, errofilemaker, logfilemaker etc.
   counter * thecounter;
   errorfilemaker * theerrorfilemaker; 
   logfilemaker * thelogfilemaker;
   multiplestructuremaker * thestructures;
 
   // connect to objects
  
   void connecttomaker(multiplestructuremaker * newmaker){thestructures=newmaker;}
   void connecttocounter(counter * newcounter){thecounter=newcounter;};
   void connecttoerrorfile(errorfilemaker* newerrorfile){theerrorfilemaker=newerrorfile;};
   void connecttologfile(logfilemaker* newlogfile){thelogfilemaker=newlogfile;};

   void printacc();
   
   //set functions
   void setoutputdirectory(std::string newdir){theoutputdirectory=newdir;}
   void setoutputformat(std::string newformat){theoutputformat=newformat;}
   void setoptimize(bool opt){optimize=opt;}
   void setphysico(bool phys){physico=phys;}
   void setforce(bool f){force=f;}
   void setchonly(bool c){chonly=c;}
   void setconfsearch(int c){confsearch=c;}
   void setmaximumsubs(int max){maximumsubs=max;}
   void setmaximumoutput(int max){maximumoutput=max;}

   //types in hydrogen array (ignore replace delete etc)
   void setacc(int value,int index){accounting[index]=value;}
   void removeacc(int index){for(int i=index;i<999;i++){accounting[i]=accounting[i+1];};}

   void substitutablehydrogen(int i); //make particular hydrogen have attribute ORIGH
   void makeallhydrogenssubstitutable();//make all hydrogens have attribute ORIGH
   void makeonlyringhydrogenssubstitutable();//make only ring hydrogens substitutable as ORIGH
   void makerandomhydrogenssubstitutable(int numberofrandom);//only 8
   void intelligentlymakehydrogenssubstitutable();
   void makehydrogenssystematic();
   void set_substupto_and_maxh();//substupto will start off at the first hydrogen to substitute
                                 //if substupto == maxh we are at the last hydrogen to substitute
   //accessor
   int getacc(int index){return accounting[index];}
   
   
   void resetall(OpenBabel::OBMol moleculewewantinobject); //by default all hydrogens IGNOH
   void ignorehydrogen(int i); //atom at this number must not be substituted
                               //make particular hydrogen have attribute IGNOH

   //void deleteallbondsinmolecule();
   void deletedeletablehydrogens();

   int count(int level);
 
   //modifiers
   result & operator= (result &other);
   void add(OpenBabel::OBMol &other, int hydrogenatomtodelete);

   //print and save
   void print();
   void printatoms();
   void printbonds();
   void save();
   //// Descriptors
   void calculatedescriptors();
   void calculateformula();
};

#endif