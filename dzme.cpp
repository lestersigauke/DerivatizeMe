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

//#define DEBUG

#include "dzme.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <iostream>
#include <utility>
#include <sstream>
#include <vector>
#include <cmath>



#include "result.h"

#define HYDROGENLIST 1
#define HYDROGENRING 2
#define HYDROGENRANDOM 3
#define HYDROGENALL 4
#define HYDROGENSYSTEMATIC 5

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
MoleculeIterator SubstituentIterator;

std::vector<OpenBabel::OBMol> substmols;

std::string coremolecule, hydrogenlist, substituentpath, outputdirectory, logfile, substituentlist, errorfile, outputfile,outputformat; 
int hydrogenswitch = HYDROGENALL;;
std::vector<std::string> hydrogenvector,substituentvector;
std::vector<int> hydrogenintvector;

bool filtering;
bool optimize=false;
bool physico=false;
bool force=true;//force substitution at all sites
bool chonly=true;

int moleculenumber=1;
int confsearch;
int numberofrandom;
int maxsubstitution;
int maxtotal;
//
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void makehydrogenssystematic(result aresult)
{
  //aresult.deleteallbondsinmolecule();
  result thisresult=aresult;
  //thisresult.print();
  // loop through Cl, NO2, F, ...
  for(int j=0;j<substmols.size();j++)
  {
    // loop through the substitutablehydrogens
    //std::cout << "***************** DEALING WITH SUBSTITUENT " << j << std::endl;
    for(int s=0;s<1000;s++)
    {
         //std::cout << "************* The atom is " << s << " but is it ORIGH?? " << std::endl;
         
         if(thisresult.getacc(s)==ORIGH)
         {
           OpenBabel::OBMol copyofsubstituent=substmols[j]; //copy of substituent 
           //at this point we have a hydrogen and a substituent to replace it with
           result newresult;     // create new result classes for each
           //std::cout << "We have created new result because we have substituent " << j << " and s is " << s << std::endl;
           newresult=thisresult;   // copy data from thisresult to each of the result classes
           //newresult.print();
           //std::cout << "copied newresult successfully" << std::endl;
           //std::cout << "about to add copy of substituent to newresult at position " << s << std::endl;
           
           newresult.add(copyofsubstituent,s);    // add the substituent to the appropriate place.
           //newresult.setacc(DELEH,s);
           OpenBabel::OBAtom * h = newresult.themolecule.GetAtom(s+1);//GetAtom uses 1 based index
           //newresult.themolecule.DeleteAtom(h);
           //newresult.print();
           //std::cout << "about to delete hydrogents" << std::endl; 
           //newresult.print();
           //newresult.deleteallbondsinmolecule();
           //newresult.themolecule.PerceiveBondOrders();
           //newresult.themolecule.DeleteHydrogens();           
           //newresult.print();
           result finalresult=newresult;
           finalresult.save();
           //std::cout << "have saved" << std::endl;
         }
        //std::cout << " End of loop " << s << std::endl;
    }
  }
}

void fix(result thisresult)
{
  //std::cout << "maxh " << thisresult.maxh << " >= subsupto " << thisresult.substupto << std::endl;
  //std::cout << "nocurrently " << thisresult.nocurrentlysubst << " < maximumsubs " << thisresult.maximumsubs << std::endl;
  //std::cout << "maximumoutput " << thisresult.maximumoutput << " > currentnumber " << thisresult.thecounter->currentnumber() << std::endl;

      if(  (thisresult.maxh>=thisresult.substupto) && 
           (thisresult.nocurrentlysubst<thisresult.maximumsubs) &&
           (thisresult.maximumoutput>thisresult.thecounter->currentnumber()))
      {
        //NB substupto
        OpenBabel::OBMolAtomIter a(thisresult.themolecule);
        while((*a).GetIndex()<thisresult.substupto){++a;};//ignore all hydrogens until substupto
        // at this point (*a) is an atom
        // e.g. if( (*a).IsHydrogen() )    // (*a).GetType()    will return the type as char * // (*a).GetIdx()     will return a number as the index of the atom

        //std::cout << "Atom : " << (*a).GetType() << " has atomic mass : " << (*a).GetAtomicMass() 
        //        << " and Index : " << (*a).GetIndex() << " and Id : " << (*a).GetId() << std::endl;

        // what about the next substupto (*a) is the current hydrogen atom
        // we have to go through the accounting as well to find a hydrogen that has ORIGH
        
        thisresult.setacc(DELEH,thisresult.substupto); // set the substupto as DELEH in each of them
        // *******************************************************
        // now we need to iterate through the substituents
        for(int j=0;j<substmols.size();j++)
        {
           //std::cout << "Substituting now in fix with " << substmols[j].GetFormula() << std::endl;
           
           result newresult;     // create new result classes for each
           newresult=thisresult;   // copy data from thisresult to each of the result classes
           OpenBabel::OBMol copyofsubstituent=substmols[j];    
           newresult.add(copyofsubstituent,thisresult.substupto);    // add the substituent to the appropriate place.
           
           newresult.nocurrentlysubst++;
           int s=newresult.substupto+1;
           while((newresult.getacc(s)!=ORIGH)&&(s<999)){s++;}; // set substupto in each of them to be the next H.
           newresult.substupto=s;
           fix(newresult);
        }
        if(!force)
        {
           //allow unsubstituted to proceed
           result newresult;     // create new result classes for each
           newresult=thisresult;   // copy data from thisresult to each of the result classes
           newresult.setacc(IGNOH,newresult.substupto);
           ///IGNOH
          int s=newresult.substupto+1;
           while((newresult.getacc(s)!=ORIGH)&&(s<999)){s++;}; // set substupto in each of them to be the next H.           
           newresult.substupto=s;           
           fix(newresult);           
        }        
      }
      else
      {
         //all hydrogens have been substituted and so we need to delete hydrogens and save etc.
         //thisresult.themolecule.DeleteHydrogens();
         //std::cout << "maximumoutput " << thisresult.maximumoutput << " > currentnumber " << thisresult.thecounter->currentnumber() << std::endl;
         if(thisresult.maximumoutput>thisresult.thecounter->currentnumber())
         {           
           //thisresult.print();
           thisresult.save();          
         }
      }


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
}

// Everything starts here
int main(int argc, char** argv)
{
  hydrogenswitch=HYDROGENRANDOM;
  numberofrandom=30;
  substituentpath="./subs";
  outputformat="none";
  logfile="logfile.txt";
  outputfile="output.sdf";
  outputdirectory="./output";
  substituentlist="cl.pdb,nh2.pdb";
  coremolecule="core.pdb";
  maxsubstitution=100;
  maxtotal=100000;
  optimize=false;
  chonly=true;
  confsearch=0;
  filtering=false;

//std::cout << "argc is " << argc << std::endl;
//////////////////////////////////////////////////////////////////////////
// this is all about the command line inputs
 if(argc>1)
 {
  std::string firststring=std::string(argv[1]);
  if (firststring.compare("-help")==0)
  {
     std::cout << "DZME help" << std::endl << std::endl;
     std::cout << "Correct usage:" << std::endl << std::endl;
     std::cout << "           ./dzme -m coremolecule -h hydrogensubs -I substituentdir -O outputdir -Of outputformat -L logfile -s listofsubstituents -x sdfoutput -m maximumsubs" << std::endl << std::endl;
     std::cout << "Example:   ./dzme -m ./molecule.xyz -h list 7,9,13 -I ./substituents -O ./outputdir -L ./logfile.txt" << std::endl;
     std::cout << "            -s cl.xyz,no2.xyz,ch3.xyz -x output.sdf -i " << std::endl << std::endl;
     std::cout << " -m coremolecule:" << std::endl;
     std::cout << "                    core molecule to be derivatized (default \"core.pdb\")" << std::endl;
     std::cout << " -h hydrogensubs:" << std::endl;
     std::cout << "    hydrogensubs:   list 7,9,13     {only allows substitution of hydrogens 7,9,13 systematic}" << std::endl;
     std::cout << "                                    {note that this must be the exact atom numbers in the entire molecule}" << std::endl;
     std::cout << "                    all             {allows substitution of all hydrogens in a systematic manner}" << std::endl;
     std::cout << "                    ring            {allows substitution of all hydrogens attached to rings a systematic manner}" << std::endl;
     std::cout << "                    random  30      {allows substitution of approximately thirty percent of all hydrogens systematically (default)" << std::endl;
     std::cout << "                    systematic      {allows monosubstitution of all hydrogens in a systematic manner}" << std::endl;
     std::cout << " -I substituentdir:" << std::endl;
     std::cout << "                    directory containing substituents (default \"./subs\")" << std::endl;
     std::cout << " -O outputdir:" << std::endl;
     std::cout << "                    directory containing individual files created (default \"./output\")" << std::endl;
     std::cout << " -Of outputformat:" << std::endl;
     std::cout << "                    file type for individual files created (default \"none\")" << std::endl;
     std::cout << "                    set as \"none\" to only create an sdf file output" << std::endl;
     std::cout << " -L logfile:" << std::endl;
     std::cout << "                    log file containing smiles/energies of all structures" << std::endl;
     std::cout << "                    default \"logfile.txt\"" << std::endl;
     std::cout << " -s substituents:" << std::endl;
     std::cout << "                    list of structures in the substituent directory that will be used for functionalization" << std::endl;
     std::cout << "                    default is \"cl.pdb,nh2.pdb\"" << std::endl;
     std::cout << " -x outputsdf:" << std::endl;
     std::cout << "                    sdf file containing all generated structures (default \"output.sdf\")" << std::endl;
     std::cout << " -i:" << std::endl;
     std::cout << "                    for all and ring, only make C-H's substitutable (default)" << std::endl;
     std::cout << "                    ignored for -h list" << std::endl;
     std::cout << " -a:" << std::endl;
     std::cout << "                    for all and ring, make all H's substitutable" << std::endl;
     std::cout << "                    ignored for -h list" << std::endl;
     std::cout << " -maxs maxsubstituents:" << std::endl;
     std::cout << "                    the maximum number of substituents in the products (default 100)" << std::endl;
     std::cout << " -maxt maxtotal:" << std::endl;
     std::cout << "                    cut off production when maxtotal is reached (default 100000)" << std::endl;
     std::cout << " -opt:" << std::endl;
     std::cout << "                    perform MMFF optimization(default off)" << std::endl;
     std::cout << " -phys:" << std::endl;
     std::cout << "                    calculate physicochemical properties(default off)" << std::endl;
     std::cout << " -noforce:" << std::endl;
     std::cout << "                    allow sites to remain unsubstituted(default force unless maxs set)" << std::endl;
     std::cout << " -conf:" << std::endl;
     std::cout << "                    search for a better conformer(default off)" << std::endl;

     exit(0);
 }
 else
 {
   //std::cout << argv[0] << std::endl;
   for (int i = 1; i < argc; i++)
      {
       //if (i + 1 != argc)
        {   
            //std::cout << "dealing with " << argv[i] << std::endl;
            std::string argstring=std::string(argv[i]);
            if (argstring.compare("-m")==0)
            {
              if(i+1 == argc){std::cout << "usage -m coremolecule" << std::endl;exit(1);}
              coremolecule = argv[i + 1]; std::cout << "found -m " << coremolecule << std::endl;
            }
            else if (argstring.compare("-I")==0)
            {
              if(i+1 == argc){std::cout << "usage -I substituentpath" << std::endl;exit(1);}
              substituentpath = argv[i+1]; std::cout << "found -I "<< substituentpath << std::endl; 
            }
            else if (argstring.compare("-O")==0)
            {
              if(i+1 == argc){std::cout << "usage -O outputdirectory" << std::endl;exit(1);}
              outputdirectory = argv[i+1]; std::cout << "found -O "<< outputdirectory << std::endl; 
            }
            else if (argstring.compare("-L")==0)
            {
              if(i+1 == argc){std::cout << "usage -L logfile.txt" << std::endl;exit(1);}
              logfile = argv[i+1];  std::cout << "found -L "<< logfile << std::endl; 
            }
            //else if (argstring.compare("-E")==0)
            //{
		        //  errorfile = argv[i+1]; std::cout << "found -E " << errorfile << std::endl;
            //}
            else if (argstring.compare("-Of")==0)
            {
              if(i+1 == argc){std::cout << "usage -Of pdb or -Of xyz etc" << std::endl;exit(1);}
		          outputformat = argv[i + 1]; //std::cout << "found -Of " << outputformat << std::endl;
            }                
            else if (argstring.compare("-x")==0)
            {
              if(i+1 == argc){std::cout << "usage -x sdfoutputfilename" << std::endl;exit(1);}
              outputfile = argv[i + 1]; //std::cout << "found -x " << outputfile << std::endl;
            }
            else if (argstring.compare("-s")==0)
            {
              if(i+1 == argc){std::cout << "usage -s listofsubst e.g -s cl.pdb,nh2.xyz" << std::endl;exit(1);}
              substituentlist = argv[i+1];  //std::cout << "found -s " << substituentlist << std::endl;
            }
            else if (argstring.compare("-maxs")==0)
            {
              if(i+1 == argc){std::cout << "usage -maxs maximumsubstitutions e.g -maxs 3" << std::endl;exit(1);}
              maxsubstitution = atoi(argv[i+1]);  //std::cout << "found -maxs " << substituentlist << std::endl;
              force=false;
            }            
            else if (argstring.compare("-maxt")==0)
            {
              if(i+1 == argc){std::cout << "usage -maxt maximumderivatives e.g -maxt 10" << std::endl;exit(1);}
              maxtotal = atoi(argv[i+1]);  //std::cout << "found -maxt " << substituentlist << std::endl;
            }            
            else if (argstring.compare("-conf")==0)
            {
              if(i+1 == argc){std::cout << "usage -conf level e.g -conf 0" << std::endl;exit(1);}
              confsearch = atoi(argv[i+1]);
              optimize=true;
            }             
            //else if (argstring.compare("-f")==0)
            //{
            //  filtering=true; //UNUSED               
            //}
           else if (argstring.compare("-opt")==0)
            {
                optimize=true; //std::cout << "found -opt " << std::endl;
            }             
           else if (argstring.compare("-phys")==0)
            {
                physico=true; //std::cout << "found -phys " << std::endl;
            }             
            else if (argstring.compare("-noforce")==0)
            {
                force=false; //std::cout << "found -noforce " << std::endl;
            }             
            else if (argstring.compare("-i")==0)
            {
                chonly=true; //std::cout << "found -i " << std::endl;
            }             
            else if (argstring.compare("-a")==0)
            {
                chonly=false; //std::cout << "found -a " << std::endl;
            }             
            else if (argstring.compare("-h")==0)
            {
              if(i+1 == argc){std::cout << "usage -h all/ring/random/list" << std::endl;exit(1);}
              std::string hydrogenswitchstring=std::string(argv[i+1]); //std::cout << "found -h " << hydrogenswitchstring << std::endl;
              if (hydrogenswitchstring.compare("all")==0) // substitute all hydrogens
              {
                hydrogenswitch = HYDROGENALL; // turns on the use of all hydrogen substitution
              }
              else if (hydrogenswitchstring.compare("ring")==0)  //substitute only ring hydrogens
              {
       			    hydrogenswitch = HYDROGENRING; //turns on the switch for ring substitutions only
			        }
              else if (hydrogenswitchstring.compare("random")==0)
              {
       			    hydrogenswitch = HYDROGENRANDOM; //turns on qthe switch for random substitutions
                if(i+2 == argc){std::cout << "usage -h random percent e.g. -h random 30" << std::endl;exit(1);} 
                numberofrandom = atoi(argv[i + 2]);
			        }
              else if (hydrogenswitchstring.compare("systematic")==0)
              {
       			    hydrogenswitch = HYDROGENSYSTEMATIC; //non recursive
			        }
              // Check whether hydrogenlist == a (all); r(rings); or 1,2,3,4 (specific hydrogens)
              else if (hydrogenswitchstring.compare("list")==0) // use of some regular expression to select for a list e.g. commainhydrogenlist
              {
                hydrogenswitch = HYDROGENLIST; // turns on the use of specific hydrogens list
                //in here we can split and get a list of hydrogens as strings --
                if(i+2 == argc){std::cout << "usage -h list listofh e.g. -h list 7,9,13" << std::endl;exit(1);} 
                hydrogenlist = argv[i + 2];
                split(hydrogenvector,hydrogenlist,",");
                for(int q=0;q<hydrogenvector.size();q++)
                {
                  int t;
                  std::istringstream (hydrogenvector[q]) >> t;
                  //std::cout << "hydrogenlist member no: " << q << " is "<< hydrogenvector[q] << " and " << t << std::endl; 
                  hydrogenintvector.push_back(t);
                }
                //at this point we can push_back into vector the hydrogens as integers --
                for(int q=0;q<hydrogenintvector.size();q++)
                {
                  //std::cout << "hydrogenintlist member no: " << q << " is "<< hydrogenintvector[q] << std::endl; 
                }
              }
            }                
            else
            { 
              //std::cout << "Please try again." << std::endl;
              //sleep(2000);
              //exit(0);
            }
           
         }
       }
  }  
 }

 //split works nicely in this way, we send substituentlist ("no.xyz,ch3.xyz,oh.xyz")
 //and we say the delimiter is "," then it sends back to us
 //substituentvector[0]=no.xyz, substituentvector[1]=ch3.xyz, substituentvector[2]=oh.xyz
 //and substituentvector.size()==3

 split(substituentvector,substituentlist,",");
 for(int j=0;j<substituentvector.size();j++)
 {
  std::string fullname;
  fullname = substituentpath + "/" + substituentvector[j];
  //std::cout << "substituent number " << j << " is " << substituentvector[j] << " ::: with path " << fullname <<std::endl;
 }

 std::cout << "---------------------------------------------------" << std::endl;                   
 std::cout << "DerivatizeME" << std::endl;                   
 std::cout << "Lester Sigauke, Kevin A. Lobb" << std::endl;                   
 std::cout << "2020, CMCDD Research Group" << std::endl;                   
 std::cout << "Chemistry Dept, Rhodes University" << std::endl;                   
 std::cout << "---------------------------------------------------" << std::endl;                   
 std::cout << "core molecule is --                -m " << coremolecule << std::endl;
 std::cout << "hydrogen list (if used) is --      -h list " << hydrogenlist << std::endl;
 std::cout << "% random hydrogens (if used) is -- -h random " << numberofrandom << std::endl;
 std::cout << "substituents are located --        -I  " << substituentpath << std::endl;
 std::cout << "the good stuff gets dumped at --   -O  " << outputdirectory << std::endl;
 std::cout << "   with a format of --             -Of " << outputformat << std::endl;
 std::cout << "the logfile is --                  -L  " << logfile << std::endl;
 //std::cout << "the errorfile is --" << errorfile << std::endl;
 std::cout << "the outputfile is --               -x  " << outputfile << std::endl;
 std::cout << "the substituents are --            -s  " << substituentlist << std::endl;
 std::cout << "maximum substitution level is --   -maxs " << maxsubstitution << std::endl;
 std::cout << "maximum number of derivatives is   -maxt " << maxtotal << std::endl;
 std::cout << "optimize is set --                 -opt " << optimize << std::endl;
 std::cout << "conformer search is set --         -conf " << confsearch << std::endl;
 std::cout << "physico is set --                  -phys " << physico << std::endl << std::endl;
 std::cout << "**** for help type ./dzme -help" << std::endl;
 std::cout << "---------------------------------------------------" << std::endl;
 //exit(0);
 //return 0;


 OpenBabel::OBConversion obconversion;
 OpenBabel::OBMol natprod;
 

 ////////////////////////////////////////////////////
 // read in the substituents into a substituent list


 // Substituents now as strings
 for(int j=0;j<substituentvector.size();j++)
 {
   std::string fullname,filetype;
   fullname = substituentpath + "/" + substituentvector[j];
   filetype = fullname.substr(fullname.find_last_of(".")+1);
   //std::cout << "file extension for subst " << fullname << " is " << filetype << std::endl;
   obconversion.SetInFormat(filetype.c_str());
   OpenBabel::OBMol newmolecule;
   obconversion.ReadFile(&newmolecule,fullname);
   //substituents going into substmols.
   //NB This is where substmols gets populated with all the substituents as openbabel mols
   // ************************************************************************************
   substmols.push_back(newmolecule);
   //std::cout << "substituent number " << j << " is " << substituentvector[j] 
   //          << " ::: with path " << fullname <<std::endl;
 }

 // at this point we have *******substmols********* (in a list)

 std::string filetype;
 filetype = coremolecule.substr(coremolecule.find_last_of(".")+1);
 //std::cout << "file extension for coremolecule " << coremolecule << " is " << filetype << std::endl;
 obconversion.SetInFormat(filetype.c_str());
 //obconversion.SetInFormat("pdb");
 //Natural product as pdb
 //NB read in the natural product
 obconversion.ReadFile(&natprod, coremolecule);
    
 std::cout << "coremolecule has formula " << natprod.GetFormula() << std::endl;
 std::cout << "---------------------------------------------------" << std::endl;
 std::cout << "Starting Generation - Please be patient: " << std::endl;

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //NB create the initial result object and two objects, one for counting and the other for logfile
 counter mycounter;             //create counter object
 logfilemaker mylogfilemaker;          //create logfile object
 errorfilemaker myerrorfilemaker;      //create an errorfile object
 multiplestructuremaker mystructuremaker;
   
  //name of xxxxfile stored in xxxxfile
 mylogfilemaker.openlogfile(logfile);
 myerrorfilemaker.openerrorfile(errorfile); 
 mystructuremaker.openfile(outputfile);  

 result firstnode;              //create result object called firstnode
 //firstnode.theoutputdirectory=outputdirectory;
 firstnode.setoutputdirectory(outputdirectory);
 firstnode.setoutputformat(outputformat);
 firstnode.setmaximumsubs(maxsubstitution);
 firstnode.setmaximumoutput(maxtotal);
 firstnode.setoptimize(optimize);
 firstnode.setphysico(physico);
 firstnode.setforce(force);
 firstnode.setchonly(chonly);
 firstnode.setconfsearch(confsearch);
 firstnode.resetall(natprod);     //we are going to have to in resetall sort out accouting etc
  
 //make connections between objects
 firstnode.connecttocounter(&mycounter);      // make connection between result object firstnode and 
 firstnode.connecttologfile(&mylogfilemaker);     // logfilemaker object mylogfilemaker
 firstnode.connecttoerrorfile(&myerrorfilemaker);    // errorfilemaker object myerrorfilemaker
 firstnode.connecttomaker(&mystructuremaker);
                                            
  //////////////////////////////////////////////
  //THESE ARE OUR CHOICES::
  //1. Go through the hydrogen list as integers and set  firstnode.substitutablehydrogen(int i)
  //for(int q=0;q<hydrogenintvector.size();q++)
  //{
  //  firstnode.substitutablehydrogen(q);
  //}
  //or
  //2. firstnode.makeallhydrogenssubstitutable()
  //or
  //3. firstnode.makerandomhsubstitutable(fraction)
  //or
  //4. 
  ////////////////////////////////////////////////
  //HYDROGENLIST HYDROGENRING HYDROGENRANDOM HYDROGENALL HYDROGENSYSTEMATIC

  if (hydrogenswitch == HYDROGENRING)
  {firstnode.makeonlyringhydrogenssubstitutable();}

  else if (hydrogenswitch == HYDROGENALL)
  {firstnode.makeallhydrogenssubstitutable();}

  else if (hydrogenswitch == HYDROGENRANDOM)
  {
   firstnode.makerandomhydrogenssubstitutable(numberofrandom);
  }
  else if (hydrogenswitch == HYDROGENSYSTEMATIC)
  {
   firstnode.makeallhydrogenssubstitutable();
  }
  else if (hydrogenswitch == HYDROGENLIST)
  {
     for(int q=0;q<hydrogenintvector.size();q++)
     	{

    	    firstnode.substitutablehydrogen(hydrogenintvector[q]);
  	}
  }
  ////////////////////////////////////////////
  //just to check the above has worked.  CONFIRMS MOLECULAR MASS
  //firstnode.printacc();
  //firstnode.setpredlog();

  firstnode.set_substupto_and_maxh();
  
  //resultinitial2.add(natprod,0,2);//we need to add openbabel moleculed ot the object
  //resultinitial2.print();//????

  ///////////////////////////////////////////////////
  // real work is done here, first is nested for, second is recursive:
  
  if (hydrogenswitch == HYDROGENSYSTEMATIC)
    makehydrogenssystematic(firstnode);
  else 
    fix(firstnode);//do the work in fully recursive mode
  
  // Remember to close the file that you opened
  mylogfilemaker.closelogfile();
  myerrorfilemaker.closeerrorfile();
  mystructuremaker.closefile();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
#ifdef DEBUG
  std::cout << "out main" << std::endl;
#endif
*/
std::cout << std::endl;
std::cout << "Complete." << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
