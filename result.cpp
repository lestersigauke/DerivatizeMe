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

#include "result.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/forcefield.h>
#include <openbabel/conformersearch.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local Manipulation Tools

void Translate(OpenBabel::OBMol * thesubst, Vector3 from, Vector3 to)
{
    //subst_h will be translated to a point on the CH bond to create a longer bond for the new CC^M
     Vector3 displacement = to-from;

     for(OpenBabel::OBMolAtomIter a(thesubst);a;a++)
     {
            float x=a->GetX()+displacement.x;
            float y=a->GetY()+displacement.y;
            float z=a->GetZ()+displacement.z;
            a->SetVector(x,y,z);           
     }

}

void Rotate(OpenBabel::OBMol * thesubst, Vector3 mainvector_norm, Vector3 substvector_norm)
{
   //vectors in opposite directions for main the CH vector, for subst the HC vector
   Vector3 main=mainvector_norm;
   Vector3 subst=substvector_norm;
   main.normalize();
   subst.normalize();
   Matrix3 mymatrix;
   
    Vector3 axis = Vector3::cross(main,subst);
        axis.normalize();
        float angle = Math::radiansToDegrees(std::acos(Vector3::dot(main,subst)));
   
   mymatrix.rotate(axis,180.0-angle);
     
   for(OpenBabel::OBMolAtomIter a(thesubst);a;a++)
     {
        Vector3 atomvector=Vector3(a->GetX(),a->GetY(),a->GetZ());
        Vector3 finalvect=atomvector*mymatrix;
        a->SetVector(finalvect.x,finalvect.y,finalvect.z);      
     }
}



void BondRotate(OpenBabel::OBMol * thesubst, Vector3 axis,float angle)
{

   Matrix3 mymatrix;
   OpenBabel::OBAtom * hs = thesubst->GetAtom(1);
   Vector3 from = Vector3(hs->GetX(),hs->GetY(),hs->GetZ());
   Vector3 to = Vector3(0.000,0.000,0.000);
   //move the substituent to the origin before rotation
   Translate(thesubst,from,to);

   axis.normalize();
   mymatrix.rotate(axis,180.0-angle);
     
   for(OpenBabel::OBMolAtomIter a(thesubst);a;a++)
     {
        Vector3 atomvector=Vector3(a->GetX(),a->GetY(),a->GetZ());
        Vector3 finalvect=atomvector*mymatrix;
        a->SetVector(finalvect.x,finalvect.y,finalvect.z);      
     }
   Translate(thesubst,to,from);

}

double AtomDistance(OpenBabel::OBAtom * atoma,OpenBabel::OBAtom * atomb)
{
   double deltax=atoma->GetX()-atomb->GetX();
	 double deltay=atoma->GetY()-atomb->GetY();
	 double deltaz=atoma->GetZ()-atomb->GetZ();
	 deltax=deltax*deltax;
	 deltay=deltay*deltay;
	 deltaz=deltaz*deltaz;
	 double total=deltax+deltay+deltaz;
	 total=std::sqrt(total);   
	 return total;
}

int badcontacts(OpenBabel::OBMol * mol1,OpenBabel::OBMol * mol2)
{
  int count=0;
  for(OpenBabel::OBMolAtomIter a(mol1);a;a++)
     {
     for(OpenBabel::OBMolAtomIter b(mol2);b;b++)
       {
           if(AtomDistance((&*a),(&*b))<1.6)count++;
       }
     }
  return count;
}

int badcontacts(OpenBabel::OBMol * mol1)
{
  int count=0;
  for(OpenBabel::OBMolAtomIter a(mol1);a;a++)
     {
     for(OpenBabel::OBMolAtomIter b(mol1);b;b++)
       {
           if((AtomDistance((&*a),(&*b))<1.1))count++;//&&(a!=b)
       }
     }
  int atoms=mol1->NumAtoms();
  count = count-atoms;
  return count;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// File Tools

void multiplestructuremaker::openfile(std::string filename)
{
    const char * filechars = filename.c_str();
    myfile.open(filechars);
}

void errorfilemaker::openerrorfile(std::string filename)
{   
    const char * filechars = filename.c_str();
    myefile.open(filechars);
}

void logfilemaker::openlogfile(std::string filename)
{
   const char * filenamechars = filename.c_str();
   myfile.open(filenamechars);
}

//close the files 
void multiplestructuremaker::closefile()
{
    myfile.close();
}

void errorfilemaker::closeerrorfile()
{
   myefile.close();
}

void logfilemaker::closelogfile()
{
   myfile.close();
}

//write to the files

void errorfilemaker::writetoerrorfile(std::string stringtowrite)
{
  myefile << stringtowrite << std::flush;
}

void logfilemaker::writetologfile(std::string stringtowrite)
{
  myfile << stringtowrite << std::flush;
}
void logfilemaker::writetologfile(int inttowrite)
{
  myfile << inttowrite << std::flush;
}
void logfilemaker::writetologfile(float floattowrite)
{
  myfile << floattowrite << std::flush;
}

void multiplestructuremaker::writetofile(std::string stringtowrite)
{
  myfile << stringtowrite << std::flush;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main class RESULT implementations

void result::save()
{
  OpenBabel::OBElementTable etable;
  //std::cout << " In Save" << std::endl;
  int moleculenumber=thecounter->requestnumber();
  if((moleculenumber % 1000)==0){std::cout << "." << std::flush;}
 //printatoms();
 //printbonds();

   for(int i=999;i>=0;i--)
   {
       if(accounting[i]==DELEH)       
       {
          OpenBabel::OBAtom * h = themolecule.GetAtom(i+1);//GetAtom uses 1 based index
          themolecule.DeleteAtom(h);
          //std::cout << "Deleting atom " << i << std::endl;
       }
   }
 //  printatoms();
 //  printbonds();
 //themolecule.ConnectTheDots();//old way of doing things - remember for comparison with old versions
 //themolecule.DeleteHydrogens();
 //themolecule.PerceiveBondOrders();
 //themolecule.AddHydrogens();

  char formattednumber[100];
  sprintf(formattednumber, "%010d", moleculenumber);
  //we want to convert this to a char* or std:string that looks like "00000027.xyz or 00000036.xyz
  // later we look at desired output format to create this filename from the formatted number
  //std::cout << " The file name is" <<  filename << std::endl;

  thelogfilemaker->writetologfile(formattednumber);
  themolecule.SetTitle(formattednumber);
  //std::cout << " About to write to log file " << std::endl;
  thelogfilemaker->writetologfile(" : ");
 
  OpenBabel::OBConversion obconversionsmi,obconversionsdf,obconversion;
  //std::cout << " After OBConversion " << std::endl;
  obconversionsmi.SetOutFormat("smi");
  obconversionsdf.SetOutFormat("sdf");
  

  std::stringstream buffersmi,buffersdf;
  obconversionsmi.Write(&themolecule,&buffersmi);
  std::string smilesstring;
  buffersmi >> smilesstring;//this puts the smiles string into filepart
  thelogfilemaker->writetologfile(smilesstring);

  //std::cout << " about to confsearch with value " << this->confsearch << std::endl;
  if(this->confsearch>0)
  {
     OpenBabel::OBConformerSearch cse;
     OpenBabel::OBForceField *pFF = OpenBabel::OBForceField::FindForceField("MMFF94");
     //OpenBabel::OBConformerScore * csc = OpenBabel::OBConformerScore();
     if(confsearch==1){cse.Setup(themolecule,10,5,5,2);}
     else if(confsearch==2){cse.Setup(themolecule,20,5,5,10);}
     else if(confsearch==3){cse.Setup(themolecule,30,5,5,25);}
     else{cse.Setup(themolecule,40,7,7,30);}
     cse.Search();     
     cse.GetConformers(themolecule);
     
     int nc=themolecule.NumConformers();
     //std::cout << "Number of conformers is " << nc << std::endl;
     float mine=10000000;
     int minen=0;
     for(int i=0;i<nc;i++)
     {
        themolecule.SetConformer(i);
        pFF->Setup(themolecule);
        float energy=pFF->Energy();
        if(energy<mine){mine=energy;minen=i;}
        //std::cout << "Energy for conformer " << i << " is " << themolecule.GetEnergy(i) << ":" << energy << std::endl;

     }
     themolecule.SetConformer(minen);
  }


  if(this->optimize)
  {
    OpenBabel::OBForceField *pFF = OpenBabel::OBForceField::FindForceField("MMFF94");
    if(pFF->Setup(themolecule))
    {
     //std::cout << "Minimizing" << std::endl;
     //pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);
     //pFF->SetLogFile(&std::cerr);
     pFF->ConjugateGradients(1000);
     double TotE=pFF->Energy();
     thelogfilemaker->writetologfile(", MMFF ENERGY: ");   
     std::stringstream buffer;
     buffer.setf(std::ios_base::fixed);
	  buffer.precision(6);
	  buffer.width(8);
     buffer << TotE;
     std::string mystring=buffer.str();
     thelogfilemaker->writetologfile(mystring);
    }
  }

  
   formula=themolecule.GetSpacedFormula(1,"",false);
   themolecule.SetFormula(formula);//this crashes everything when OBMols are destroyed, would like to uncomment but can't.
                                      //tried static objects, dynamic object creation and destruction, just live with wrong formula
                                      //in sdf, and print right stochiometry to log file
   thelogfilemaker->writetologfile(", Formula: "+formula);   
  
   if(this->physico)
   {
     this->calculatedescriptors();
     thelogfilemaker->writetologfile(", HBA: ");
     thelogfilemaker->writetologfile(this->HBA);
     thelogfilemaker->writetologfile(", HBD: ");
     thelogfilemaker->writetologfile(this->HBD);
     thelogfilemaker->writetologfile(", MW: ");
     thelogfilemaker->writetologfile(this->MW);
     thelogfilemaker->writetologfile(", LogP: ");
     thelogfilemaker->writetologfile(this->LOGP);
     thelogfilemaker->writetologfile(", TPSA: ");
     thelogfilemaker->writetologfile(this->TPSA);
   }
   //printatoms();
   //printbonds();

  if(theoutputformat.compare("none")!=0)
  {
    obconversion.SetOutFormat(theoutputformat.c_str());
    std::string filename=theoutputdirectory+"/"+ formattednumber + "." + theoutputformat;
    obconversion.WriteFile(&themolecule,filename);
  }
  
  thestructures->writetofile(obconversionsdf.WriteString(&themolecule));
  thelogfilemaker->writetologfile("\n");

  //could manually write in xyz format  
  /*int numberofatoms=themolecule.NumAtoms();
  char numatoms[100];
  sprintf(numatoms,"%d\n",numberofatoms);
  thestructures->writetofile(numatoms);
  thestructures->writetofile(formattednumber);
  thestructures->writetofile("\n");

  for( OpenBabel::OBMolAtomIter a(&themolecule); a; ++a )
  {
      std::stringstream buffer;
	  buffer.setf(std::ios_base::fixed);
	  buffer.precision(6);
	  buffer.width(2);
     buffer.unsetf(std::ios_base::fixed);
	  buffer << etable.GetSymbol(a->GetAtomicNum());
     buffer.setf(std::ios_base::fixed);
	  buffer << " " ;buffer.width(9);buffer << a->x();buffer << " ";buffer.width(9);buffer << a->y();buffer << " ";buffer.width(9);buffer << a->z() << std::endl;
     std::string mystring=buffer.str();
     thestructures->writetofile(mystring);
  }*/
}


void result::add(OpenBabel::OBMol &other, int hydrogenatomtodelete)
{
   OpenBabel::OBAtom * c;
   OpenBabel::OBAtom * h;

  //std::cout << "in add -- NB the hydrogen we are substituting is int " << hydrogenatomtodelete << std::endl;
  //std::cout << "lets see the substituent " << other.GetFormula() << std::endl;  

   //rotation/translation tricky
   if((themolecule.NumAtoms()>0)&&(other.NumAtoms()>0))
   {
      //we need to know the CH (or NH vector that is proceeding)
      h = themolecule.GetAtom(hydrogenatomtodelete+1);//GetAtom uses 1 based index
 
      //h is easy, because we specify hydrogenatomtodelete
     
      //go through the rest of the atoms in mol, and if not hydrogen and close to the
      //hydrogen we know the other atom on the bond
      for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
      {
        if(!((&*a)->IsHydrogen())&&(AtomDistance(&*a,h)<1.2)){c=(&*a);}
      }
      //now we know C and H that is going to be added to
      OpenBabel::OBBond *thebond = c->GetBond(h);
      int bondno = thebond->GetIdx();
      //printbonds();
      //std::cout << "deleting hl " << c->GetIdx() << "-" << h->GetIdx() << std::endl;
      themolecule.DeleteBond(thebond);
      //printbonds();

      //for the substituent you made it easy for me. The H is atom 1
      //now look for the next atom
      OpenBabel::OBAtom * hs = other.GetAtom(1);
      //std::cout << "We have the hydrogen and its index (zero based) is " << hs->GetIndex() 
      //          << "and to check its mass is " << hs->GetAtomicMass()  << std::endl;
      //std::cout << "hs has coords :" << hs->GetX() << "," << hs->GetY() << "," << hs->GetZ() << std::endl; 
      OpenBabel::OBAtom * cs = other.GetAtom(2);
  
      //std::cout << "We have the heavy atom and its index (zero based) is " << cs->GetIndex() 
      //          << "and to check its mass is " << cs->GetAtomicMass()  << std::endl;
      //std::cout << "cs has coords :" << cs->GetX() << "," << cs->GetY() << "," << cs->GetZ() << std::endl;

      //At this point we have our substituent XH vector   
      //create vector from, and vector to
      //int contacts=badcontacts(&themolecule,&other);
      //std::cout << " we have bad contacts " << contacts << std::endl;
      Vector3 from = Vector3(hs->GetX(),hs->GetY(),hs->GetZ());
      //std::cout << "after trying to access hs" << std::endl;
      Vector3 to = Vector3(0.000,0.000,0.000);
      //move the substituent to the origin before rotation
      Translate(&other,from,to);

      //Rotate the substituent appropriately using the CH and XH vectors
      Vector3 main=Vector3(c->GetX()-h->GetX(),c->GetY()-h->GetY(),c->GetZ()-h->GetZ());

      //in normalize
      //main.print();
      main.normalize();

      //std::cout << "about to create vector subst" << std::endl;
      //std::cout << "we create subst from hs and cs x y and z coords" << std::endl;

      //std::cout << "cs has coords :" << cs->GetX() << "," << cs->GetY() << "," << cs->GetZ() << std::endl;
      //std::cout << "hs has coords :" << hs->GetX() << "," << hs->GetY() << "," << hs->GetZ() << std::endl;            
      Vector3 subst=Vector3(hs->GetX()-cs->GetX(),hs->GetY()-cs->GetY(),hs->GetZ()-cs->GetZ());

      subst.normalize();
  
      //Rotate to the appropriate orientation
      Rotate(&other,main,-subst);
      //std::cout << "cs has coords :" << cs->GetX() << "," << cs->GetY() << "," << cs->GetZ() << std::endl;
      //std::cout << "hs has coords :" << hs->GetX() << "," << hs->GetY() << "," << hs->GetZ() << std::endl;   

      //translate to the appropriate place
      to = Vector3(c->GetX()+0.297*(h->GetX()-c->GetX()),c->GetY()+0.297*(h->GetY()-c->GetY()),c->GetZ()+0.297*(h->GetZ()-c->GetZ()));
      from = Vector3(-hs->GetX(),-hs->GetY(),-hs->GetZ());
      Translate(&other,from,to);

      //int minimumangle=0;
      //int minimumbad=10000;
   }

   //NB this is to add substituent to the molecule   
   //NB other is the substituent, and we only will add from atom number 1
   //we included atom number 0 as H to get the H-X vector for alignment of fragment.
   int currentatom=0;
   float massdifference=0.0;
   float donordifference=0.0;
   float acceptordifference=0.0;
   //printbonds();
   //OpenBabel::addFragment(themolecule,other)
   
   //OpenBabel::OBElementTable etable;
   OpenBabel::OBMol newfrag=other;
   newfrag.DeleteAtom(newfrag.GetAtom(1));
   int noatoms=themolecule.NumAtoms(); 
   //printatoms();

   themolecule+=newfrag;

   //printatoms()
   //int newnoatoms=themolecule.NumAtoms(); 
   //std::cout << "new number of atoms " << newnoatoms << std::endl;
   //printbonds();

   //explicitly create the new bond
   OpenBabel::OBBond * newbond;
   newbond=themolecule.NewBond();
   newbond->SetParent(&themolecule);
   newbond->SetBondOrder(1);
   newbond->SetBegin(c);
   OpenBabel::OBAtom * endatom=themolecule.GetAtom(noatoms+1);
  
   //std::cout << "start atom is " << etable.GetSymbol(c->GetAtomicNum()) << std::endl;
   //std::cout << "end atom is " << etable.GetSymbol(endatom->GetAtomicNum()) << std::endl;
   //std::cout << "now bonding " << c->GetIdx() << "-" << noatoms+1 << std::endl;
   newbond->SetEnd(endatom);
   //printbonds();
   //

   //alternative is to work atom by atom
/* for(OpenBabel::OBMolAtomIter a(other);a;a++)
    {
      //need to duplicate all atoms from the subtituent and add them to the current molecule
      //note that we drew the subst with H first, so that we could ignore the first hydrogen
      //e.g. drew HCl to substitute Cl
      //std::cout << "Adding atoms to the molecule, atom number " << currentatom << std::endl;
      if(currentatom>=1)
      {
        OpenBabel::OBAtom newatom;
        newatom.Duplicate(&*a);
        themolecule.InsertAtom(newatom);
        massdifference+=newatom.GetAtomicMass();
        if(currentatom==1)
        {
           OpenBabel::OBBond * newbond;
           newbond=themolecule.CreateBond();
           newbond->SetBondOrder(1);
           newbond->SetBegin(c);
           newbond->SetEnd(&newatom);
           std::cout << "now bonding " << c->GetIndex() << "-" << newatom.GetIndex() << std::endl;
        }

      }
      
      
       currentatom++;
    }   
*/

//for dynamic filtering
//massdifference-=1.0079;
//addapproxmass(massdifference);
/*
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("HBD");
    if(pDescr)
    {
      HBD+=pDescr-> Predict(&other);
    }
    pDescr = OpenBabel::OBDescriptor::FindType("HBA1");
    if(pDescr)
    {
      HBA1+=pDescr-> Predict(&other);
    }
    pDescr = OpenBabel::OBDescriptor::FindType("HBA2");
    if(pDescr)
    {
      HBA2+=pDescr-> Predict(&other);
    }
*/

//h = themolecule.GetAtom(hydrogenatomtodelete+1);//GetAtom uses 1 based index
accounting[hydrogenatomtodelete]=DELEH;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main class RESULT minor functions

void result::print()
   {
      std::cout << "molecular mass " << themolecule.GetExactMass() << std::endl;
      for(int i=0;i<50;i++)std::cout << accounting[i];
      std::cout << std::endl;
   }

void result::printatoms()
{
   OpenBabel::OBElementTable etable;
  for(OpenBabel::OBMolAtomIter a(themolecule);a;a++)
  {
    std::cout << a->GetIdx() << etable.GetSymbol(a->GetAtomicNum()) << "|";         
  }
  std::cout << std::endl;
  //std::cout << "number of atoms " << noatoms << std::endl;
}  

void result::printbonds()
{
   for(OpenBabel::OBMolBondIter b(themolecule);b;b++)
   {
    int begina=b->GetBeginAtomIdx();        
    int enda=b->GetEndAtomIdx();
    std::cout << begina << ":" << enda << "|" ;        
   }
   std::cout << std::endl;


}
 
void result::resetall(OpenBabel::OBMol moleculewewantinobject)
   {
      themolecule=moleculewewantinobject;
      for(int j=0;j<1000;j++)setacc(NOTAH,j);
      for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
      {      
        int i =(&*a)->GetIndex();
        if((&*a)->IsHydrogen()){setacc(IGNOH,i);}
      }

     //FOR_BONDBFS_OF_MOL(b, themolecule)
     // {
     //   OpenBabel::OBBond * mybond = &*b;
     //   themolecule.DeleteBond(mybond,false); 	
         // The variable b behaves like OBBond* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

     // }

   }

void result::set_substupto_and_maxh()
{
      substupto=-1;
      maxh=-1;
      int found=0;
      OpenBabel::OBElementTable etable;
      for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
      {
        int i = (&*a)->GetIndex();        
        //std::cout << "Index " << i << " atom type " << etable.GetSymbol(a->GetAtomicNum()) << " accounting " << accounting[i] << std::endl;
        if(((&*a)->IsHydrogen())&&(found==0)&&(accounting[i]==ORIGH)){substupto=i;found=1;}
        if(((&*a)->IsHydrogen())&&(accounting[i]==ORIGH)){maxh=i;}
      }
      if(substupto==-1)
      { 
         std::cout << "problem -- no hydrogens found " << std::endl;
         std::cout << " when choosing random this is a possible outcome. If this is the case then simply rerun" << std::endl;
         exit(1);
      }
      else
      {
         //std::cout << "substupto is at value " << substupto << std::endl;
         //std::cout << "maxh is at value " << maxh << std::endl;
      }
      //at this point we have set substupto (being the first hydrogen to substituted
      //and maxh as being the last hydrogen to substitute.
 
}

void result::ignorehydrogen(int i)
{
        setacc(IGNOH,i);
}

void result::substitutablehydrogen(int i)
{
        setacc(ORIGH,i);
}

int result::count(int level)
   {
      int thevalue=0;
      for(int i=0;i<1000;i++){if(accounting[i]==level)thevalue++;};
      return thevalue;
   }

/////////
void result::calculatedescriptors()
{
    OpenBabel::OBDescriptor* pDescr_hbd = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDescr_hba = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDescr_logp = OpenBabel::OBDescriptor::FindType("logP");
    OpenBabel::OBDescriptor* pDescr_mw = OpenBabel::OBDescriptor::FindType("MW");
    OpenBabel::OBDescriptor* pDescr_tpsa = OpenBabel::OBDescriptor::FindType("TPSA");
    //all or nothing
    if(pDescr_hbd && pDescr_hba && pDescr_logp && pDescr_mw && pDescr_tpsa)
    {
       HBD = pDescr_hbd-> PredictAndSave(&themolecule);
       HBA = pDescr_hba-> PredictAndSave(&themolecule);
       LOGP = pDescr_logp-> PredictAndSave(&themolecule);
       MW = pDescr_mw-> PredictAndSave(&themolecule);
       TPSA = pDescr_tpsa-> PredictAndSave(&themolecule);
       //std::cout << "****** " <<"The no of Hydrogen Bond Donors -HBD- is " << HBD << " ****** "<<std::endl;
    }
}

void result::calculateformula()
{
   //themolecule.SetFormula("C10H10O10N10");
   std::cout << "about to calculate formula" << std::endl;
   std::string formula=themolecule.GetSpacedFormula(1,"",false);
   std::cout << formula << std::endl;
   //themolecule.SetFormula(themolecule.GetSpacedFormula(1,"",false));
   //themolecule.GetFormula();
   std::cout << "finished calculating and setting formula" << std::endl;
     //OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("formula");
    //if(pDescr)
    //{
        //themolecule.SetFormula();

      // std::cout << "****** " <<"formula" << pDescr->Predict(&themolecule) << " ****** "<<std::endl;
   // }
}

/*void result::setHBD()
   { 
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("HBD");
    if(pDescr)
    {
       HBD = pDescr-> Predict(&themolecule);
       std::cout << "****** " <<"The no of Hydrogen Bond Donors -HBD- is " << HBD << " ****** "<<std::endl;
    }
   }
*/
/*int result::addnoofhydrobondons(int a)
   {
    noofhydrobondons+=a;
    return noofhydrobondons;
   }


void result::setHBA1()
   { 
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("HBA1");
    if(pDescr)
    {
       HBA1 = pDescr-> Predict(&themolecule);
       std::cout << "****** " <<"The no of Hydrogen Bond Acceptors 1 -HBA1- is " << HBA1 << " ****** "<<std::endl;
    }
   }

void result::setHBA2()
   { 
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("HBA2");
    if(pDescr)
    {
       HBA1 = pDescr-> Predict(&themolecule);
       std::cout << "****** " <<"The no of Hydrogen Bond Acceptors 2 -HBA2- is " << HBA2 << " ****** "<<std::endl;
    }
   }

void result::addHBD(float additionaldonors)
{
    HBD+=additionaldonors;
}

void result::addHBA1(float additionalacceptors)
{
   HBA1+=additionalacceptors;
}

void result::addHBA2(float additionalacceptors)
{
   HBA2+=additionalacceptors;
}


void result::addapproxmass(float additionalmass)
{
   approxmass=additionalmass;
}
*/

//   *****predlogp*****

/*void result::setpredlogp()
   { 
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("logP");
    if(pDescr)
    {
       logP = pDescr-> Predict(&themolecule);
       std::cout << "****** " <<"The octanol/water partition coefficient -logP- is " << logP << " ****** "<<std::endl;
    }
   }
*/

/*void result::deleteallbondsinmolecule()
{
   //delete all bonds in themolecule
  for(OpenBabel::OBMolBondIter b(themolecule);b;++b)
  {  
    themolecule.DeleteBond(&*b);
  }      
}*/

void result::deletedeletablehydrogens()
{
   themolecule.DeleteHydrogens();
   //todo change deletion based on accounting DELEH
   //for(int t=999;t>=0;t--)
   //{
     //delete all hydrogens (marked in accounting as DELEH)
   //}
}


void result::printacc()
{
     for (int q=0;q<50;q++)
     {
       std::cout << accounting[q] << " ";
     }
     std::cout << std::endl;
}

void result::makeonlyringhydrogenssubstitutable() //make only ring hydrogens substitutable as ORIGH
{
   
   // This is a trial to use the OBMol ring iterator based on thisresult?
   // provided that #include <openbabel/obiter.h> and #include <openbabel/mol.h> have been included 
   for (OpenBabel::OBMolRingIter r(themolecule);r;++r)
   {

     for (OpenBabel::OBMolAtomIter a(themolecule);a;++a)
     {
            int i=(*a).GetIdx();
            if ((*r).IsInRing(i))
            {
                for(OpenBabel::OBAtomAtomIter h(*a);h;++h)
                {
                  if((*h).IsHydrogen())
                   {
                      int j = (*h).GetIndex();
                      OpenBabel::OBAtom * c;
                      //go through the rest of the atoms in mol, and if not hydrogen and close to the
                      //hydrogen we know the other atom on the bond
                      for(OpenBabel::OBMolAtomIter b(themolecule);b;b++) 
                      {
                         if(!(b->IsHydrogen())&&(AtomDistance(&*h,&*b)<1.2)){c=(&*b);}
                      }
                      if(this->chonly)
                      {
                        if(c->IsCarbon()){setacc(ORIGH,j);};
                      }
                      else
                      {
                        setacc(ORIGH,j);
                      }         
                      
                   }
               }
            }
     }
    }

}

void result::makeallhydrogenssubstitutable()
{
     int i=-1;
      for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
      {
        i++;
        if((&*a)->IsHydrogen())
        {
             OpenBabel::OBAtom * c;
            //go through the rest of the atoms in mol, and if not hydrogen and close to the
            //hydrogen we know the other atom on the bond
            for(OpenBabel::OBMolAtomIter b(themolecule);b;b++) 
            {
              if(!((&*b)->IsHydrogen())&&(AtomDistance(&*a,&*b)<1.2)){c=(&*b);}
            }
            if(this->chonly)
            {
              if(c->IsCarbon()){setacc(ORIGH,i);};
            }
            else
            {
              setacc(ORIGH,i);
            }         
        }
      } 
}

void result::intelligentlymakehydrogenssubstitutable()
{
     int i=-1;
      for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
      {
          i++;
          if((&*a)->IsHydrogen())
          {
             OpenBabel::OBAtom * c;
            //go through the rest of the atoms in mol, and if not hydrogen and close to the
            //hydrogen we know the other atom on the bond
            for(OpenBabel::OBMolAtomIter b(themolecule);b;b++) 
            {
              if(!((&*b)->IsHydrogen())&&(AtomDistance(&*a,&*b)<1.2)){c=(&*b);}
            }

            if(c->IsCarbon()){setacc(ORIGH,i);};
          }
      } 
}


void result::makerandomhydrogenssubstitutable(int numberofrandom)
{
  //? how many hydrogens are available?? =N
  int N=0;
  for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
  {    
        if((&*a)->IsHydrogen()){N++;}
  } 
  //std::cout<< "We have a total of " << N << " hydrogen atoms" << std::endl;
  //? proportion= 1000*8/N //oldversion
  int proportion= numberofrandom*10;//convert % to be out of 1000
  //int proportion = (8000)/N;

  int i=-1;
  srand(time(0));
  for(OpenBabel::OBMolAtomIter a(themolecule);a;a++) 
  {
         i++;
         int random = rand() % 1000;        
         //std::cout << "we have " << N << " hydrogens and proportion is " << proportion << " and random is " << random << std::endl;
         if((&*a)->IsHydrogen() && (random<proportion))
         {
            setacc(ORIGH,i);
         }
  } 
  //generate random number between 1 and 1000
  // if number < proportion, then it must be substitutable.
}

result & result::operator= (result &other)
   {
        themolecule=other.themolecule;
        for(int i=0;i<1000;i++)
        {
           accounting[i]=other.getacc(i);
        }
        substupto=other.substupto;
        maxh=other.maxh;
        // copy the classes into the result class

        thecounter=other.thecounter;
        thelogfilemaker=other.thelogfilemaker;
        theerrorfilemaker=other.theerrorfilemaker;
        thestructures=other.thestructures;
        
        
        //noofhydrobondons=other.noofhydrobondons;
        //noofhydrobonaccs=other.noofhydrobonaccs;
        //logP=other.logP;
        //HBA1=other.HBA1;
        //HBA2=other.HBA2;
        //HBD=other.HBD;
        //approxmass=other.approxmass;
        theoutputdirectory=other.theoutputdirectory;
        theoutputformat=other.theoutputformat;
        maximumsubs=other.maximumsubs;
        maximumoutput=other.maximumoutput;
        optimize=other.optimize;
        physico=other.physico;
        nocurrentlysubst=other.nocurrentlysubst;
        force=other.force;
        confsearch=other.confsearch;
        
        return *this;
   };
