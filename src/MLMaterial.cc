//###################################################################################
//
// (C) Copyright European Space Agency, 2019
// 
// This file is subject to the terms and conditions defined in file 'LICENCE.txt', 
// which is part of this source code package. No part of the package, including 
// this file, may be copied, modified, propagated, or distributed except 
// according to the terms contained in the file ‘LICENCE.txt’.“ 
//
// This file forms part of the Mulassis application, available from the 
// European Space Software Repository (ESSR): https://essr.esa.int/
//
//###################################################################################


////////////////////////////////////////////////////////////////////////////////
//
#include "MLMaterial.hh"
#include "MLMaterialData.hh"
#include "MLMaterialMessenger.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"

#include <vector>
#include <iomanip>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLMaterial::MLMaterial ()
{
  Material.clear();
  Element.clear();
  Isotope.clear();
  // some default materials vacuum (0), air (1)  and aluminium (2)  defined here
  // examples of vacuum
  //
  G4double a,z;

  static G4bool bmat = false ;

  if (!bmat) {

    // vacuum
    G4double  density    = universe_mean_density;    //from PhysicalConstants.h
    G4double pressure    = 3.e-18*pascal;
    G4double temperature = 2.73*kelvin;
    G4Material* aMaterial= new G4Material("Vacuum", z=1., a=1.01*g/mole, 
      density, kStateGas,temperature,pressure);
    Material.push_back(aMaterial);

    // air
    a                   = 14.01*g/mole;
    G4Element* aElement = new G4Element(" N"," N", 7., a);
    Element.push_back(aElement);
    a        = 16.00*g/mole;
    aElement = new G4Element(" O", " O", 8., a);
    Element.push_back(aElement);
    
    density   = 1.290*mg/cm3;
    aMaterial = new G4Material("Air", density, 2);
    aMaterial->AddElement(Element[0],0.78);
    aMaterial->AddElement(Element[1],0.22);
    Material.push_back(aMaterial);
    
    // aluminium
    G4Element* Al = new G4Element("Al", "", z= 13., a=26.98*g/mole);
    aMaterial = new G4Material("Aluminium", density=2.700*g/cm3, 1);
    aMaterial->AddElement(Al, 1);
    Material.push_back(aMaterial);
    
     //silicon
    G4Element* Si = new G4Element("Si", "", z= 14., a=28.0855*g/mole);
    aMaterial = new G4Material("Silicon", density=2.3290*g/cm3, 1);
    aMaterial->AddElement(Si, 1);
    Material.push_back(aMaterial);

    G4Element* C = new G4Element(" C"," C", 6., a=12.0107*g/mole);
    Element.push_back(C);

    G4Element* H = new G4Element(" H"," H", 1., a=1.00794*g/mole);
    Element.push_back(H);

    // create commands for interactive definition of the calorimeter  
    materialMessenger = new MLMaterialMessenger(this);
    bmat              = true;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
MLMaterial::~MLMaterial ()
{
  delete materialMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterial::AddMaterial (G4String name, G4String formula, G4double density,
			      G4String state, G4double tem, G4double pres)
{
  //
  //  G4cout << " Nb of Isotope  = " << G4Isotope::GetNumberOfIsotopes() << endl;
  //  G4cout << " Nb of Element  = " << G4Element::GetNumberOfElements() << endl;
  //  G4cout << " Nb of Material = " << G4Material::GetNumberOfMaterials() << endl;
//  G4cout << " Material properties  = name:" << name << " formula:" <<
//  formula <<  " density:" << density << " state:" << state << " tem:" << tem << " pres:" << pres << G4endl;
  G4int isotope, Z;
  size_t i;
  for (i = 0; i<Material.size(); i++) {
    if (Material[i]->GetName() == name) {
      G4cerr <<" AddMaterial : material " <<name
             <<" already exists." <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
      return;
    }
  }

  char *tokenPtr1 = NULL;
  char *sname     = NULL;
  G4String str, s1(".0123456789");
  G4String element, isotopename;
  G4int ncomponents, natoms;
  G4double atoms;
  size_t ls, id, ll, lr;
  ncomponents = 0;

  sname       = new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  while (tokenPtr1 != NULL) {
    ncomponents++;
    tokenPtr1 = strtok( NULL, "-");
  }
  delete[] sname;

  G4Material* aMaterial = 0;

  if (state == "") {
    aMaterial = new G4Material(name, density*g/cm3, ncomponents);
  } else if (state == "solid" && tem > 0.) {
    aMaterial = new G4Material(name, density*g/cm3, ncomponents, 
					   kStateSolid, tem * kelvin );
  } else if (state == "gas" && pres > 0.) {
    aMaterial = new G4Material(name, density*g/cm3, ncomponents, 
					   kStateSolid, tem * kelvin, pres*pascal );
  }
  if (aMaterial == 0) {
    G4cerr <<" AddMaterial : Name " <<name <<"." <<G4endl;
    G4cerr <<"--> Command failed." <<G4endl;
    return;
  }

  sname=new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  if (tokenPtr1 == NULL) {
    G4cerr <<" AddMaterial : No Elements specified in Formula or Empty Formula: \"" 
	   << formula << "\"" << G4endl;
    G4cerr <<" --> Command failed."<<G4endl;
    return;
  }

  while (tokenPtr1 != NULL) {
    isotope = 0;
    natoms = 0;
    atoms = 0.;
    //          G4cout << tokenPtr1 << G4endl;
    str     = G4String(tokenPtr1);
    ls      = str.length();
    ll      = str.find("(");
    lr      = str.find(")");
    if (ll == lr) {
      id = str.find_first_of(s1);
      element = str.substr(0,id);
      
      if (element.length() == 1) element.insert(0," ");
      for (i = 0; i<110; i++) {
        if (element == ELU[i]) break;
      }
      if (i == 110) {
        for (i = 0; i<110; i++) {
          if (element == ELL[i]) break;
        }
        if (i == 110) {
          for (i = 0; i<110; i++) {
            if (element == EUU[i]) break;
          }
        }
      }
      
      if (i == 110) {
        G4cerr <<"AddMaterial : Invalid element in material formula."
               <<element <<G4endl;
        G4cerr <<"--> Command rejected." <<G4endl;
//        delete aMaterial;
//	Material[NbMat] = NULL;
        return;
      }

      Z       = i+1;
      element = ELU[i];
      if (id == std::string::npos) {
        natoms = 1;
      } else {
        natoms = atoi((str.substr(id,ls-id)).c_str());
	atoms = atof((str.substr(id,ls-id)).c_str());
      }
//  G4cout << "   Elements = " << element << G4endl;
//  G4cout << "   Nb of natoms = " << natoms << " atoms=" <<atoms <<G4endl;
    } else {
      element = str.substr(0,ll);
      isotope = atoi((str.substr(ll+1,lr-ll)).c_str());
      if (element.length() == 1) element.insert(0," ");
      for (i = 0; i<110; i++) {
        if (element == ELU[i]) break;
      }
      if (i == 110) {
        for (i = 0; i<110; i++) {
          if (element == ELL[i]) break;
        }
        if (i == 110) {
          for (i = 0; i<110; i++) {
            if (element == EUU[i]) break;
          }
        }
      }
      if (i == 110) {
        G4cerr <<"AddMaterial : Invalid element in material formula."
               <<element <<G4endl;
        G4cerr <<"--> Command rejected." <<G4endl;
//        delete aMaterial;
//	Material[NbMat] = NULL;
        return;
      }

      Z           = i+1;
      element     = ELU[i];
      isotopename = element+str.substr(ll+1,lr-ll);
      if (lr == std::string::npos ) {
        natoms = 1;
      } else {
        natoms = atoi((str.substr(lr+1,ls-lr)).c_str());
	atoms = atof((str.substr(lr+1,ls-lr-1)).c_str());
      }  
      if (natoms == 0 && atoms == 0.) natoms = 1;
//  G4cout << "   Elements = " << element << G4endl;
//  G4cout << "   Isotope Nb = " << isotope << G4endl;
//  G4cout << "   Nb of atoms = " << natoms << G4endl;
    }
    if (isotope != 0) {
      if (G4Isotope::GetIsotope(isotopename) == NULL) {
	//        G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, A[Z-1]*g/mole);
        G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, isotope*g/mole);
        G4Element* aElement = new G4Element(isotopename, element, 1);
        aElement->AddIsotope(aIsotope, 100.*perCent);
        Isotope.push_back(aIsotope);
        if (natoms > 0) {
	        	aMaterial->AddElement(aElement, natoms);
		} else {
			aMaterial->AddElement(aElement, atoms);
	}	
        Element.push_back(aElement);
      } else {
	if (natoms > 0) {
	        aMaterial->AddElement( G4Element::GetElement(isotopename,false) , natoms);
	} else {
		aMaterial->AddElement( G4Element::GetElement(isotopename,false) , atoms);
	}
      }

    } else {
//        G4cout << "   element= " << element << " Z=" << Z << " A[Z-1]*g/mole=" << A[Z-1] << G4endl;
        if ( G4Element::GetElement(element) == NULL) {
        G4Element* aElement = new G4Element(element, element, Z, A[Z-1]*g/mole);
	if (natoms > 0) {
	        aMaterial->AddElement(aElement, natoms);
	} else {
		aMaterial->AddElement(aElement, atoms);
	}        
	Element.push_back(aElement);
//          G4cout << "   new_added_aElement= " << aElement  << G4endl;
      } else {
//            G4cout << "   existed_aElement= " << G4Element::GetElement(element)  << G4endl;
            if (natoms > 0) {
	        aMaterial->AddElement( G4Element::GetElement(element,false) , natoms);
	} else {
		aMaterial->AddElement( G4Element::GetElement(element,false) , atoms);
	}
      }
    }
    tokenPtr1 = strtok( NULL, "-");
    //      str.empty();
    //element.erase();
    //
  }

  delete[] sname;
  Material.push_back(aMaterial);
  G4cout <<" Material:" <<name <<" with formula: " <<formula <<" added! "
         <<G4endl;
  G4cout <<"     Nb of Material = " <<Material.size() <<G4endl;
  G4cout <<"     Nb of Isotope =  " <<Isotope.size() <<G4endl;
  G4cout <<"     Nb of Element =  " <<Element.size() <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterial::AddNISTMaterial (G4String nist_name)
{
  G4NistManager *man = G4NistManager::Instance();
  G4Material    *aMaterial = 0;

  // Test whether the material is already added to the tables.
  size_t i;
  for (i = 0; i<Material.size(); i++) {
    if (Material[i]->GetName() == nist_name) {
      G4cerr <<" AddMaterial : material " << nist_name
             <<" already exists." <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
      return;
    }
  }


  aMaterial = man->FindOrBuildMaterial( nist_name);
  
  if (!aMaterial ){
    G4cerr << "G4NIST material " << nist_name << " is unknown to the G4NistManager."<<G4endl;
  } else {
    Material.push_back( aMaterial);
    G4cout <<" Material:" << nist_name <<" added from G4NistManager tables." << G4endl
	   <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterial::DeleteMaterial (G4int j)
{
  size_t i(j-1);
  if (i > Material.size()) {
    G4cerr <<"DeleteMaterial : Invalid material index " <<j <<"." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterial::DeleteMaterial (G4String )
{
  G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
         <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLMaterial::GetMaterialIndex (G4String name)
{
  size_t i ;
  for (i = 0; i < Material.size(); i++) {
    if (Material[i]->GetName() == name) break;
  }
  G4int k = G4int(i);
  if (i == Material.size()) k = -1;
  return k;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterial::ListMaterial ()
{
  G4cout <<" There are" <<std::setw(3) <<Material.size()
         <<" materials defined."  <<G4endl;
  for (size_t i = 0; i< Material.size(); i++) 
    G4cout <<"     Material Index " <<std::setw(3) <<i+1 <<" "
           <<std::setw(14) <<Material[i]->GetName()
           <<"  density: " <<std::setw(6) <<std::setprecision(3)
           <<G4BestUnit(Material[i]->GetDensity(),"Volumic Mass") <<G4endl;
  G4cout <<G4endl;
  
}
////////////////////////////////////////////////////////////////////////////////
