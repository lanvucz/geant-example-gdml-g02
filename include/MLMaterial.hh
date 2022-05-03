//###################################################################################

// (C) Copyright European Space Agency, 2019
// 
// This file is subject to the terms and conditions defined in file 'LICENCE.txt', 
// which is part of this source code package. No part of the package, including 
// this file, may be copied, modified, propagated, or distributed except 
// according to the terms contained in the file ‘LICENCE.txt’.“ 

//###################################################################################
#ifndef MLMaterial_HH
#define MLMaterial_HH
////////////////////////////////////////////////////////////////////////////////
//
#include "G4Material.hh"
#include <vector>

class MLMaterialMessenger;
////////////////////////////////////////////////////////////////////////////////
//
class MLMaterial
{
public:

  MLMaterial ();
  ~MLMaterial ();

public:

  void  AddMaterial (G4String, G4String, G4double, G4String, G4double, G4double );
  void  AddNISTMaterial(G4String);
  G4Material* GetMaterial (G4int i)  {return Material[i];};
  G4Material* GetMaterial (G4String name)
    {return G4Material::GetMaterial(name);} ;
  G4int GetMaterialIndex (G4String);
  G4int GetNbOfMaterial () {return Material.size();};
  void  DeleteMaterial (G4int);
  void  DeleteMaterial (G4String);

  void  ListMaterial();

private:

  MLMaterialMessenger         *materialMessenger;

  std::vector<G4Material*>   Material;
  std::vector<G4Element*>    Element;
  std::vector<G4Isotope*>    Isotope;

private:
  static const G4String        ELU[110];
  static const G4String        ELL[110];
  static const G4String        EUU[110];
  static const G4double        A[110];
       
};
////////////////////////////////////////////////////////////////////////////////
#endif
