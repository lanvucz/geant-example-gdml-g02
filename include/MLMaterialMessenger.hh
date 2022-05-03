//###################################################################################

// (C) Copyright European Space Agency, 2019
// 
// This file is subject to the terms and conditions defined in file 'LICENCE.txt', 
// which is part of this source code package. No part of the package, including 
// this file, may be copied, modified, propagated, or distributed except 
// according to the terms contained in the file ‘LICENCE.txt’.“ 

//###################################################################################
#ifndef MLMaterialMessenger_h
#define MLMaterialMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

class MLMaterial;
////////////////////////////////////////////////////////////////////////////////
//
class MLMaterialMessenger: public G4UImessenger
{
public:
  MLMaterialMessenger(MLMaterial* );
  ~MLMaterialMessenger();

  void SetNewValue (G4UIcommand*, G4String);

private:

  MLMaterial                *materialsManager;

  G4UIdirectory             *MaterialDir;
  G4UIcmdWithoutParameter   *ListCmd;
  G4UIcmdWithAnInteger      *DeleteIntCmd;
  G4UIcmdWithAString        *DeleteNameCmd;
  G4UIcommand               *AddCmd;
  G4UIcmdWithAString        *AddNISTCmd;
  G4UIcmdWithAString        *ListNISTCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
