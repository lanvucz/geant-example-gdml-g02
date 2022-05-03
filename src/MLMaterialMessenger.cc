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
#include "MLMaterialMessenger.hh"

#include "MLMaterial.hh"
#include "G4NistManager.hh"

#include <sstream>

////////////////////////////////////////////////////////////////////////////////
//
MLMaterialMessenger::MLMaterialMessenger (MLMaterial * MLMat)
  :materialsManager(MLMat)
{ 
  MaterialDir = new G4UIdirectory("/geometry/material/");
  MaterialDir->SetGuidance(" Controls for defining geometry materials" );

  AddCmd = new G4UIcommand("/geometry/material/add",this);
  AddCmd->SetGuidance(
    "  add a mateial by name, composition formula and density");
  AddCmd->SetGuidance("  name: e.g. water ");
  AddCmd->SetGuidance("  formula (e.g. H2-O for water");
  AddCmd->SetGuidance("  density (in units of g/cm3) : den>0.");
  G4UIparameter* MatName = new G4UIparameter("material",'s',false);
  MatName->SetGuidance("material name");
  AddCmd->SetParameter(MatName);
  //
  G4UIparameter* MatForm = new G4UIparameter("formula",'s',false);
  MatForm->SetGuidance("material formula");
  AddCmd->SetParameter(MatForm);
  //    
  G4UIparameter* DenPrm = new G4UIparameter("density",'d',false);
  DenPrm->SetGuidance("density of the material");
  DenPrm->SetParameterRange("density >0.");
  AddCmd->SetParameter(DenPrm);
  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* StatePrm = new G4UIparameter("state",'s',true);
  StatePrm->SetGuidance("state of the material (optional): gas | solid");
  AddCmd->SetParameter(StatePrm);
  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* TempPrm = new G4UIparameter("temp",'d',true);
  TempPrm->SetGuidance("temperature of the material in Kelvin (optional)");
  AddCmd->SetParameter(TempPrm);
  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* PresPrm = new G4UIparameter("pres",'d',true);
  PresPrm->SetGuidance("pressure of the gas material in Pascal (optional)");
  AddCmd->SetParameter(PresPrm);
  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  AddNISTCmd = new G4UIcmdWithAString("/geometry/material/addNIST",this);
  AddNISTCmd->SetGuidance("Add a predefined Geant4 NIST material");
  AddNISTCmd->SetParameterName("NISTmaterial",false);
  AddNISTCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  //
  ListNISTCmd = new G4UIcmdWithAString("/geometry/material/listNIST",this);
  ListNISTCmd->SetGuidance("List the predefined Geant4 NIST materials");
  ListNISTCmd->SetParameterName("type", true, true);
  ListNISTCmd->SetDefaultValue( "all");
  ListNISTCmd->SetCandidates( "all simple compound hep");
  ListNISTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  DeleteIntCmd = new G4UIcmdWithAnInteger("/geometry/material/delete",this);
  DeleteIntCmd->SetGuidance("Delete material by its index");
  DeleteIntCmd->SetParameterName("matIdx",false);
  DeleteIntCmd->SetRange("matIdx>=0 && matIdx<100");
  DeleteIntCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  DeleteNameCmd = new G4UIcmdWithAString("/geometry/material/deleteName",this);
  DeleteNameCmd->SetGuidance("Delete material by its name.");
  DeleteNameCmd->SetParameterName("DeleteName",false);
  DeleteNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //  
  ListCmd = new G4UIcmdWithoutParameter("/geometry/material/list",this);
  ListCmd->SetGuidance("List the materials defined");
  ListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}
////////////////////////////////////////////////////////////////////////////////
//
MLMaterialMessenger::~MLMaterialMessenger ()
{
  delete MaterialDir;
  delete AddCmd;
  delete AddNISTCmd;
  delete DeleteIntCmd;
  delete DeleteNameCmd;
  delete ListCmd;
  delete ListNISTCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLMaterialMessenger::SetNewValue (G4UIcommand* command,G4String newValue)
{    
  if (command == DeleteIntCmd) {
    materialsManager->DeleteMaterial(DeleteIntCmd->GetNewIntValue(newValue));

  } else if (command == DeleteNameCmd) {
    materialsManager->DeleteMaterial(newValue);

  } else if (command == ListCmd) {
    materialsManager->ListMaterial();

  } else if (command == AddCmd) {
    G4double den, tem, pres ;
    G4String state;
    char mat[80], form[80], stat[10];
    stat[0] = ' ';
    tem = pres = -1.;
    const char* t = newValue;
    std::istringstream is((char*)t);
    is >>mat >>form >>den >>stat >> tem >> pres ;
    G4String material=mat;
    G4String formula=form;
    if (pres == -1.) { 
      state = "";
    } else {
      state = stat;
    }
    //    G4cout<< "stat = " <<state<< "tem = " << tem<< " pre = " << pres << G4endl;
    //     tick *= G4UIcommand::ValueOf(unt);
    materialsManager->AddMaterial(material,formula,den,state,tem,pres);
  } else if (command == AddNISTCmd) {
    materialsManager->AddNISTMaterial( newValue);
  } else if (command == ListNISTCmd) {
    G4String n = newValue;
    if ( n.length() == 0) n = "all";
    G4NistManager::Instance()->ListMaterials( n);
  }
}
////////////////////////////////////////////////////////////////////////////////
