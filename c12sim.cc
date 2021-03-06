#include "c12sim.h"

#include <cstdio>
#include <fstream>
#include <vector>
#include <memory>
#include <unistd.h>
#include "getopt.h"
#include "dollar.hpp"

#include "G4SystemOfUnits.hh"
//#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif
#include "G4UImanager.hh"
#include "G4UIQt.hh"
#include "G4UIterminal.hh"
#include "qmainwindow.h"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "FTFP_BERT.hh"
#include "QBBC.hh"
#include "Randomize.hh"
#include "TRandom3.h"

#include "SimulationManager.h"
#include "B1OpticalPhysics.hh"
#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"
#include "B1ParallelWorldConstruction.hh"
#include "G4ParallelWorldPhysics.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "SimplePrimaryGeneratorAction.hh"
//______________________________________________________________________________

int main(int argc,char** argv)
{ $

   int          run_number        = 0;
   int          number_of_events  = -1;
   int          event_start_offset = 0;
   std::string  input_file_name   = "";
   std::string  output_file_name  = "";
   std::string  output_tree_name  = "";
   std::string  output_dir        = "";
   std::string  theRest           = "";
   bool         run_manager_init  = false;
   bool         use_gui           = true;
   bool         use_vis           = true;
   bool         is_interactive    = true;
   bool         has_macro_file    = false;
   int          set_rand_seed     = -1;
   double       torus_field_sign     = 1.0;
   double       solenoid_field_sign  = 1.0;
   double       torus_field_scale    = 1.0;
   double       solenoid_field_scale = 1.0;
   int          simple_eg            = 0;

   //---------------------------------------------------------------------------

   int index = 0;
   int iarg  = 0;
   opterr    = 1; //turn off getopt error message
   const struct option longopts[] =
   {
      {"run",         required_argument,  0, 'r'},
      {"batch",       no_argument,        0, 'b'},
      {"input",       required_argument,  0, 'i'},
      {"output",      required_argument,  0, 'o'},
      {"gui",         required_argument,  0, 'g'},
      {"vis",         required_argument,  0, 'V'},
      {"interactive", no_argument,        0, 'I'},
      {"dir",         required_argument,  0, 'D'},
      {"treename",    required_argument,  0, 't'},
      {"events",      required_argument,  0, 'n'},
      {"offset",      required_argument,  0, 'O'},
      {"help",        no_argument,        0, 'h'},
      {"init",        no_argument,        0, 'N'},
      {"field-dir",   no_argument,        0, 'f'},
      {"dl-field-maps",no_argument,       0, 'F'},
      {"rand",        required_argument,  0, 'R'},
      {"solenoid-field", required_argument,  0, 'S'},
      {"toroid-field", required_argument,  0, 'T'},
      {"torus-field", required_argument,  0, 'T'},
      {"simple-eg", required_argument,  0, 's'},
      {"simple", required_argument,  0, 's'},
      {0,0,0,0}
   };
   bool should_exit = false;
   while(iarg != -1) {
      iarg = getopt_long(argc, argv, "o:h:g:t:D:r:V:i:R:n:O:S:T:s:bhINfF", longopts, &index);

      switch (iarg)
      {
         case 'i':
            input_file_name = optarg;
            break;

         case 'b':
            is_interactive = false;
            use_gui = false;
            use_vis = false;
            break;

         case 'I':
            is_interactive = true;
            break;

         case 'V':
            if( atoi(optarg) == 0 ){
               use_vis = false;
            } else  {
               use_vis = true;
            }
            break;

         case 'r':
            run_number = atoi( optarg );
            break;

         case 'n':
            number_of_events = atoi( optarg );
            break;

         case 'O':
            event_start_offset = atoi( optarg );
            if(event_start_offset < 0 ) {
               std::cout << "Error: -O, --offset is negative! Aborting.\n";
               std::exit(EXIT_FAILURE);
            }
            break;

         case 'g':
            if( atoi(optarg) == 0 ){
               use_gui = false;
            } else  {
               use_gui = true;
            }
            break;

         case 't':
            output_tree_name = optarg;
            break;

         case 'D':
            output_dir = optarg;
            break;

         case 'N':
            run_manager_init = true;
            break;

         case 'o':
            output_file_name = optarg;
            if( fexists(output_file_name) ) {
               std::cout << "Error : " << output_file_name << " already exist"  << std::endl;
               exit(EXIT_FAILURE);
            }
            break;

         case 'h':
            print_help();
            should_exit = true;
            break;

         case 'f':
            print_field_dir();
            exit(0);
            break;

         case 'F':
            download_field_maps();
            std::quick_exit(0);
            break;

         case 'R':
            set_rand_seed = atoi( optarg );
            break;

         case 'T':
            if( std::string(optarg) == "in-bend") {
               torus_field_sign = -1.0;
            } else if( std::string(optarg) == "out-bend") {
               torus_field_sign = 1.0;
            } else {
               try {
                  //std::cout << " argument " << optarg << std::endl;
                  torus_field_scale = std::stod( optarg );
                  if(torus_field_scale < 0 ) torus_field_sign = -1.0;
                  torus_field_scale = std::abs(torus_field_scale);
               } catch (const std::invalid_argument& ia) {
                  std::cerr << "Invalid argument for torus field : " << optarg << '\n';
                  exit(0);
               }
            }
            break;

         case 'S':
            if( std::string(optarg) == "in-bend") {
               solenoid_field_sign = -1.0;
            } else if( std::string(optarg) == "out-bend") {
               solenoid_field_sign = 1.0;
            } else {
               try {
                  solenoid_field_scale = std::stod( optarg );
                  if(solenoid_field_scale < 0 ) solenoid_field_sign = -1.0;
                  solenoid_field_scale = std::abs(solenoid_field_scale);
               } catch (const std::invalid_argument& ia) {
                  std::cerr << "Invalid argument for solenoid field : " << optarg << '\n';
                  exit(0);
               }
            }
            break;

         case 's':
            simple_eg = atoi( optarg );
            break;

         case '?':
            print_help();
            exit(EXIT_FAILURE);
            break;
      }
   }
   if(should_exit) return 0;


   // here we assume the last argument is a macro file 
   if( optind < argc ) {
      has_macro_file = true;
   }
   for (int i = optind; i < argc; i++) {
      theRest        += argv[i];
   }

   // Get the piped commands
   std::vector<std::string> piped_commands;
   if(!isatty(STDIN_FILENO)) {
      std::cout << "Reading piped commands...\n";
      std::string lineInput;
      while(std::getline(std::cin,lineInput)) {
         piped_commands.push_back(lineInput);
      }
   }

   check_field_maps();

   //std::cout << " the rest of the arguments: " << theRest << std::endl;
   //std::cout << "input  : " << input_file_name << std::endl;
   //std::cout << "output : " << output_file_name << std::endl;
   //std::cout << "  tree : " << output_tree_name << std::endl;
   //std::cout << "  dir  : " << output_dir << std::endl;

   //---------------------------------------------------------------------------

   // Detect interactive mode (if no arguments) and define UI session
   // Note third argument of G4UIExecutive can be ("qt", "xm", "win32", "gag", "tcsh", "csh")
   G4UIExecutive* ui = 0;
   char* ex_name = "c12sim \0";
   if( is_interactive || use_gui ) {
      if( use_gui ) {
         ui = new G4UIExecutive(1, &ex_name, "qt");
      } else {
         ui = new G4UIExecutive(1, &ex_name, "tcsh");
      }
   }

   // Choose the Random engine
   TRandom3 rand(0);
   G4Random::setTheEngine(new CLHEP::RanecuEngine);
   G4long seed = rand.GetSeed();
   if( set_rand_seed != -1 ) {
      seed = long(set_rand_seed);
      rand.SetSeed(seed);
   }
   CLHEP::HepRandom::setTheSeed(seed);

   // Construct the default run manager
//#ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//#else
   G4RunManager* runManager = new G4RunManager;
//#endif

   SimulationManager * simManager = SimulationManager::GetInstance();
   if( input_file_name.size() != 0 ) {
      simManager->SetInputFileName( input_file_name );
   }
   if( output_file_name.size() != 0 ) {
      simManager->SetOutputFileName( output_file_name );
   }
   if( output_tree_name.size() != 0 ) {
      simManager->SetOutputTreeName( output_tree_name );
   }
   if( output_dir.size() != 0 ) {
      simManager->SetOutputDirectoryName( output_dir );
   }

   simManager->SetSolenoidFieldScale( solenoid_field_sign*solenoid_field_scale );
   simManager->SetToroidFieldScale( torus_field_sign*torus_field_scale );

   simManager->SetEventNumber(event_start_offset);


   B1DetectorConstruction      * realWorld     = new B1DetectorConstruction();

   // Detector construction with parallel world 
   // Not using the parallel world because there is a bug in it. 
   // http://bugzilla-geant4.kek.jp/show_bug.cgi?id=1800
   //G4String                      paraWorldName = "ParallelWorld";
   //B1ParallelWorldConstruction * parallelWorld = new B1ParallelWorldConstruction(paraWorldName);
   //realWorld->RegisterParallelWorld(parallelWorld);

   runManager->SetUserInitialization(realWorld);

   // Physics list
   G4PhysListFactory     factory;
   //QGSP_BIC_EMY QGSP_BERT_HP_PEN QGSP_BIC_LIV
   G4VModularPhysicsList * physicsList =  factory.GetReferencePhysList("FTFP_BERT");
   //G4VModularPhysicsList* physicsList = new CLAS12_QGSP_BIC(paraWorldName,false);
   //G4VModularPhysicsList* physicsList = new QBBC;
   physicsList->DumpCutValuesTable();

   // This connects the phyics to the parallel world (and sensitive detectors)
   //physicsList->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName,/*layered_mass=*/true));
   //physicsList->ReplacePhysics(new G4IonQMDPhysics());
   //physicsList->SetDefaultCutValue(0.005*um);

   physicsList->RegisterPhysics(new B1OpticalPhysics());
   // This is needed to make use of the G4UserLimits applied to logical volumes.
   //physicsList->RegisterPhysics(new G4StepLimiterPhysics());

   physicsList->SetVerboseLevel(1);
   runManager->SetUserInitialization(physicsList);

   // User action initialization
   runManager->SetUserInitialization(new B1ActionInitialization(run_number));
  
   // Primary generator
   if( simple_eg == 0 ) {
      runManager->SetUserAction( new B1PrimaryGeneratorAction(event_start_offset) );
   } else {
      runManager->SetUserAction( new SimplePrimaryGeneratorAction(simple_eg) );
   }

   // Initialize G4 kernel
   runManager->Initialize();
   // This constructs  the geometry and many other things. It might be best to leave this 
   // to a UI command...

   // Initialize visualization
   //
   G4VisManager* visManager = nullptr;//new G4VisExecutive;
   if(use_vis) {
      visManager = new G4VisExecutive("Quiet");
      // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
      // G4VisManager* visManager = new G4VisExecutive("Quiet");
      visManager->Initialize();
   }

   // Get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();
   G4UIQt* qui = static_cast<G4UIQt*> (UImanager->GetG4UIWindow());
   if (qui) {
      qui->GetMainWindow()->setVisible(true);
   }

   // Process macro or start UI session

   // Execute macro file if provoded otherwise run the default macro
   if( has_macro_file ) {
      G4String command = "/control/execute ";
      G4String fileName = argv[optind];
      UImanager->ApplyCommand(command+fileName);

   } else {

      G4String command   = "/control/macroPath ";
      G4String mac_dir   = C12SIM_MACRO_DIR;
      G4String fileName = "init_default.mac";

      //std::cout << " executing " << command+mac_dir << std::endl;
      UImanager->ApplyCommand(command+mac_dir);
      command = "/control/execute ";
      if( use_vis ) { 
         //std::cout << " executing " << command+fileName << std::endl;
         UImanager->ApplyCommand(command+fileName);
      }
   }

   // run the piped commands
   for(auto cmd : piped_commands){
      std::cout << cmd << std::endl;
      UImanager->ApplyCommand(G4String(cmd));
   }


   // run the set number of events
   if( number_of_events > 0) {
      G4String command = "/run/beamOn " + std::to_string(number_of_events);
      UImanager->ApplyCommand( command );
   }

   // interactive mode
   if( is_interactive )  {
      ui->SessionStart();
      delete ui;
   }

   if(use_vis) delete visManager;
   delete runManager;
   std::string profiler_log_filename = "data/log/prof_c12sim_";
   profiler_log_filename += std::to_string(run_number);
   profiler_log_filename += ".md";
   std::ofstream prof_log(profiler_log_filename);
   std::cout << " ================================================\n";
   std::cout << " dollar profiling:\n";
   std::cout << " ================================================\n";
   dollar::text(std::cout);
   dollar::markdown(prof_log);
   std::cout << " ================================================\n";

   // Job termination
   // Free the store: user actions, physics_list and detector_description are
   // owned and deleted by the run manager, so they should not be deleted 
   // in the main() program !

}
//______________________________________________________________________________

std::string exec(const char* cmd) {
   std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
   if (!pipe) return "ERROR";
   char buffer[128];
   std::string result = "";
   while (!feof(pipe.get())) {
      if (fgets(buffer, 128, pipe.get()) != NULL)
         result += buffer;
   }
   return result;
}
//______________________________________________________________________________

bool fexists(const std::string& filename) {
   std::ifstream ifile(filename.c_str());
   if( ifile ) return true;
   return false;
}
//______________________________________________________________________________

void print_field_dir()
{
   std::cout << "c12sim field map directory : " C12SIM_DATA_DIR << std::endl;

}
//______________________________________________________________________________

void check_field_maps()
{
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat
   bool failed = false;

   if( ! fexists(C12SIM_DATA_DIR"/clas12SolenoidFieldMap.dat") ) {
      std::cerr << "Error: Solenoid field map missing!" << std::endl;
      std::cerr << "file \"" C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat\" not found. " << std::endl;
      failed = true;
   }
   if( ! fexists(C12SIM_DATA_DIR"/clas12TorusOriginalMap.dat") ) {
      std::cerr << "Error: Torus field map missing!" << std::endl;
      std::cerr << "file \"" C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat\" not found. " << std::endl;
      failed = true;
   }
   if(failed) {
      std::cerr << "Use \"c12sim --dl-field-maps\" to download and install the field maps" << std::endl;
      std::quick_exit(0);
   }
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}
void download_field_maps()
{
   const int dir_err = system("mkdir -p "C12SIM_DATA_DIR);
   if (-1 == dir_err) {
     printf("Error creating directory!n");
     exit(1);
   }
   std::cout << "c12sim field map directory : " C12SIM_DATA_DIR << std::endl;
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat
   bool dl_solenoid = false;
   bool dl_toroid   = false;
   const unsigned int solenoid_file_size = 31037896;
   const unsigned int toroid_file_size   = 31037896;

   //std::cout << " solenoid file size : " << filesize(C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat") << "\n";
   //std::cout << " toroid   file size : " << filesize(C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat") << "\n";

   if ((!fexists(C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat")) ||
       (filesize(C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat")) < solenoid_file_size) {
     dl_solenoid = true;
   }
   if ((!fexists(C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat")) ||
       (filesize(C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat")) < toroid_file_size) {
     dl_toroid = true;
   }

   if (dl_solenoid) {
     std::cout << " Downloading  "
                  "http://clasweb.jlab.org/12gev/field_maps/"
                  "clas12SolenoidFieldMap.dat"
               << std::endl;
     exec(" mkdir -p " C12SIM_DATA_DIR " ; cd " C12SIM_DATA_DIR " ; wget "
          "http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat "
          "-O clas12SolenoidFieldMap.dat ");
   }

   if (dl_toroid) {
     std::cout << " Downloading  "
                  "http://clasweb.jlab.org/12gev/field_maps/"
                  "clas12TorusOriginalMap.dat"
               << std::endl;
     // std::cout << " cd " C12SIM_DATA_DIR "\n";
     exec(" cd " C12SIM_DATA_DIR " ; wget "
          "http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat "
          "-O clas12TorusOriginalMap.dat ");
   }
   std::cout << " solenoid file size : " << filesize(C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat") << "\n";
   std::cout << " toroid   file size : " << filesize(C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat") << "\n";
}
//______________________________________________________________________________

void print_help() {

   std::cout << "usage: c12sim [options] [macro file]    \n";
   std::cout << "Options:                               \n";
   std::cout << "    -r, --run=NUMBER         Set simulation \"run\" NUMBER\n";
   std::cout << "    -i, --input=FILE         Set the input file from which events will be read\n";
   std::cout << "    -o, --output=NAME        Set the output file name which will have the run number appended\n";
   std::cout << "    -n, --events=#           Causes the execution of the ui command \"/run/beamOn #\".\n";
   std::cout << "                             This happens just before exiting (in the case of batch mode) or returning to UI prompt.\n";
   std::cout << "                             Default is \"clas12sim\". Note this is just the file basename; use -D to set the directory \n";
   std::cout << "    -O, --offset=#           Starting event offset.\n";
   std::cout << "    -D, --dir=NAME           Set the output directory. The default is \"data/rootfiles/\"\n";
   std::cout << "    -t, --treename=NAME      Set the output tree name. The default is \"clasdigi_hits\"\n";
   std::cout << "    -g, --gui=#              Set to 1 (default) to use qt gui or\n";
   std::cout << "                             0 to use command line\n";
   std::cout << "    -V, --vis=#              set to 1 (default) to visualization geometry and events\n";
   std::cout << "                             0 to turn off visualization\n";
   std::cout << "    -I, --interactive        run in interactive mode (default)\n";
   std::cout << "    -N, --init               run without initializing G4 kernel\n";
   std::cout << "    -b, --batch              run in batch mode\n"; 
   std::cout << "    -f, --field-dir          print the field dir\n"; 
   std::cout << "    -F, --dl-field-maps      downloads the field maps into the field dir\n"; 
   std::cout << "    -R, --rand=NUMBER        set the random number seed\n";
   std::cout << "    -S, --solenoid-field=B   Scale or turn off solenoid field\n";
   std::cout << "    -T, --toroid-field=B     Scale or turn off solenoid field\n";
   //   {"solenoid-field", required_argument,  0, 'S'},
   //   {"toroid-field", required_argument,  0, 'T'},
   //   {"torus-field", required_argument,  0, 'T'},
}
//______________________________________________________________________________

