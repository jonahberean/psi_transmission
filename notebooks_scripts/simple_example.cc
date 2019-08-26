//c++ code using ROOT for data analysis.
// v 1.0 (12/11/2017), written by Edgard PIERRE, TRIUMF (epierre@triumf.ca).

//compilation:
//gcc `root-config --cflags` simple_example.cc -o simple_example `root-config --libs`

// ./executable file_name  
// ./simple_example /path_to_the_file/run001

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TRint.h"
#include "TNtupleD.h"
#include "TStopwatch.h"

//--------------------------------------------------//

using namespace std;

int main (int argc, char* argv[]) {


    TString trajfile = argv[1];

    if ((argc!=2) || trajfile=="-h" || trajfile=="-H" || trajfile=="--help" || trajfile=="-help"){
        cout << endl;
        cout << "USAGE:" << endl << endl;
        cout << "-------------------------------------------------------------------------------------" << endl;
        cout << "   Enjoy !!!!!!!!!!" << endl;
        cout << "   e.g.: ./simple_example  run001" << endl;
        cout << "   where the file run001 contains the UCN counts and other informations" << endl;
        cout << "-------------------------------------------------------------------------------------" << endl << endl;
        return 0;
    }

   TRint app("plot",NULL,NULL);

    cout.precision(15);
    cout.width(20);

    TFile rootfile("UCN_count.root","RECREATE");
    TNtupleD *ntuple  = new TNtupleD("ntuple", "UCN_count","t:UCN_cnt:SEC"); 

    Double_t Bin_time(1000);
    Double_t TA_input(8);
    Double_t KP_count(1);

    Double_t UCN_cnt(0), useless(0), SEC(0), time(0);

    char line1[256];
    TString line = "";

    int linecntr(0);
    int nb_lines_header(6);
    ifstream in1;
    in1.open(trajfile);

    // Checks to see if the file is opened:
    if (!(in1.is_open())){

        // is_open() returned false and there is a problem:
        cout << "Can't open the file!" << endl;
        return 0;

    } else {

        // While there are still lines in the file, keep reading:
        while (! in1.eof() ) {

            // place the line from myfile into the line variable:
            in1.getline(line1,256);

            // Displays the line we gathered:
            //cout << line1 << endl;

            // Counts lines
            if (line1[0] != '#') linecntr++;

        }// end of while loop over lines
        cout << "The trajectory contains " << linecntr-nb_lines_header << " points." << endl << endl;

        // Closes the stream:
        in1.close();

    }// end if
    const int npts = linecntr;

    // Starts total computation time monitoring
    TStopwatch watch;
    watch.Reset();
    watch.Start();
    in1.open(trajfile);
    for (Int_t ipth=1; ipth<nb_lines_header; ipth++){
        in1.getline(line1,256);
    }

   for (Int_t ipt=1; ipt<linecntr-nb_lines_header; ipt++){
        in1.getline(line1,256);

        while (line1[0] == '#'){// while loop over comment lines
            in1.getline(line1,256);
        }// end of while loop over line1
        sscanf(line1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&time,&UCN_cnt,&useless,&useless,&useless,&SEC,&useless,&useless,&useless,&useless,&useless,&useless,&useless,&useless,&useless,&useless,&useless,&useless);
	//time=ipt*Bin_time/1000.;
	
	ntuple->Fill(time,UCN_cnt,SEC); 

   }// end of for loop over ipt
    in1.close();
    watch.Stop();
    watch.Print();
    rootfile.Write();
// We're done
cout << "Finished" << endl;


   app.Run();

    return 0;
}
