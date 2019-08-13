/* srepFolder: The directory where s-reps stored.
 * command: ./tps_srep template.vtk target.vtk template.m3d target_srep_folder
 * Liyun Tu
 * Mar 20, 2014
*/

#include "thinplatesplinepdmtosrep.h"



int main( int argc, char* argv[] ){

    if( argc != 5 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <template PDM, target PDM, template s-rep and target s-rep folders>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

//    const char* srepFolder = argv[1];

    thinplatesplinepdmtosrep obj;
    obj.tps_to_srep(argv[1], argv[2], argv[3], argv[4]);


    return 0;
}


