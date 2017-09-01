#include "project1_problems.h"
#include <iostream>
using namespace std;
int main(int argc, char* argv[])
{
    char *outfilename;
    int n;
    if( argc<=2){
        cout <<"Error: provide output file name and matrix dimension on same line"<<endl;
        exit(1);
    }
    else{
        outfilename = argv[1];
        n = atoi(argv[2]);
    }
    //problem_1b(n,outfilename);
    //problem_1cd(n,outfilename);
    problem_1e(n);
}
