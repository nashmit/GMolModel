/** @file
 * This file checks the basic functionality of Context class.
 */

#include <string>
#include <iostream>
#include <sstream>
#include "simmain.hpp"

int main(int argc, char **argv){
    
    // Declarations    

    InitialContext initialContext;

    // Do whatever

    std::cout << "Initial context's atom list:" << std::endl;
    initialContext.PrintAtomList();

    return 0;
}
