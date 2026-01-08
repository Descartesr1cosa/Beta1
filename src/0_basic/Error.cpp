#include "0_basic/Error.h"

void ERROR::Abort(std::string Message)
{
    std::cout << "Error ! !\t" << Message << "\n"
              << std::flush;
    exit(-1);
}