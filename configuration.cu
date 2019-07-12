/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include <fstream>
#include <string>
#include "configuration.h"
using namespace std;


__host__ void Configuration::read(char *path)
{   
    // Open configuration file.
    ifstream config_file;
    config_file.open(path);

    // Ignore header line.
    string header;
    getline(config_file, header);

    // Read in values.
    config_file >> Lx >> Ly >> dt >> nt >> nprey >> mass >> kappa >> g >> catt >> latt >> crep >> lrep >> gamma >> cavoid >> omega >> R;

    // Close the file.
    config_file.close();

    return;
}
