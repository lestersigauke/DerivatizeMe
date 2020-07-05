#ifndef CMCDD_DZME_H
#define CMCDD_DZME_H
//-----------------------------------------------------------------------------
// Copyright (c) 2020 CMCDD Research Group, Lester Sigauke, Kevin Lobb
//                    Chemistry Dept, Rhodes University
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//-----------------------------------------------------------------------------

#include <openbabel/mol.h>

#include <string>


//typedef std::list<OpenBabel::OBMol> MoleculeList;
typedef std::list<OpenBabel::OBMol>::iterator MoleculeIterator;

//string tools for command line input
void split( std::vector<std::string> & theStringVector,
       const  std::string  & theString,
       const  std::string  & theDelimiter);
//creates a vector of words from a string

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
void
split( std::vector<std::string> & theStringVector, 
       const  std::string  & theString,
       const  std::string  & theDelimiter)
{
    size_t  start = 0, end = 0;

    while ( end != std::string::npos)
    {
        end = theString.find( theDelimiter, start);
        theStringVector.push_back( theString.substr( start,
                       (end == std::string::npos) ? std::string::npos : end - start));
        start = (   ( end > (std::string::npos - theDelimiter.size()) )
                  ?  std::string::npos  :  end + theDelimiter.size());
    }
}
#endif