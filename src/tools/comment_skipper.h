/*************************************************************************
 * comment_skipper.h - skips all comments from the input.                *
 *                                                                       *
 * Copyright (C) 2014-2015 by Laboratory of Medicinal Chemistry of       *
 * Department of Chemistry, Lomonosov Moscow State University            *
 *                                                                       *
 * Programmed by Oleg. I. Titov                                          *
 *                                                                       *
 * This file is part of the Electrostatic Tools project.                 *
 * For more information, see <http://qsar.chem.msu.ru/elec_tools/>       *
 *                                                                       *
 * Based on work:                                                        *
 * O.I. Titov, D.A. Shulga, V.A. Palyulin, N.S. Zefirov                  *
 * Description of halogen bonding on the basis of multicenter multipole  *
 * expansion // Doklady Chemistry, 2013, vol. 450, p. 139                *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation version 2 of the License.                *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 ************************************************************************/

#ifndef COMMENT_SKIPPER_HEADER
#define COMMENT_SKIPPER_HEADER

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <ios>

// warning!! reads all stream contents to buffer

class CommentSkipper
{
public:
  CommentSkipper(std::string open = std::string("/*"), std::string close = std::string("*/"), std::string line = ("//"), bool skip_blank_lines = false)
      : comm_open(open), comm_close(close), line_begin(line), blank_lines(skip_blank_lines) { }

  const std::string& skip(const std::string& text)
  {
    // removing c++-style comments
    std::ostringstream os;
    std::istringstream strm(text);
    std::string s;
    while (std::getline(strm, s))
    {
      std::string::iterator found = std::search(s.begin(), s.end(), line_begin.begin(), line_begin.end());
      os << std::string(s.begin(), found) << std::endl; 
    }
        
    // removing c-style comments
    std::string text_ = os.str();
    std::string::const_iterator first = text_.begin(), last = text_.end();
    data = std::string();
    while (first != last)
    {
      std::string::const_iterator found = std::search(first, last, comm_open.begin(), comm_open.end());
      data += std::string(first, found);
      found = std::search(found, last, comm_close.begin(), comm_close.end());
      for (int i = 0; i < comm_close.length(); ++i) if (found != last) found++; else break;
      first = found;
    }
    
    // removing blank lines
    if (blank_lines)
    {
      std::ostringstream ostrm;
      std::istringstream strm(data);
      std::string l;
      while (strm.good())
      {
        std::getline(strm, l);
        if (l.find_first_of(" \t\n\r") != std::string::npos) ostrm << l << std::endl;
      }
      data = ostrm.str();
    }
    
    return data;
  }

  const std::string& skip(std::istream& s)
  {
    std::string buffer, line;
    s >> std::noskipws;
    while (s.good())
    {
      std::getline(s, line);
      buffer += line + "\n";
    }
    return skip(buffer);
  }

  const std::string& skipped() { return data; }

private:
  std::string data, comm_open, comm_close, line_begin;
  bool blank_lines;
};

#endif // COMMENT_SKIPPER_HEADER
