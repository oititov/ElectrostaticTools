#ifndef PAR_TABLE_H
#define PAR_TABLE_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <iomanip>

//
//  Reads parameters from a table of the format
//
//                  header1   header2    header3   ....
//  type1           par       par        fixed_par*
//  type2           par       par        par_as_above-
//  fixed_type      par       par        par            *
//  type_as_above   par       par        par            -
//  ....
// 
//

class ParameterTable
{
public:
  ParameterTable(std::istream& strm) : cell_width(15)
  {
    readMe(strm); 
  }
  bool exists(const std::string& row, const std::string& col) const
  {
    return par_map.find(make_pair(row,col)) != par_map.end();
  }
  bool fixed(const std::string& row, const std::string& col) const 
  { 
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() ) return unsigned(i->second) >= fixed_begin; else throw std::domain_error("Unknown parameters requested from parameter table");    
  }
  double get(const std::string& row, const std::string& col) 
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() ) 
    {
      used_parameters[i->second]++;
      return raw_data[i->second]; 
    } else {
      std::map<std::string, double>::const_iterator j = col_defaults.find(col);
      if ( j != col_defaults.end() ) return j->second;
      else throw std::domain_error("Unknown parameters requested from parameter table");
    }
  }
  double getUnmarking(const std::string& row, const std::string& col) const
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() ) 
    {
      return raw_data[i->second]; 
    } else {
      std::map<std::string, double>::const_iterator j = col_defaults.find(col);
      if ( j != col_defaults.end() ) return j->second;  
      else throw std::domain_error("Unknown parameters requested from parameter table");
    }
  }
  void mark(const std::string& row, const std::string& col)
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() )
    {
      used_parameters[i->second]++;
    } else throw std::domain_error("Unknown parameters requested from parameter table");
  }
  double* pGet(const std::string& row, const std::string& col)
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() ) 
    {
      used_parameters[i->second]++;
      return &(raw_data[i->second]); 
    } else throw std::domain_error("Unknown parameters requested from parameter table");
  }
  int getUseCount(const std::string& row, const std::string& col) const
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() )
    {
      return used_parameters[i->second];
    } else throw std::domain_error("Unknown parameters requested from parameter table");
  }
  void set(const std::string& row, const std::string& col, const double value) 
  {
    std::map<Sp, int>::iterator i = par_map.find(make_pair(row,col));
    if ( i == par_map.end() ) throw std::domain_error("Unknown parameters requested from parameter table");
    if ( unsigned(i->second) < fixed_begin ) raw_data[i->second] = value;
  }
  void setDefault(const std::string& col, const double val)
  {
    col_defaults[col] = val;
  }
  void unsetDefault(const std::string& col)
  {
    col_defaults.erase(col);
  }
  void markUnused(const std::string& row, const std::string& col)
  {
    std::map<Sp, int>::const_iterator i = par_map.find(make_pair(row,col));
    if ( i != par_map.end() ) used_parameters[i->second] = 0;
  }
  std::vector<double> getActive() const 
  { 
    std::vector<double> result;
    for (unsigned int i = 0; i < fixed_begin; ++i) if (used_parameters[i]) result.push_back(raw_data[i]);
    return result;     
  }
  int countActive() const { int result = 0;  for (unsigned int i = 0; i < fixed_begin; ++i) if (used_parameters[i]) result++; return result; }
  void setActive(const std::vector<double>& v, const int n = -1)
  { 
    unsigned int l;
    if (n < 0) l = v.size(); else l = n;
    if ( l > fixed_begin || l > v.size()) throw std::domain_error("Wrong length vector specified to Parameter Table");
    for (unsigned int i = 0, j = 0; i < fixed_begin && j < l; ++i)  if (used_parameters[i])  { raw_data[i] = v[j]; j++; }
  }
  void setCellWidth(unsigned int w) { cell_width = w; }
  void readMe(std::istream& strm) 
  {
    std::string line;
    // reading header
    do std::getline(strm, line); while (line.find_first_not_of(" \t\n\r") == std::string::npos && strm.good()) ;
    std::istringstream ss(line);
    std::string t;
    while (ss.good()) 
    {
      t = "";
      ss >> t; 
      if (!t.empty())
      {
        char f = ' ';
        if (t.substr(t.size()-1, t.size()) == "*") 
        { 
          f = '*';
          t = t.substr(0, t.size()-1);
        }
        col_flags.push_back(f);
        headers.push_back(t); 
      }
    }
    //reading table
    
    std::vector<std::pair<double, char> > pairs; 
    while ( strm.good() )
    {
      // read the line
      std::getline(strm, line);
      if (line.find_first_not_of(" \t\n\r") == std::string::npos) continue;
      std::istringstream ss(line);
      std::string type;
      ss >> type;
      rows.push_back(type);
      char flag = ' ';
      for ( unsigned int i = 0; i < headers.size(); ++i)
      {
        ss >> t;
        std::istringstream tss(t);
        double v; char c = ' ';
        tss >> v; if (tss.good()) tss >> c;
        pairs.push_back(std::make_pair(v, c));
        cell_flags.push_back(c);
      }
      if (ss.good()) ss >> flag;
      row_flags.push_back(flag);
    }
    
    // processing the table
    if (row_flags[0] == '-') throw std::domain_error("Reading parameter table: '-' row flag specified in first row.");
    
    std::vector<double> tmp_fixed; // put all fixed values here and record minus position in par_map for further rearrangement
    tmp_fixed.push_back(0.0); // shift so the first real index is 1
    for (unsigned int i = 0; i < rows.size(); ++i)
    {
      for (unsigned int j = 0; j < headers.size(); ++j)
      {
        std::pair<double, char> pair = pairs[i*headers.size()+j];
        switch (col_flags[j])
        {
          case ' ':
            switch (row_flags[i])
            {
              case ' ':
                switch (pair.second)
                {
                  case ' ':
                    par_map[std::make_pair(rows[i], headers[j])] = raw_data.size();
                    raw_data.push_back(pair.first);
                    break;
                  case '-':
                    if (i == 0) throw std::domain_error("Reading parameter table: '-' cell flag specified in first row.");
                    par_map[std::make_pair(rows[i], headers[j])] = par_map[std::make_pair(rows[i-1], headers[j])];
                    break;
                  case '*':
                    par_map[std::make_pair(rows[i], headers[j])] = -tmp_fixed.size();
                    tmp_fixed.push_back( pair.first );
                    break;
                  default:
                    throw std::domain_error("Reading parameter table: unknow cell flag specified.");  
                }
                break;
              case '-':
                par_map[std::make_pair(rows[i], headers[j])] = par_map[std::make_pair(rows[i-1], headers[j])];
                break;
              case '*':
                par_map[std::make_pair(rows[i], headers[j])] = -tmp_fixed.size();
                tmp_fixed.push_back( pair.first );
                break;
              default:
                throw std::domain_error("Reading parameter table: unknow row flag specified.");            
            }
            break;
          case '*':
            par_map[std::make_pair(rows[i], headers[j])] = -tmp_fixed.size();
            tmp_fixed.push_back( pair.first );
            break;
          default:
            throw std::domain_error("Reading parameter table: unknow col flag specified.");
        }
      }
    }
    // final rearrangement
    fixed_begin = raw_data.size();
    raw_data.reserve(raw_data.size() + tmp_fixed.size() - 1);
    for (unsigned int i = 1; i < tmp_fixed.size(); ++i) raw_data.push_back(tmp_fixed[i]);
    for (std::map<Sp, int>::iterator i = par_map.begin(), ie = par_map.end(); i != ie; ++i)
      if (i->second < 0) i->second = fixed_begin - i->second - 1;
    used_parameters.resize(raw_data.size(), 0);
  }
  void printMe(std::ostream& strm) const
  {
    unsigned int width = 0; for (unsigned int i = 0; i < rows.size(); ++i) width = (width < rows[i].size()) ? rows[i].size() : width;
    strm << std::setw(width) << ' ';
    for (unsigned int i = 0; i < headers.size(); ++i) strm << " " << std::setw(cell_width) << headers[i] << col_flags[i];
    strm << std::endl;
    std::vector<char>::const_iterator cell_flag = cell_flags.begin();
    std::vector<char>::const_iterator row_flag = row_flags.begin();
    
    for (unsigned int i = 0; i < rows.size(); ++i)
    {
      strm << rows[i];
      for (unsigned int j = 0; j < headers.size(); ++j) strm << " " << std::setw(cell_width) << this->getUnmarking(rows[i], headers[j]) << (*(cell_flag++));
      strm << " " << (*(row_flag++)) << std::endl;
    }
  }
  void printUseCount(std::ostream& strm) const
  {
    strm << "Parameter use counters: " << std::endl;
    unsigned int width = 0; for (unsigned int i = 0; i < rows.size(); ++i) width = (width < rows[i].size()) ? rows[i].size() : width;
    strm << std::setw(width) << ' ';
    for (unsigned int i = 0; i < headers.size(); ++i) strm << " " << std::setw(cell_width) << headers[i];
    strm << std::endl;
    std::vector<char>::const_iterator cell_flag = cell_flags.begin();
    std::vector<char>::const_iterator row_flag = row_flags.begin();

    for (unsigned int i = 0; i < rows.size(); ++i)
    {
      strm << rows[i];
      for (unsigned int j = 0; j < headers.size(); ++j) strm << " " << std::setw(cell_width) 
                                                             << this->getUseCount(rows[i],headers[j]);
      strm << std::endl;
    }
  }
private:
  std::vector<std::string> rows; // to store type ordering
  std::vector<std::string> headers; // to store column ordering
  std::vector<char> cell_flags; // keep for output
  std::vector<char> row_flags; // keep for output
  std::vector<char> col_flags; // keep for output
  typedef std::pair<std::string, std::string> Sp; // type - header pairs
  std::map<Sp, int> par_map; 
  std::vector<double> raw_data;
  std::vector<int> used_parameters;
  unsigned int fixed_begin;  
  unsigned int cell_width;

  std::map<std::string, double> col_defaults;
};

inline std::ostream& operator << (std::ostream& strm, const ParameterTable& table)
{
  table.printMe(strm); 
  return strm;
}

#endif // PAR_TABLE_H
