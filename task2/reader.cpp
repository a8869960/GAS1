#include "reader.h"

reader_t::reader_t (const char *filename)
{
  std::ifstream file;
  file.open (filename);
  auto find_word = [this] (const std::string word)
  {
      std::string tmp;
      while (file >> tmp)
        {
          if (word == tmp)
              return true;
        }
      return false;
  };


}
