#pragma once

#include <fstream>
#include <vector>
#include <string>

#include "bgmg_log.h"

// reader must know the type

template<typename T>
void save_value(std::ofstream& os, T value) {
  os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<>
inline void save_value<std::string>(std::ofstream& os, std::string str) {
  size_t numel = str.size();
  os.write(reinterpret_cast<const char*>(&numel), sizeof(size_t));
  os.write(reinterpret_cast<const char*>(&str[0]), numel * sizeof(char));
}

template<typename T>
void load_value(std::ifstream& is, T* value) {
  is.read(reinterpret_cast<char*>(value), sizeof(T));
}

template<>
inline void load_value<std::string>(std::ifstream& is, std::string* str) {
  size_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(size_t));
  str->resize(numel);
  is.read(reinterpret_cast<char*>(&(*str)[0]), numel * sizeof(char));
}

template<typename T>
void load_vector(std::ifstream& is, std::vector<T>* vec) {
  size_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(size_t));
  vec->resize(numel);
  is.read(reinterpret_cast<char*>(&(*vec)[0]), numel * sizeof(T));
}

template<>
inline void load_vector<std::string>(std::ifstream& is, std::vector<std::string>* vec) {
  size_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(size_t));
  vec->resize(numel);
  for (int i = 0; i < numel; i++) load_value(is, &(*vec)[i]);
}

template<typename T>
void save_vector(std::ofstream& os, const std::vector<T>& vec) {
  size_t numel = vec.size();
  os.write(reinterpret_cast<const char*>(&numel), sizeof(size_t));
  os.write(reinterpret_cast<const char*>(&vec[0]), numel * sizeof(T));
}

template<>
inline void save_vector<std::string>(std::ofstream& os, const std::vector<std::string>& vec) {
  size_t numel = vec.size();
  os.write(reinterpret_cast<const char*>(&numel), sizeof(size_t));
  for (int i = 0; i < numel; i++) save_value(os, vec[i]);
}

enum SerializeDirection {
  SerializeDirection_Save = 0,
  SerializeDirection_Load = 1,
};

class BgmgSerializer {
 public:
  BgmgSerializer(std::string filename, SerializeDirection direction) : filename_(filename), direction_(direction) {
    if (direction_==SerializeDirection_Save) {
      os_.open(filename, std::ofstream::binary);
      if (!os_) BGMG_THROW_EXCEPTION(std::runtime_error(::std::runtime_error("can't open" + filename)));
    }

    if (direction_==SerializeDirection_Load) {
      is_.open(filename, std::ifstream::binary);
      if (!is_) BGMG_THROW_EXCEPTION(::std::runtime_error("can't open" + filename));
    }
  }

  ~BgmgSerializer() {
    if (direction_==SerializeDirection_Save) os_.close();
    if (direction_==SerializeDirection_Load) is_.close();
  }

  bool is_save() const { return direction_==SerializeDirection_Save; }
  bool is_load() const { return direction_==SerializeDirection_Load; }

template<typename T>
void dump_value(T* value) {
  if (direction_==SerializeDirection_Save) {
    save_value<T>(os_, *value);
  }

  if (direction_==SerializeDirection_Load) {
    load_value<T>(is_, value);
    if (!is_) BGMG_THROW_EXCEPTION(::std::runtime_error("can't read from " + filename_));
  }
}

template<typename T>
void dump_vector(std::vector<T>* vec) {
  if (direction_==SerializeDirection_Save) {
    save_vector<T>(os_, *vec);
  }
  
  if (direction_==SerializeDirection_Load) {
    load_vector<T>(is_, vec);
    if (!is_) BGMG_THROW_EXCEPTION(::std::runtime_error("can't read from " + filename_));
  }
}

private:
  const std::string filename_;
  const SerializeDirection direction_;
  std::ofstream os_;
  std::ifstream is_;
};
