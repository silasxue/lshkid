/* 
    Copyright (C) 2008 Wei Dong <wdong@princeton.edu>. All Rights Reserved.
  
    This file is part of LSHKIT.
  
    LSHKIT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LSHKIT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LSHKIT.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * \file matrix-io.h
 * \brief Matrix I/O routines, inlined implementation (Matrix::load & Matrix::save).
 */

#ifndef LSHKIT_MATRIX_IO
#define LSHKIT_MATRIX_IO

#include <lshkit/common.h>
#include <cassert>
#include <fstream>

#ifdef MATRIX_MMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#endif

namespace lshkit{

template <typename T>
void Matrix<T>::peek (const std::string &path, int *elem_size, int *size, int *dim)
{
    unsigned header[3]; /* entry size, row, col */
    assert(sizeof header == 3*4);
    std::ifstream is(path.c_str(), std::ios::binary);
    verify(is);
    is.read((char *)header, sizeof header);
    verify(is);
    *elem_size = header[0];
    *size = header[1];
    *dim = header[2];
}

template <class T>
void Matrix<T>::load (std::istream &is)
{
    unsigned header[3]; /* entry size, row, col */
    assert(sizeof header == 3*4);
    is.read((char *)header, sizeof header);
    verify(is);
    verify(header[0] == sizeof(T));
    reset(header[2], header[1]);
    unsigned sz = sizeof(T) * dim * N;
    is.read((char *)dims, sz);
    verify(is);
}

template <class T>
void Matrix<T>::save (std::ostream &os)
{
    unsigned header[3];
    header[0] = sizeof(T);
    header[1] = N;
    header[2] = dim;
    os.write((char *)header, sizeof header);
    verify(os);
    unsigned sz = sizeof(T) * dim * N;
    os.write((char *)dims, sz);
    verify(os);
}

template <class T>
void Matrix<T>::load (const std::string &path)
{
    std::ifstream is(path.c_str(), std::ios::binary);
    load(is);
}

template <class T>
void Matrix<T>::save (const std::string &path)
{
    std::ofstream os(path.c_str(), std::ios::binary);
    save(os);
}

template <class T>
void Matrix<T>::loadMeta (std::istream &is)
{
	unsigned header[3]; /* entry size, row, col */
	assert(sizeof header == 3*4);
	is.read((char *)header, sizeof header);
	verify(is);
	verify(header[0] == sizeof(T));
	reset(header[2], header[1]);
	//unsigned sz = sizeof(T) * dim * N;
	//is.read((char *)dims, sz);
	verify(is);
}

template <class T>
void Matrix<T>::readIth(std::istream &is, int i)
{
	int k = sizeof(T) * dim;
	is.seekg(3*4+k*i, std::ios_base::beg);
	is.read((char*)dims, k);
}

template <class T>
void Matrix<T>::loadMeta(FILE* fdata)
{
	unsigned header[3]; /* entry size, row, col */
	assert(sizeof header == 3*4);
	fread(header, sizeof(unsigned), 3, fdata);
	verify(header[0] == sizeof(T));
	reset(header[2], header[1]);
}

template <class T>
void Matrix<T>::readIth(FILE* fdata, int i)
{
	int k = sizeof(T) * dim;
	fseek(fdata, 3*4+k*i, SEEK_SET);
	fread(dims, sizeof(T), dim, fdata);
}


#ifdef MATRIX_MMAP
template <class T>
void Matrix<T>::map (const std::string &path) {

    if (dims != NULL) delete[] dims;

    int size = 0;

    fd = open(path.c_str(), O_RDONLY);
    verify(fd >= 0);

    if (read(fd, &size, 4) != sizeof(size)) verify(0);
    if (read(fd, &N, 4) != sizeof(N)) verify(0);
    if (read(fd, &dim, 4) != sizeof(dim)) verify(0);

    verify(size == sizeof(T));

    size_t sz = sizeof(T) * size_t(dim) * size_t(N);
    
    char *start = (char *)mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (start == MAP_FAILED) { verify(0); }
    dims = (T *)(start + 3 * 4);

}

template <class T>
void Matrix<T>::unmap () {
    char *start = (char *)dims - 3 * 4;
    if (munmap(start, sizeof(T) * size_t(dim) * size_t(N)) != 0) verify(0);
    close(fd);

    dims = 0;
    dim = 0;
    N = 0;
}

#endif

}

#endif

