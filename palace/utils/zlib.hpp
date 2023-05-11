// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_ZLIB_HPP
#define PALACE_UTILS_ZLIB_HPP

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <mfem.hpp>

#if defined(MFEM_USE_ZLIB)
#include <zlib.h>
#endif

namespace palace::utils
{

//
// String compression using zlib (https://panthema.net/2007/0328-ZLibString.html).
//

// Compress a STL string using zlib with given compression level and return the binary
// data.
std::string CompressString(const std::string &str, int level = Z_BEST_COMPRESSION)
{
#if defined(MFEM_USE_ZLIB)
  z_stream zs;
  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, level) != Z_OK)
  {
    throw(std::runtime_error("deflateInit failed while compressing."));
  }

  zs.next_in = (Bytef *)str.data();
  zs.avail_in = str.size();  // Set the z_stream's input

  int ret;
  char outbuffer[32768];
  std::string outstring;

  // Retrieve the compressed bytes blockwise.
  do
  {
    zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = deflate(&zs, Z_FINISH);

    if (outstring.size() < zs.total_out)
    {
      // Append the block to the output string.
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  } while (ret == Z_OK);

  deflateEnd(&zs);

  if (ret != Z_STREAM_END)  // An error occurred that was not EOF
  {
    std::ostringstream oss;
    oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }
  return outstring;
#else
  return str;
#endif
}

// Decompress an STL string using zlib and return the original data.
std::string DecompressString(const std::string &str)
{
#if defined(MFEM_USE_ZLIB)
  z_stream zs;
  memset(&zs, 0, sizeof(zs));

  if (inflateInit(&zs) != Z_OK)
  {
    throw(std::runtime_error("inflateInit failed while decompressing."));
  }

  zs.next_in = (Bytef *)str.data();
  zs.avail_in = str.size();

  int ret;
  char outbuffer[32768];
  std::string outstring;

  // Get the decompressed bytes blockwise using repeated calls to inflate.
  do
  {
    zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = inflate(&zs, 0);

    if (outstring.size() < zs.total_out)
    {
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  } while (ret == Z_OK);

  inflateEnd(&zs);

  if (ret != Z_STREAM_END)  // An error occurred that was not EOF
  {
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }
  return outstring;
#else
  return str;
#endif
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_ZLIB_HPP
