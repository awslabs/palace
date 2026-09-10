// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include "utils/interfaceownership.hpp"

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    std::cerr << "ownership.csv queries.csv coordinate_scale\n";
    return 2;
  }
  mfem::Mpi::Init(argc, argv);
  const double scale = std::stod(argv[3]);
  std::ifstream input(argv[2]);
  std::string line;
  if (!std::getline(input, line) || line != "Group,Slot,X,Y,Z")
  {
    return 3;
  }
  std::map<int, std::unique_ptr<palace::InterfaceOwnershipPartition>> groups;
  int count = 0, mismatches = 0;
  while (std::getline(input, line))
  {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int group, expected;
    mfem::Vector point(3);
    if (!(row >> group >> expected >> point(0) >> point(1) >> point(2)))
    {
      return 4;
    }
    auto &partition = groups[group];
    if (!partition)
    {
      partition =
          std::make_unique<palace::InterfaceOwnershipPartition>(argv[1], group, scale);
    }
    point /= scale;
    const int actual = partition->SelectSlot(point);
    if (actual != expected)
    {
      if (++mismatches <= 5)
      {
        std::cerr << "group=" << group << " expected=" << expected << " actual=" << actual
                  << '\n';
      }
    }
    count++;
  }
  std::cout << "queries=" << count << " mismatches=" << mismatches << '\n';
  return count > 0 && mismatches == 0 ? 0 : 1;
}
