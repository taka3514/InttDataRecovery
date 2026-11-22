#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include <cdbobjects/CDBTTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace
{

constexpr double BCO_SEC = 106e-9;

struct Key
{
  int layer;
  int phi;
  int sn;

  bool operator<(const Key& other) const
  {
    return std::tie(layer, phi, sn) <
           std::tie(other.layer, other.phi, other.sn);
  }
};

struct Bin
{
  ULong64_t begin;
  ULong64_t end;
  ULong64_t count;

  double phi;
  double rate;

  int fs;
  int fc;
};

struct Range
{
  ULong64_t begin;
  ULong64_t end;
};

double median(std::vector<double> values)
{
  if (values.empty())
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::sort(values.begin(), values.end());

  const size_t n = values.size();

  return (n % 2)
             ? values[n / 2]
             : 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

}  // namespace

// Usage:
// root -l -b -q 'MakeCDB_official.C(79516,0,40000,1.5,1.0)'

void MakeCDB_official(
    int run = 79516,
    int seg_min = 0,
    int seg_max = 40000,
    double phi_size_cut = 1.5,
    double connect_gap_sec = 1.0,
    const char* output_file = "")
{
  if (seg_min < 0 ||
      seg_max < seg_min ||
      phi_size_cut <= 0 ||
      connect_gap_sec < 0)
  {
    std::cerr << "ERROR: invalid arguments" << std::endl;
    return;
  }

  const std::string dir = Form(
      "/sphenix/tg/tg01/commissioning/INTT/work/takakiku/"
      "GitHub_OutPut/ScanTrkrCluster/%d",
      run);

  const std::string output_dir =
      "/sphenix/tg/tg01/commissioning/INTT/work/takakiku/"
      "GitHub_OutPut/WordShiftCDB/Official";

  const std::string out =
      (output_file && *output_file)
          ? output_file
          : Form("%s/calib_%08d.root", output_dir.c_str(), run);

  std::map<Key, std::vector<Bin>> bins;
  std::map<Key, std::pair<int, int>> mapping;

  ULong64_t first_bco =
      std::numeric_limits<ULong64_t>::max();

  ULong64_t last_bco = 0;

  int nfiles = 0;
  int invalid = 0;

  // ==================================================
  // Collect existing ROOT files
  // ==================================================

  std::vector<std::pair<int, std::string>> input_files;

  const std::string prefix = Form("stats_%08d_", run);
  const std::string suffix = ".root";

  if (!fs::exists(dir) || !fs::is_directory(dir))
  {
    std::cerr
        << "ERROR: input directory not found: "
        << dir << std::endl;
    return;
  }

  try
  {
    for (const auto& entry : fs::directory_iterator(dir))
    {
      if (!entry.is_regular_file())
      {
        continue;
      }

      const std::string name =
          entry.path().filename().string();

      if (name.size() != prefix.size() + 5 + suffix.size())
      {
        continue;
      }

      if (name.compare(0, prefix.size(), prefix) != 0 ||
          name.compare(
              name.size() - suffix.size(),
              suffix.size(),
              suffix) != 0)
      {
        continue;
      }

      const std::string seg_str =
          name.substr(prefix.size(), 5);

      if (!std::all_of(
              seg_str.begin(),
              seg_str.end(),
              [](unsigned char c)
              {
                return c >= '0' && c <= '9';
              }))
      {
        continue;
      }

      const int seg = std::stoi(seg_str);

      if (seg < seg_min || seg > seg_max)
      {
        continue;
      }

      input_files.emplace_back(
          seg,
          entry.path().string());
    }
  }
  catch (const fs::filesystem_error& e)
  {
    std::cerr
        << "ERROR: directory scan failed: "
        << e.what() << std::endl;
    return;
  }

  std::sort(
      input_files.begin(),
      input_files.end());

  std::cout
      << "Found " << input_files.size()
      << " ROOT files"
      << " | run=" << run
      << " | segment range="
      << seg_min << "-" << seg_max
      << std::endl;

  if (input_files.empty())
  {
    std::cerr
        << "ERROR: no input ROOT files found in "
        << dir << std::endl;
    return;
  }

  // ==================================================
  // Read segment statistics
  // ==================================================

  for (const auto& [seg, filename] : input_files)
  {
    std::unique_ptr<TFile> file(
        TFile::Open(filename.c_str(), "READ"));

    if (!file || file->IsZombie())
    {
      ++invalid;
      continue;
    }

    TTree* tree = nullptr;
    file->GetObject("Single", tree);

    if (!tree)
    {
      ++invalid;

      std::cerr << "Missing Single: "
                << filename << std::endl;

      continue;
    }

    const char* required[] = {
        "layer",
        "phi",
        "sn",
        "felix_server",
        "felix_channel",
        "begin_bco",
        "end_bco",
        "cluster_count",
        "mean_phi_size",
        "cluster_rate_tb0",
        "hit_rate_tb0"};

    bool valid = true;

    for (const char* name : required)
    {
      if (!tree->GetBranch(name))
      {
        std::cerr << "Missing "
                  << name
                  << " in "
                  << filename
                  << std::endl;

        valid = false;
      }
    }

    if (!valid)
    {
      ++invalid;
      continue;
    }

    int layer = 0;
    int phi = 0;
    int sn = 0;

    int fs = 0;
    int fc = 0;

    ULong64_t begin = 0;
    ULong64_t end = 0;
    ULong64_t count = 0;

    double mean = 0;
    double rate = 0;
    double hit_rate = 0;

    tree->SetBranchAddress("layer", &layer);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("sn", &sn);

    tree->SetBranchAddress("felix_server", &fs);
    tree->SetBranchAddress("felix_channel", &fc);

    tree->SetBranchAddress("begin_bco", &begin);
    tree->SetBranchAddress("end_bco", &end);

    tree->SetBranchAddress("cluster_count", &count);
    tree->SetBranchAddress("mean_phi_size", &mean);

    tree->SetBranchAddress("cluster_rate_tb0", &rate);
    tree->SetBranchAddress("hit_rate_tb0", &hit_rate);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i)
    {
      if (tree->GetEntry(i) <= 0)
      {
        continue;
      }

      if (layer < 3 ||
          layer > 6 ||
          phi < 0 ||
          phi >= ((layer < 5) ? 12 : 16) ||
          sn < 0 ||
          sn > 1 ||
          end <= begin)
      {
        continue;
      }

      const Key key{layer, phi, sn};

      const auto mapping_it = mapping.find(key);

      if (mapping_it != mapping.end() &&
          mapping_it->second != std::make_pair(fs, fc))
      {
        std::cerr
            << "ERROR: inconsistent FELIX mapping in "
            << filename << std::endl;

        return;
      }

      mapping[key] = {fs, fc};

      bins[key].push_back({
          begin,
          end,
          count,
          mean,
          rate,
          fs,
          fc});

      first_bco = std::min(first_bco, begin);
      last_bco = std::max(last_bco, end);
    }

    ++nfiles;

    if (nfiles % 100 == 0)
    {
      std::cout
          << "Processed " << nfiles
          << " / " << input_files.size()
          << " files (current segment = "
          << seg << ")"
          << std::endl;
    }
  }

  if (!nfiles ||
      first_bco == std::numeric_limits<ULong64_t>::max())
  {
    std::cerr
        << "ERROR: no valid stats input in "
        << dir << std::endl;

    return;
  }

  // ==================================================
  // Prepare CDB output
  // ==================================================

  if (gSystem->mkdir(output_dir.c_str(), true) != 0 &&
      gSystem->AccessPathName(output_dir.c_str()))
  {
    std::cerr
        << "ERROR: cannot create output directory "
        << output_dir << std::endl;

    return;
  }

  CDBTTree cdb(out);

  int period_id = 0;
  long long nbad = 0;

  const ULong64_t gap_bco =
      static_cast<ULong64_t>(
          std::llround(connect_gap_sec / BCO_SEC));

  // ==================================================
  // Analyze each half-ladder
  // ==================================================

  for (int layer = 3; layer <= 6; ++layer)
  {
    for (int phi = 0;
         phi < ((layer < 5) ? 12 : 16);
         ++phi)
    {
      for (int sn = 0; sn < 2; ++sn)
      {
        const Key key{layer, phi, sn};

        auto bins_it = bins.find(key);

        if (bins_it == bins.end() ||
            bins_it->second.empty())
        {
          continue;
        }

        auto& values = bins_it->second;

        std::sort(
            values.begin(),
            values.end(),
            [](const Bin& a, const Bin& b)
            {
              return std::tie(a.begin, a.end) <
                     std::tie(b.begin, b.end);
            });

        const auto mapping_it = mapping.find(key);

        if (mapping_it == mapping.end())
        {
          std::cerr
              << "ERROR: missing FELIX mapping"
              << " layer=" << layer
              << " phi=" << phi
              << " sn=" << sn
              << std::endl;

          return;
        }

        const int fs = mapping_it->second.first;
        const int fc = mapping_it->second.second;

        if (!((fs >= 0 && fs <= 7) ||
              (fs >= 3001 && fs <= 3008)) ||
            fc < 0 || fc > 13)
        {
          std::cerr
              << "ERROR: invalid FELIX mapping"
              << " fs=" << fs
              << " fc=" << fc
              << std::endl;

          return;
        }

        std::vector<double> good_rates;
        std::vector<double> bad_rates;
        std::vector<Range> periods;

        bool active = false;

        ULong64_t bad_start = 0;
        ULong64_t bad_end = 0;

        // ------------------------------------------
        // Find bad periods
        // ------------------------------------------

        for (const auto& bin : values)
        {
          if (bin.count == 0 ||
              !std::isfinite(bin.phi) ||
              bin.phi <= 0)
          {
            continue;
          }

          if (bin.rate > 0 &&
              std::isfinite(bin.rate))
          {
            if (bin.phi < phi_size_cut)
            {
              bad_rates.push_back(bin.rate);
            }
            else
            {
              good_rates.push_back(bin.rate);
            }
          }

          if (bin.phi < phi_size_cut)
          {
            if (!active)
            {
              bad_start = bin.begin;
              bad_end = bin.end;
              active = true;
            }
            else if (
                bin.begin <= bad_end ||
                bin.begin - bad_end <= gap_bco)
            {
              bad_end = std::max(bad_end, bin.end);
            }
            else
            {
              periods.push_back({
                  bad_start,
                  bad_end});

              bad_start = bin.begin;
              bad_end = bin.end;
            }
          }
          else if (active)
          {
            periods.push_back({
                bad_start,
                bad_end});

            active = false;
          }
        }

        if (active)
        {
          periods.push_back({
              bad_start,
              bad_end});
        }

        // ------------------------------------------
        // Determine failure type
        // ------------------------------------------

        const double normal_center = median(good_rates);
        const double bad_center = median(bad_rates);

        double hit_rate_ratio =
            std::numeric_limits<double>::quiet_NaN();

        int type = -1;

        if (std::isfinite(normal_center) &&
            normal_center > 0 &&
            std::isfinite(bad_center))
        {
          hit_rate_ratio =
              bad_center / normal_center;

          if (bad_center > normal_center)
          {
            type = 0;  // Enhancement
          }
          else if (bad_center < normal_center)
          {
            type = 1;  // Dropping
          }
          else
          {
            type = -1;  // Equal / unknown
          }
        }

        // ------------------------------------------
        // Store one bad period per CDB entry
        // ------------------------------------------

        for (const auto& period : periods)
        {
          const int id = period_id++;

          cdb.SetIntValue(
              id, "felix_server", fs);

          cdb.SetIntValue(
              id, "felix_channel", fc);

          cdb.SetIntValue(
              id, "layer", layer);

          cdb.SetIntValue(
              id, "phi", phi);

          cdb.SetIntValue(
              id, "sn", sn);

          cdb.SetUInt64Value(
              id, "bad_begin", period.begin);

          cdb.SetUInt64Value(
              id, "bad_end", period.end);

          cdb.SetIntValue(
              id, "bad_type", type);

          cdb.SetDoubleValue(
              id, "normal_center", normal_center);

          cdb.SetDoubleValue(
              id, "bad_center", bad_center);

          cdb.SetDoubleValue(
              id, "hit_rate_ratio", hit_rate_ratio);

          ++nbad;
        }
      }
    }
  }

  // ==================================================
  // Write CDBTTree
  // ==================================================

  if (period_id == 0)
  {
    std::cerr
        << "ERROR: no bad periods found; "
        << "CDB file was not written"
        << std::endl;

    return;
  }

  cdb.Commit();
  cdb.WriteCDBTTree();

  std::cout
      << "Saved CDB: " << out
      << " | files=" << nfiles
      << " | invalid=" << invalid
      << " | bad period entries=" << period_id
      << " | bad periods=" << nbad
      << " | first BCO=" << first_bco
      << " | last BCO=" << last_bco
      << std::endl;
}