#include "InttWordShift.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

bool InttWordShift::LoadCDB(const std::string& filename)
{
  m_config.clear();

  if (filename.empty())
  {
    return true;
  }

  std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));

  if (!file || file->IsZombie())
  {
    std::cerr << "InttWordShift: cannot open " << filename << std::endl;
    return false;
  }

  TTree* tree = nullptr;
  file->GetObject("Multiple", tree);

  if (!tree)
  {
    std::cerr << "InttWordShift: Multiple tree missing" << std::endl;
    return false;
  }

  for (const char* name : {
    "Ifelix_server",
    "Ifelix_channel",
    "gbad_begin",
    "gbad_end",
    "Ibad_type"})
  {
    if (!tree->GetBranch(name))
    {
      std::cerr << "InttWordShift: missing branch " << name << std::endl;
      return false;
    }
  }
  Int_t felix_server = -1;
  Int_t felix_channel = -1;

  ULong64_t begin = 0;
  ULong64_t end = 0;
  Int_t bad_type = -1;

  if (tree->SetBranchAddress("Ifelix_server", &felix_server) < 0 ||
    tree->SetBranchAddress("Ifelix_channel", &felix_channel) < 0 ||
    tree->SetBranchAddress("gbad_begin", &begin) < 0 ||
    tree->SetBranchAddress("gbad_end", &end) < 0 ||
    tree->SetBranchAddress("Ibad_type", &bad_type) < 0)
  {
    std::cerr << "InttWordShift: branch setup failed" << std::endl;
    return false;
  }

  Config loaded;
  size_t nperiods = 0;

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0)
    {
      return false;
    }

    // Accept CDB server indices (0..7) and packet IDs (3001..3008).
    const int packet_id = (felix_server >= 0 && felix_server <= 7)
                              ? felix_server + 3001 : felix_server;
    if (packet_id < 3001 || packet_id > 3008 ||
        felix_channel < 0 || felix_channel > 13)
    {
      std::cerr << "InttWordShift: invalid FELIX server/channel" << std::endl;
      return false;
    }


    const uint64_t first = begin;
    const uint64_t last = end;
    const int type = bad_type;
    if (first > last || (type != -1 && type != 0 && type != 1))
    {
      std::cerr << "InttWordShift: invalid period" << std::endl;
      return false;
    }

    // Unknown classification has no reliable correction direction.
    if (type == -1)
    {
      continue;
    }

    // CDB: 0 = enhancement, 1 = dropping.
    // Real-data alignment: enhancement -> next chip/channel, dropping -> next ADC/timing.
    const int direction = (type == 0) ? -1 : +1;

    // CDB intervals use the full GL1 clock; raw INTT headers carry 40 bits.
    constexpr uint64_t bco_modulus = uint64_t{1} << 40U;
    constexpr uint64_t bco_mask = bco_modulus - 1;
    if (last - first >= bco_modulus)
    {
      std::cerr << "InttWordShift: period spans more than one 40-bit cycle"
                << std::endl;
      return false;
    }
    const uint64_t first_raw = first & bco_mask;
    const uint64_t last_raw = last & bco_mask;
    auto& periods = loaded[packet_id][felix_channel];
    if (first_raw <= last_raw)
    {
      periods.push_back({first_raw, last_raw, direction});
    }
    else
    {
      // An interval crossing a raw-counter rollover is two inclusive ranges.
      periods.push_back({first_raw, bco_mask, direction});
      periods.push_back({0, last_raw, direction});
    }
    ++nperiods;
    
  }

  for (auto& [packet, fees] : loaded)
  {
    for (auto& [channel, periods] : fees)
    {
      std::sort(periods.begin(), periods.end(),[](const Period& a, const Period& b){return a.begin_bco < b.begin_bco;});

      for (size_t i = 1; i < periods.size(); ++i)
      {
        if (periods[i].begin_bco <=periods[i - 1].end_bco)
        {
          std::cerr << "InttWordShift: overlapping periods" << " packet=" << packet << " fee=" << channel << std::endl;
          return false;
        }
      }
    }
  }
  m_config = std::move(loaded);
  std::cout << "InttWordShift: loaded " << nperiods << " periods from " << filename << std::endl;
  return true;
}
int InttWordShift::GetShift(int packet_id, int fee, uint64_t bco) const
{
  auto packet_it = m_config.find(packet_id);
  if (packet_it == m_config.end())
  {
    return 0;
  }

  auto fee_it = packet_it->second.find(fee);
  if (fee_it == packet_it->second.end())
  {
    return 0;
  }

  for (const auto& period : fee_it->second)
  {
    if (period.begin_bco <= bco && bco <= period.end_bco)
    {
      return period.shift;
    }
  }
  return 0;
}
uint64_t InttWordShift::CalcBCO(unsigned int header, unsigned int bco_word)
{
  uint64_t bco = 0;
  bco |= (static_cast<uint64_t>((header >> 16U) & 0xffU) << 32U);
  bco |= (static_cast<uint64_t>(bco_word & 0xffffU) << 16U);
  bco |= (bco_word >> 16U) & 0xffffU;
  return bco;
}

unsigned int InttWordShift::ShiftHit(unsigned int current, unsigned int next, int shift)
{
  constexpr unsigned int ADC_TIMING = 0xE000FFFFU;
  constexpr unsigned int CHIP_CHAN  = 0x1FFF0000U;
  if (shift == +1)
  {
    // ADC + Timingを次のhitから取得
    return (next & ADC_TIMING) |
           (current & CHIP_CHAN);
  }
  if (shift == -1)
  {
    // Chip + Channelを次のhitから取得
    return (current & ADC_TIMING) |
           (next & CHIP_CHAN);
  }
  return current;
}
std::vector<InttWordShift::HitLocation>
InttWordShift::FindHits(const std::vector<unsigned int>& raw)
{
  std::vector<HitLocation> hits;
  size_t i = 0;
  while (i < raw.size())
  {
    if ((raw[i] & 0xff00ffffU) != 0xad00cadeU)
    {
      ++i;
      continue;
    }
    if (raw.size() - i < 3)
    {
      break;
    }
    const uint64_t bco = CalcBCO(raw[i], raw[i + 1]);
    i += 3;
    while (i < raw.size())
    {
      if ((raw[i] & 0xff00ffffU) == 0xad00cadeU || raw[i] == 0xcafeff80U)
      {
        break;
      }
      hits.push_back({i, bco});
      ++i;
    }
  }
  return hits;
}
InttWordShift::Result InttWordShift::Process(int packet_id, int fee, const std::vector<unsigned int>& raw) const
{
  Result result;
  result.corrected = raw;
  int last_footer = 0;
  for (int i = static_cast<int>(raw.size()) - 1; i > 0; --i)
  {
    if (raw[i] == 0xcafeff80U)
    {
      last_footer = i;
      break;
    }
  }
  result.end_here = last_footer;
  auto packet_it = m_config.find(packet_id);
  if (packet_it == m_config.end())
  {
    return result;
  }
  if (packet_it->second.find(fee) == packet_it->second.end())
  {
    return result;
  }
  const auto hits = FindHits(raw);
  result.end_here = 0;
  result.last_hit = -1;
  for (int i = static_cast<int>(raw.size()) - 1; i > 0; --i)
  {
    if (raw[i] != 0xcafeff80U)
    {
      continue;
    }
    for (const auto& hit : hits)
    {
      if (hit.index > static_cast<size_t>(i))
      {
        result.end_here = i;
        result.last_hit = static_cast<int>(hit.index);
        break;
      }
    }
    if (result.last_hit >= 0)
    {
      break;
    }
  }
  if (result.last_hit < 0)
  {
    return result;
  }
  for (size_t k = 0; k + 1 < hits.size(); ++k)
  {
    const size_t current_index = hits[k].index;
    const size_t next_index = hits[k + 1].index;
    if (current_index >= static_cast<size_t>(result.end_here))
    {
      break;
    }
    const int shift = GetShift(packet_id, fee, hits[k].bco);
    if (shift == 0)
    {
      continue;
    }
    result.corrected[current_index] = ShiftHit(raw[current_index], raw[next_index], shift);
  }
  return result;
}