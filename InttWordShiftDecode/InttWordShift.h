#ifndef INTTWORDSHIFT_H
#define INTTWORDSHIFT_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

class InttWordShift
{
 public:
  struct Result
  {
    std::vector<unsigned int> corrected;
    int end_here = 0;  // 最後に処理可能なFooter
    int last_hit = -1; // 先読み用hit（decodeしない）
  };

  InttWordShift() = default;
  ~InttWordShift() = default;

  bool LoadCDB(const std::string& filename);

  Result Process(int packet_id, int fee, const std::vector<unsigned int>& raw) const;

 private:
  struct Period
  {
    uint64_t begin_bco = 0;
    uint64_t end_bco = 0;
    int shift = 0;
  };

  struct HitLocation
  {
    size_t index = 0;
    uint64_t bco = 0;
  };

  // packet_id -> fee -> 補正区間
  using Config = std::map<int, std::map<int, std::vector<Period>>>;

  Config m_config;

  // CDBから該当する補正方向を取得
  int GetShift(int packet_id, int fee, uint64_t bco) const;

  // 既存intt_poolと同じBCO復元
  static uint64_t CalcBCO(unsigned int header, unsigned int bco_word);

  // FEE別rawからhitの位置と所属BCOを抽出
  static std::vector<HitLocation> FindHits(const std::vector<unsigned int>& raw);

  // 指定されたhitを補正する
  static unsigned int ShiftHit(unsigned int current, unsigned int next, int shift);
};

#endif